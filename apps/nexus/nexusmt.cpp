#include <GL/glew.h>
#include "nexusmt.h"
#include <map>
#include <queue>


using namespace nxs;
using namespace vcg;
using namespace std;

void Policy::Visit(Node *node, std::queue<Node *> &qnode) {    
  std::vector<Node *>::iterator n;
  for(n = node->in.begin(); n != node->in.end(); n++) 
    if(!(*n)->visited) 
      Visit(*n, qnode);  
  
  node->visited = true;
  qnode.push(node);
}

bool FrustumPolicy::Expand(unsigned int patch, Nexus::Entry &entry) {
  if(entry.error == 0) return false;
  float dist = Distance(entry.sphere, frustum.ViewPoint());
  /*Point3f line = frustum.viewPoint() - cell->sphere.center;
    float dist = line.Norm() - cell->sphere.radius;  */
  if(dist < 0) return true;
  //  dist = pow(dist, 1.2);
  /*  cerr << "Dist: " << dist << endl;
  cerr << "Entry error: " << entry.error << endl;
  cerr << "error: " << error << endl;
  cerr << "resolution at dist: " << frustum.Resolution(dist) << endl;*/
  return entry.error > error * frustum.Resolution(dist);
}


NexusMt::NexusMt(): vbo(VBO_AUTO), vbo_size(0),
		    policy(NULL), error(4), realtime(true),
		    mode(SMOOTH) {
  policy = new FrustumPolicy();
}

bool NexusMt::Load(const string &filename) {
  if(!Nexus::Load(filename)) return false;
  LoadHistory();

  use_colors = false;
  use_normals = false;
  use_textures = false;
  use_data = false;

  SetComponent(COLOR, true);
  SetComponent(NORMAL, true);
  SetComponent(TEXTURE, true);
  SetComponent(DATA, true);

  return true;
}

bool NexusMt::InitGL() {
  GLenum ret = glewInit();
  if(ret != GLEW_OK) return false;
  if(!GLEW_ARB_vertex_buffer_object)
    vbo = VBO_OFF;
  return true;
}

void NexusMt::Render() {
  Frustumf frustum;
  frustum.GetView();

  vector<unsigned int> cells;
  if(policy) {
    policy->GetView();
    Extract(cells, policy);
  } else {
    ExtractFixed(cells, error);
  }

  glEnableClientState(GL_VERTEX_ARRAY);
  if(use_colors)
    glEnableClientState(GL_COLOR_ARRAY);
  if(use_normals)
    glEnableClientState(GL_NORMAL_ARRAY);
  //TODO textures and data.

  for(unsigned int i = 0; i < cells.size(); i++) {
    unsigned int cell = cells[i];
    Nexus::Entry &entry = index[cell];
    //frustum culling
    //    if(frustum.Outside(entry.sphere.center, entry.sphere.radius))
    //      continue;
    Patch patch = GetPatch(cell);
    glVertexPointer(3, GL_FLOAT, 0, patch.VertBegin());
    if(use_colors)
      glColorPointer(4, GL_UNSIGNED_BYTE, 0, patch.ColorBegin());
    if(use_normals)
      glNormalPointer(GL_SHORT, 8, patch.Norm16Begin());
    switch(mode) {
    case POINTS:
      glDrawArrays(GL_POINTS, 0, patch.nv); break;
    case DEBUG:
      glColor3ub((cell * 27)%255, (cell * 37)%255, (cell * 87)%255);
    case SMOOTH:
      if(signature & NXS_FACES)
	glDrawElements(GL_TRIANGLES, patch.nf * 3, 
		      GL_UNSIGNED_SHORT, patch.FaceBegin());
      else if(signature & NXS_STRIP)
	glDrawElements(GL_TRIANGLE_STRIP, patch.nf, 
		      GL_UNSIGNED_SHORT, patch.FaceBegin());
      break;
    default: 
      cerr << "Unsupported rendering mode sorry\n";
      exit(0);
      break;
    }
  }
  glDisableClientState(GL_VERTEX_ARRAY);
  glDisableClientState(GL_COLOR_ARRAY);
  glDisableClientState(GL_NORMAL_ARRAY);
}

void NexusMt::SetPolicy(Policy *_policy, bool _realtime) {
  policy = _policy;
  realtime = _realtime;
}

void NexusMt::SetPolicy(PolicyKind kind, float _error, bool _realtime) {
  if(policy) delete policy;
  switch(kind) {
  case FRUSTUM:  policy = new FrustumPolicy(error); break;
  case GEOMETRY: policy = NULL; break;
  default:       policy = NULL; break;
  }
  error = _error;
  realtime = _realtime;
}

void NexusMt::SetVbo(Vbo _vbo, unsigned int _vbo_size) {
  vbo = _vbo;
  if(!GLEW_ARB_vertex_buffer_object)
    vbo = VBO_OFF;
  vbo_size = _vbo_size;
}

bool NexusMt::SetMode(Mode _mode) {
  mode = _mode;
}

bool NexusMt::SetComponent(Component c, bool on) {
  if(c == COLOR && (signature & NXS_COLORS)) 
    use_colors = on;
  if(c == NORMAL && (signature & NXS_NORMALS_SHORT)) 
    use_normals = on;
  if(c == TEXTURE && (signature & NXS_TEXTURES_SHORT)) 
    use_textures = on;
  if(c == DATA && (signature & NXS_DATA32)) 
    use_data = on;
  
  components = COLOR * use_colors + NORMAL * use_normals +
               TEXTURE * use_textures + DATA * use_data;
}

bool NexusMt::SetComponents(unsigned int mask) {
  SetComponent(COLOR, mask & COLOR);
  SetComponent(NORMAL, mask & NORMAL);
  SetComponent(TEXTURE, mask & TEXTURE);
  SetComponent(DATA, mask & DATA);
  
  components = mask;
  
  if( ((mask & COLOR) && !(signature & NXS_COLORS)) ||
      ((mask & NORMAL) && !(signature & NXS_NORMALS_SHORT)) ||
      ((mask & TEXTURE) && !(signature & NXS_TEXTURES_SHORT)) ||
      ((mask & DATA) && !(signature & NXS_DATA32)) )
    return false;
  return true;
}

void NexusMt::LoadHistory() {
  //The last update erases everything.
  assert(history[0].erased.size() == 0);
  
  //maps cell -> node containing it
  map<unsigned int, unsigned int> cell_node;   
  nodes.resize(history.size());

  //building fragments and nodes.
  unsigned int current_node = 0;
  vector<Update>::iterator u;
  for(u = history.begin(); u != history.end(); u++) {      
    
    Node &node = nodes[current_node];
    node.error = 0;
    
    //created cells belong to this node, we look also for max error.      
    for(unsigned int i = 0; i < (*u).created.size(); i++) {
      unsigned int cell = (*u).created[i];
      if(index[cell].error > node.error)
	node.error = index[cell].error;
      
      cell_node[cell] = current_node;        
    }
    
    //Every erased cell already belonged to a node.
    //we record for each node its cells.
    map<unsigned int, vector<unsigned int> > node_erased;      
    
    for(unsigned int i = 0; i < (*u).erased.size(); i++) {
      unsigned int cell = (*u).erased[i];
      assert(cell_node.count(cell));
      node_erased[cell_node[cell]].push_back(cell);               
    }      
    
    //for every node with erased cells we build a frag and 
    //put the corresponding cells in it.      
    map<unsigned int, vector<unsigned int> >::iterator e;
    for(e = node_erased.begin(); e != node_erased.end(); e++) {
      //Build a new Frag.
      Frag fr;
      float max_err = -1;
      
      //Fill it with erased cells.
      vector<unsigned int> &cells = (*e).second;
      vector<unsigned int>::iterator k;
      for(k = cells.begin(); k != cells.end(); k++) {
	unsigned int cell = (*k);
	fr.push_back(cell);
	if(index[cell].error > max_err)
	  max_err = index[cell].error;
      }
      
      //Add the new Frag to the node.
      unsigned int floor_node = (*e).first;
      Node &oldnode = nodes[floor_node];
      oldnode.frags.push_back(fr);
      if(node.error < max_err)
	node.error = max_err;
      
      //Update in and out of the nodes.
      node.in.push_back(&oldnode);
      oldnode.out.push_back(&node);
    }
    current_node++;
  }
}

void NexusMt::ClearHistory() {
  nodes.clear();
}

void NexusMt::ExtractFixed(vector<unsigned int>  &selected, float error) {
  std::vector<Node>::iterator n;
  for(n = nodes.begin(); n != nodes.end(); n++)
    (*n).visited = false;
    
  std::queue<Node *> qnodo;
  qnodo.push(&nodes[0]);
  nodes[0].visited = true;
    
  for( ; !qnodo.empty(); qnodo.pop()) {
    Node &node = *qnodo.front();
      
    std::vector<Frag>::iterator fragment;
    std::vector<Node *>::iterator on;
    for(on = node.out.begin(), fragment = node.frags.begin(); 
	on != node.out.end(); ++on, ++fragment) {

      if((*on)->visited) continue;
	
      if(error < (*on)->error) { //need to expand this node.
	qnodo.push(*on);
	(*on)->visited = 1;
      } else {
	vector<unsigned int>::iterator cell;
	for(cell=(*fragment).begin(); cell != (*fragment).end(); ++cell) 
	  selected.push_back(*cell);                   
      }
    }
  }  
}    

void NexusMt::Extract(std::vector<unsigned int> &selected, Policy *policy) {
  std::vector<Node>::iterator n;
  for(n = nodes.begin(); n != nodes.end(); n++)
    (*n).visited = false;
  
  std::queue<Node *> qnodo;
  qnodo.push(&nodes[0]);
  nodes[0].visited = true;
  
  for( ; !qnodo.empty(); qnodo.pop()) {
    Node &node = *qnodo.front();   
    
    std::vector<Frag>::iterator i;
    std::vector<Node *>::iterator on;
    for(i = node.frags.begin(), on = node.out.begin(); 
	i != node.frags.end(); i++, on++) {
      if((*on)->visited) continue;
      Frag &frag = (*i);
      std::vector<unsigned int>::iterator cell;
      for(cell = frag.begin(); cell != frag.end(); cell++) {              
	if(policy->Expand(*cell, index[*cell]))
	  policy->Visit(*on, qnodo);          
      }
    }
  }
  Select(selected);
}

void NexusMt::Select(vector<unsigned int> &selected) {
  selected.clear();
  std::vector<Node>::iterator i;
  for(i = nodes.begin(); i != nodes.end(); i++) {
    Node &node = *i;
    if(!node.visited)       
      continue;                  
    
    std::vector<Node *>::iterator n;
    std::vector<Frag>::iterator f;        
    for(n = node.out.begin(), f = node.frags.begin(); 
	n != node.out.end(); n++, f++) {
      if(!(*n)->visited || (*n)->error == 0) {
	vector<unsigned int>::iterator c;
	Frag &frag = (*f);
	for(c = frag.begin(); c != frag.end(); c++)            
	  selected.push_back(*c);        
      }      
    } 
  }
}
