#include <GL/glew.h>
#include "nexusmt.h"
#include <map>
#include <queue>


using namespace nxs;
using namespace vcg;
using namespace std;

/*void Policy::Visit(Node *node, std::queue<Node *> &qnode) {    
  std::vector<Node *>::iterator n;
  for(n = node->in.begin(); n != node->in.end(); n++) 
    if(!(*n)->visited) 
      Visit(*n, qnode);  
  
  node->visited = true;
  qnode.push(node);
}

bool FrustumPolicy::Expand(unsigned int patch, Nexus::PatchInfo &entry) {
  if(entry.error == 0) return false;
  float dist = Distance(entry.sphere, frustum.ViewPoint());
  if(dist < 0) return true;
  //  dist = pow(dist, 1.2);
  return entry.error > error * frustum.Resolution(dist);
} */

float Metric::GetError(Node *node) {
  float max_error = 0;
  vector<Frag>::iterator frag;
  for(frag = node->frags.begin(); frag != node->frags.end(); frag++) {
    float error = GetError(*frag);
    if(max_error < error) max_error = error;
  }
  return max_error;
}

float Metric::GetError(Frag &frag) {
  float max_error = 0;
  vector<unsigned int>::iterator cell;
  for(cell = frag.begin(); cell != frag.end(); cell++) {              
    float error = GetError(*cell);
    if(max_error < error) max_error = error;
  }
  return max_error;
}

float FrustumMetric::GetError(unsigned int cell) {
  float max_error = 0;
  Nexus::PatchInfo &entry = (*index)[cell];    
  Sphere3f &sphere = entry.sphere;
  float dist = Distance(sphere, frustum.ViewPoint());
  if(dist < 0) return 1e40;
  float error = entry.error/frustum.Resolution(dist);
  if(frustum.IsOutside(sphere.Center(), sphere.Radius()))
    error /= 4;
  return error;
}

void Policy::Init() {
  ram_used = 0;
}

bool Policy::Expand(TNode &node) {
  //expand if node error > target error
  if(ram_used >= ram_size) return false;
  //cerr << "Error: " << error << " node.error: " << node.error << endl;
  return node.error > error; 
}

void Policy::NodeVisited(Node *node) {  
  //TODO write this a bit more elegant.
  //first we process arcs removed:
  
  for(unsigned int i = 0; i < node->out.size(); i++) {
    assert(!(node->out[i]->visited));
    Frag &frag = node->frags[i];
    for(unsigned int k = 0; k < frag.size(); k++) {
      PatchEntry &entry = (*entries)[frag[k]];
      ram_used += entry.ram_size;      
    }
  }

  vector<Node *>::iterator from;
  for(from = node->in.begin(); from != node->in.end(); from++) {
    assert((*from)->visited);
    vector<Frag> &frags = (*from)->frags;
    for(unsigned int i = 0; i < frags.size(); i++) {
      if((*from)->out[i] == node) {
        vector<unsigned int> &frag = frags[i];
        for(unsigned int k = 0; k < frag.size(); k++) {
          PatchEntry &entry = (*entries)[frag[k]];          
          ram_used -= entry.ram_size;          
        }
      }
    }
  }  
}

Prefetch::Prefetch(): thread(true), nexus(NULL) {}

Prefetch::~Prefetch() { signal(); waitfor(); }

void Prefetch::execute() {
  assert(nexus != NULL);
  while(1) {
    if(get_signaled()) return;
    if(cells.size() == 0 || 
       nexus->patches.ram_used > nexus->patches.ram_size) {
      relax(100);
      continue;
    }
    cells_mx.lock();

    pop_heap(cells.begin(), cells.end());
    unsigned int cell = cells.back();
    cells.pop_back();

    patch_mx.lock();
    Patch &patch = nexus->GetPatch(cell, false);
    patch_mx.unlock();
    
    cells_mx.unlock();
  }
}

NexusMt::NexusMt(): vbo_mode(VBO_AUTO), 
                    metric(NULL), mode(SMOOTH), prefetching(true) {    
  metric = new FrustumMetric();
  metric->index = &index;
  policy.error = 4;
  policy.ram_size = 64000000;

  prefetch.nexus = this;
}

NexusMt::~NexusMt() {}

bool NexusMt::Load(const string &filename, bool readonly) {
  if(!Nexus::Load(filename, readonly)) return false;
  LoadHistory();

  policy.entries = &patches.patches;    

  use_colors = false;
  use_normals = false;
  use_textures = false;
  use_data = false;

  SetComponent(COLOR, true);
  SetComponent(NORMAL, true);
  SetComponent(TEXTURE, true);
  SetComponent(DATA, true);
  
  SetPrefetching(prefetching);
  return true;
}

bool NexusMt::InitGL(Vbo mode, unsigned int vbosize) {
  GLenum ret = glewInit();
  if(ret != GLEW_OK) return false;
  if(!GLEW_ARB_vertex_buffer_object) {
    cerr << "No vbo available!" << endl;
    vbo_mode = VBO_OFF;
  }
  patches.vbo_size = vbosize / patches.chunk_size;
  if(vbo_mode == VBO_OFF) 
    patches.vbo_size = 0;
  return true;
}

void NexusMt::Render() {
  patches.Flush();

  vector<unsigned int> cells;
  metric->GetView();
  policy.Init();
  tri_total = 0;
  tri_rendered = 0;

  Extract(cells);
  Draw(cells);
}

void NexusMt::Draw(vector<unsigned int> &cells) {
  Frustumf frustum;
  frustum.GetView();

  glEnableClientState(GL_VERTEX_ARRAY);
  if(use_colors)
    glEnableClientState(GL_COLOR_ARRAY);
  if(use_normals)
    glEnableClientState(GL_NORMAL_ARRAY);
  //TODO textures and data.

  for(unsigned int i = 0; i < cells.size(); i++) {    
    unsigned int cell = cells[i];
    Nexus::PatchInfo &entry = index[cell];    
    tri_total += entry.nface;
    //frustum culling
    if(frustum.IsOutside(entry.sphere.Center(), entry.sphere.Radius()))
      continue;

    tri_rendered += entry.nface;
    
    Patch &patch = GetPatch(cell, false);
    char *fstart;
    char *vstart;
    char *cstart;
    char *nstart;

    if(vbo_mode != VBO_OFF) {
      unsigned int vbo_array;
      unsigned int vbo_element;
      patches.GetVbo(cell, vbo_element, vbo_array);
      assert(vbo_element);
      assert(vbo_array);

      glBindBufferARB(GL_ARRAY_BUFFER_ARB, vbo_array);
      glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER_ARB, vbo_element);

      fstart = NULL;
      vstart = NULL;
      cstart = (char *)(sizeof(float) * patch.cstart);
      nstart = (char *)(sizeof(float) * patch.nstart);
    } else {
      fstart = (char *)patch.FaceBegin();
      vstart = (char *)patch.VertBegin();
      cstart = (char *)patch.ColorBegin();
      nstart = (char *)patch.Norm16Begin();
    }

    glVertexPointer(3, GL_FLOAT, 0, vstart);
    if(use_colors)
      glColorPointer(4, GL_UNSIGNED_BYTE, 0, cstart);
    if(use_normals)
      glNormalPointer(GL_SHORT, 8, nstart);

    switch(mode) {
    case POINTS:
      glDrawArrays(GL_POINTS, 0, patch.nv); break;
    case DEBUG:
      glColor3ub((cell * 27)%255, (cell * 37)%255, (cell * 87)%255);
    case SMOOTH:
      if(signature & NXS_FACES)
	glDrawElements(GL_TRIANGLES, patch.nf * 3, 
		       GL_UNSIGNED_SHORT, fstart);
      else if(signature & NXS_STRIP)
	glDrawElements(GL_TRIANGLE_STRIP, patch.nf, 
		       GL_UNSIGNED_SHORT, fstart);
      break;
    case FLAT:
      if(signature & NXS_FACES) {
	glBegin(GL_TRIANGLES);
	for(int i = 0; i < patch.nf; i++) {
	  unsigned short *f = patch.Face(i);
	  Point3f &p0 = patch.Vert(f[0]);
	  Point3f &p1 = patch.Vert(f[1]);
	  Point3f &p2 = patch.Vert(f[2]);
	  Point3f n = ((p1 - p0) ^ (p2 - p0));
	  glNormal3f(n[0], n[1], n[2]);
	  glVertex3f(p0[0], p0[1], p0[2]);
	  glVertex3f(p1[0], p1[1], p1[2]);
	  glVertex3f(p2[0], p2[1], p2[2]);
	}
	glEnd();
      } else if(signature & NXS_STRIP) {
	cerr << "Unsupported rendering mode sorry\n";
	exit(0);
      }
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

void NexusMt::SetRamExtractionSize(unsigned int r_size) {    
  policy.ram_size = r_size/patches.chunk_size;
}

void NexusMt::SetMetric(NexusMt::MetricKind kind) {
  //do nothing at the moment.
 
}

void NexusMt::SetError(float error) {
   policy.error = error;
}

void NexusMt::SetVboSize(unsigned int _vbo_size) {
  patches.vbo_size = _vbo_size;
}

void NexusMt::SetPrefetching(bool on) {
  if(on && prefetch.Running()) return;
  if(!on && !prefetch.Running()) return;
  if(on) {
    prefetch.cells.clear();
    prefetch.start();
  } else {
    prefetch.signal();
    prefetch.waitfor();
  }
  prefetching = on;
}

bool NexusMt::SetMode(Mode _mode) {
  mode = _mode;
  return true;
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
  return true;
}

bool NexusMt::SetComponents(unsigned int mask) {
  SetComponent(COLOR, (mask & COLOR) != 0);
  SetComponent(NORMAL, (mask & NORMAL) != 0);
  SetComponent(TEXTURE, (mask & TEXTURE) != 0);
  SetComponent(DATA, (mask & DATA) != 0);
  
  components = mask;
  
  if( ((mask & COLOR) && !(signature & NXS_COLORS)) ||
      ((mask & NORMAL) && !(signature & NXS_NORMALS_SHORT)) ||
      ((mask & TEXTURE) && !(signature & NXS_TEXTURES_SHORT)) ||
      ((mask & DATA) && !(signature & NXS_DATA32)) )
    return false;
  return true;
}

//TODO: nodes and fragment sholuld be kept in another way,
// so we can save that structure instead of the history!

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
        assert(cell < index.size());
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

/*void NexusMt::ExtractFixed(vector<unsigned int>  &selected, float error) {
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
} */   

void NexusMt::Extract(std::vector<unsigned int> &selected) {
  std::vector<Node>::iterator n;
  for(n = nodes.begin(); n != nodes.end(); n++) {
    (*n).visited = false;
    (*n).pushed = false;
  }
  
  std::vector<TNode> heap;
  Node *root = &nodes[0];  
  VisitNode(root, heap);
  
  while(heap.size()) {
    pop_heap(heap.begin(), heap.end());
    TNode tnode = heap.back();
    heap.pop_back();

    Node *node = tnode.node;
    if(node->visited) continue;

    bool expand = policy.Expand(tnode);

    if(expand)  
      VisitNode(node, heap);                
  }
  Select(selected);   
}

/*void NexusMt::Extract(std::vector<unsigned int> &selected, Policy *policy) {
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
}*/

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

void NexusMt::VisitNode(Node *node, vector<TNode> &heap) {
  if(node->visited) return;
  
  vector<Node *>::iterator i;
  for(i = node->in.begin(); i != node->in.end(); i++)     
    VisitNode(*i, heap);
  
  for(unsigned int k = 0; k < node->out.size(); k++) {
    Node *outnode = node->out[k];
    float error = metric->GetError(node->frags[k]);
    //  if(node->pushed) continue
    heap.push_back(TNode(outnode, error));
    push_heap(heap.begin(), heap.end());    
  }
  
  node->visited = true;
  policy.NodeVisited(node);
}

