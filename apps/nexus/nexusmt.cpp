#include <map>
#include <queue>

#include <GL/glew.h>

#include <ptypes/pasync.h>

#include "nexusmt.h"

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
  PatchInfo &entry = (*index)[cell];    
  Sphere3f &sphere = entry.sphere;
  float dist = Distance(sphere, frustum.ViewPoint());
  if(dist < 0) return 1e20f;
  float error = entry.error/frustum.Resolution(dist);
  if(frustum.IsOutside(sphere.Center(), sphere.Radius()))
    error /= 32;  
  return error;
}

bool NexusMt::Expand(TNode &tnode) {
  //expand if node error > target error
  if(extraction_used >= extraction_max) return false;    
  if(draw_used >= draw_max) return false;
  if(disk_used >= disk_max) return false;  
  
  return tnode.error > target_error; 
}

void NexusMt::NodeVisited(Node *node) {  
  //TODO write this a bit more elegant.
  //first we process arcs removed:
  
  for(unsigned int i = 0; i < node->out.size(); i++) {
    assert(!(node->out[i]->visited));
    Frag &frag = node->frags[i];
    for(unsigned int k = 0; k < frag.size(); k++) {
      unsigned int patch = frag[k];
      PServer::Entry &entry = patches.entries[patch];
      unsigned int ram_size = entry.ram_size; 
      extraction_used += ram_size;     
      PatchInfo &info = index[patch];
      if(!frustum.IsOutside(info.sphere.Center(), info.sphere.Radius())) 
          draw_used += ram_size;
      if(!patches.entries[patch].patch)
          disk_used += ram_size;
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
          unsigned int patch = frag[k];
          PServer::Entry &entry = patches.entries[patch];          
          extraction_used -= entry.ram_size; 
          PatchInfo &info = index[patch];
          if(!frustum.IsOutside(info.sphere.Center(), info.sphere.Radius())) 
            draw_used -= entry.ram_size;
          if(!patches.entries[patch].patch) 
            disk_used -= entry.ram_size;
        }
      }
    }
  }  
}

NexusMt::NexusMt(): vbo_mode(VBO_AUTO), 
                    metric(NULL), mode(SMOOTH) {
  metric = new FrustumMetric();
  metric->index = &index;
  target_error = 4.0f;
  extraction_max = 64000000000;
  draw_max =  12000;   //1M triangles (ZSN)
  disk_max = 300;

}

NexusMt::~NexusMt() {
  if(metric)
    delete metric;
  prefetch.signal();
  prefetch.waitfor();
}

bool NexusMt::Load(const string &filename) {
  index_file = fopen((filename + ".nxs").c_str(), "rb+");
  if(!index_file) return false;
  
  unsigned int readed;
  readed = fread(&signature, sizeof(unsigned int), 1, index_file);
  if(!readed) return false;
  readed = fread(&totvert, sizeof(unsigned int), 1, index_file);
  if(!readed) return false;
  readed = fread(&totface, sizeof(unsigned int), 1, index_file);
  if(!readed) return false;
  readed = fread(&sphere, sizeof(Sphere3f), 1, index_file);
  if(!readed) return false;
  readed = fread(&chunk_size, sizeof(unsigned int), 1, index_file);
    if(!readed) return false;

  unsigned int size; //size of index
  readed = fread(&size, sizeof(unsigned int), 1, index_file);
  if(!readed) return false;

  index.resize(size);
  readed = fread(&index[0], sizeof(PatchInfo), size, index_file);
  if(readed != size) return false;

  patches.ReadEntries(index_file);
  borders.ReadEntries(index_file);
  
  //history size;
  fread(&size, sizeof(unsigned int), 1, index_file);
  vector<unsigned int> buffer;
  buffer.resize(size);
  fread(&(buffer[0]), sizeof(unsigned int), size, index_file);

  //number of history updates
  size = buffer[0];
  history.resize(size);

  unsigned int pos = 1;
  for(unsigned int i = 0; i < size; i++) {
    unsigned int erased = buffer[pos++];
    unsigned int created = buffer[pos++];
    history[i].erased.resize(erased);
    history[i].created.resize(created);
    for(unsigned int e = 0; e < erased; e++) 
      history[i].erased[e] = buffer[pos++];
    for(unsigned int e = 0; e < created; e++) 
      history[i].created[e] = buffer[pos++];
  }
  
  fclose(index_file);
  index_file = NULL;

  if(!patches.Load(filename + ".nxp", signature, chunk_size, true)) 
    return false;
  
  LoadHistory();
  
  use_colors = false;
  use_normals = false;
  use_textures = false;
  use_data = false;

  SetComponent(COLOR, true);
  SetComponent(NORMAL, true);
  SetComponent(TEXTURE, true);
  SetComponent(DATA, true);
  
  SetPrefetchSize(patches.ram_max/2);
  cerr << "Start!\n";
  
  cerr << "Started\n";
  return true;
}

void NexusMt::Close() {
  patches.Close();
}

bool NexusMt::InitGL(Vbo mode, unsigned int vbosize) {
  GLenum ret = glewInit();
  if(ret != GLEW_OK) return false;
  if(!GLEW_ARB_vertex_buffer_object) {
    cerr << "No vbo available!" << endl;
    vbo_mode = VBO_OFF;
  }
  patches.vbo_max = vbosize / chunk_size;
  if(vbo_mode == VBO_OFF)
    patches.vbo_max = 0;

  prefetch.start();
  return true;
}

void NexusMt::Render() {

  vector<unsigned int> cells;
  metric->GetView();
  
  Extract(cells);
  Draw(cells);
}

void NexusMt::Draw(vector<unsigned int> &cells) {
  tri_total = 0;
  tri_rendered = 0;

  frustum.GetView();

  vector<unsigned int> draw;
  for(unsigned int i = 0; i < cells.size(); i++) {
    PatchInfo &entry = index[cells[i]];    
    tri_total += entry.nface;

    if(!frustum.IsOutside(entry.sphere.Center(), entry.sphere.Radius())) {
      draw.push_back(cells[i]);
      tri_rendered += entry.nface;
    }
  }
  prefetch.init(this, draw, visited);

  glEnableClientState(GL_VERTEX_ARRAY);
  if(use_colors)
    glEnableClientState(GL_COLOR_ARRAY);
  if(use_normals)
    glEnableClientState(GL_NORMAL_ARRAY);
  //TODO textures and data.

  unsigned int count = draw.size();
  while(count > 0 || prefetch.draw.get_count()) {
    if(todraw.size()) {
      QueuePServer::Data *data = todraw.back();
      todraw.pop_back();
      Draw((unsigned int)(data->patch), *data);          
    } else {
    // cerr << "Getting message: " << count << endl;
      pt::message *msg = prefetch.draw.getmessage();
      QueuePServer::Data *data = (QueuePServer::Data *)(msg->param);
      if(msg->id == QueuePServer::FLUSH) {        
        if(data->vbo_element) {
	        glDeleteBuffersARB(1, &(data->vbo_element));
	        glDeleteBuffersARB(1, &(data->vbo_array));
        }
        delete data;
        delete msg;
        continue;
      }

      if(msg->id != QueuePServer::DRAW) {
        cerr << "Unknown message!\n";
        continue;
      }

      unsigned int cell = msg->result;        
      Draw(cell, *data);    
      
      delete msg;
    }
    count--;
  }
  glDisableClientState(GL_VERTEX_ARRAY);
  glDisableClientState(GL_COLOR_ARRAY);
  glDisableClientState(GL_NORMAL_ARRAY);
}

void NexusMt::Draw(unsigned int cell, QueuePServer::Data &data) {

  Patch &patch = *(data.patch);
  char *fstart;
  char *vstart;
  char *cstart;
  char *nstart;
  
  if(vbo_mode != VBO_OFF) {
    if(!data.vbo_element) {
      patches.LoadVbo(data);
    }
    
    glBindBufferARB(GL_ARRAY_BUFFER_ARB, data.vbo_array);
    glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER_ARB, data.vbo_element);
    
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
  case PATCHES:
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

void NexusMt::SetExtractionSize(unsigned int r_size) {    
  extraction_max = r_size/patches.chunk_size;
}

void NexusMt::SetMetric(NexusMt::MetricKind kind) {
  //do nothing at the moment.
 
}

void NexusMt::SetError(float error) {
   target_error = error;
}

void NexusMt::SetVboSize(unsigned int _vbo_max) {
  patches.vbo_max = _vbo_max;
}

void NexusMt::SetPrefetchSize(unsigned int size) {
  //TODO do something reasonable with this.
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
  extraction_used = 0;
  draw_used = 0;
  disk_used = 0;
  visited.clear();

  float bestcost = -1;
  int best = -1;
  float cost = 0;

  std::vector<Node>::iterator n;
  for(n = nodes.begin(); n != nodes.end(); n++) 
    (*n).visited = false;    
    
  std::vector<TNode> heap;
  Node *root = &nodes[0];  
  VisitNode(root, heap);

  sequence.clear();
  //vector<Node *> sequence;
  while(heap.size()) {
    pop_heap(heap.begin(), heap.end());
    TNode tnode = heap.back();
    heap.pop_back();

    Node *node = tnode.node;
    if(node->visited) continue;

    bool expand = Expand(tnode);

    if(expand) //{
      VisitNode(node, heap);   
          
    /*  float cost = disk_used * 100;
      if(draw_used > draw_max)
        cost += 1 * (draw_used - draw_max);      
      if(tnode.error > target_error)
        cost += 1000 * (tnode.error - target_error) * (tnode.error - target_error);
      if(best == -1 || cost < bestcost) {
        bestcost = cost;
        best = sequence.size();
      }
    } */
  }
  
  //cerr << "best: " << best << " tot: " << sequence.size() << endl;
  /*for(n = nodes.begin(); n != nodes.end(); n++) 
    (*n).visited = false;    

  for(unsigned int i = 0; i < best; i++) {
    Node *node = sequence[i];
    node->visited = true;
    for(unsigned int n = 0; n < node->out.size(); n++) {
      Node *outnode = node->out[n];
      for(unsigned int k = 0; k < outnode->frags.size(); k++) {      
        for(unsigned int j = 0; j < outnode->frags[k].size(); j++) {
          unsigned int patch = outnode->frags[k][j];
          float error = metric->GetError(patch);      
          visited.push_back(PServer::Item(patch, fabs(error - target_error)));      
        }
      }
    }    
  } */

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
}  */

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
/*bool NexusMt::TestCurrent(Node *node) {
  if(!visited) return;
  for(unsigned int i = 0; i < node->out.size(); i++)
    if(!node->out[k].visited) {
      node->current = true;
      return true;
    }
  node->current = false;
  return false;
} */

void NexusMt::VisitNode(Node *node, vector<TNode> &heap) {
  //TestCurrent(*i);  

  if(node->visited) return;  
  node->visited = true;


  vector<Node *>::iterator i;
  for(i = node->in.begin(); i != node->in.end(); i++) {    
    VisitNode(*i, heap);    
  }
  
  for(unsigned int k = 0; k < node->out.size(); k++) {
    Node *outnode = node->out[k];
    float max_error = metric->GetError(outnode);
  /*  for(unsigned int j = 0; j < node->frags[k].size(); j++) {
      unsigned int patch = node->frags[k][j];
      PServer::Entry &entry = patches.entries[patch];
      if(!entry.patch) 
        max_error /= 2;
    }*/

    for(unsigned int j = 0; j < node->frags[k].size(); j++) {
      unsigned int patch = node->frags[k][j];
      float error = metric->GetError(patch);
      if(max_error < error) max_error = error;
      visited.push_back(PServer::Item(patch, fabs(error - target_error)));
      //      push_heap(visited.begin(), visited.end());
    } 
    
    heap.push_back(TNode(outnode, max_error));
    push_heap(heap.begin(), heap.end());    
  }
  
  sequence.push_back(node);  
  NodeVisited(node);
}

/*void NexusMt::UnvisitNode(Node *node, vector<TNode> &heap) {
  node->current = false;
  vector<Node *>::iterator i;
  for(i = node->in.begin(); i != node->in.end(); i++) 
    if(TestCurrent(*i)) {
      float error = metric->GetError(*i);
      heap.push_back(TNode(*i, error));
      push_heap(heap.begin(), heap.end());      
    }
} */