#include <map>
#include <queue>

#include <GL/glew.h>

#include <ptypes/pasync.h>

#include "nexusmt.h"

using namespace nxs;
using namespace vcg;
using namespace std;

void Stats::Init() {
  ktri = 0;
  kdisk = 0;
  if(count == 25) count = 0;
  if(!count) {
    fps = 25/watch.Time();
    watch.Start();
  }
  count++;
}

/*float ExtractContest::GetError(PatchInfo &entry) {
  Sphere3f &sphere = entry.sphere;
  float dist = Distance(sphere, frustum.ViewPoint());
  if(dist < 0) 
    return 1e20f;
  if(frustum.IsOutside(sphere.Center(), sphere.Radius()))
    return -1;
  return entry.error/frustum.Resolution(dist);
  }*/

/*
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
  if(disk_used >= disk_max) {
    float ratio = disk_max / (float)disk_used;
    float error = tnode.error * ratio;
    return error > target_error;
  }
  return tnode.error > target_error; 
  }*/

/*void NexusMt::NodeVisited(Node *node) {  
  //TODO write this a bit more elegant.
  //first we process arcs removed:
  
  for(unsigned int i = 0; i < node->out.size(); i++) {
    //assert(!(node->out[i]->visited));
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
  }*/

NexusMt::NexusMt() {
  preload.mt = this;
  preload.start();
  
  cerr << "Ramsize: " << ram_max << endl;
  //  SetPrefetchSize(patches.ram_max/2);
  //  prefetch.start();
}

NexusMt::~NexusMt() {
  preload.signal();
  preload.waitfor();
}

bool NexusMt::Load(const string &filename) {
  if(!Nexus::Load(filename, true)) return false;
  if(!history.IsQuick() && !history.UpdatesToQuick())
    return false;
  return true;
}

//void NexusMt::Close() {
  //  prefetch.signal();
  //  prefetch.waitfor();
  //  patches.Close();
//}

bool NexusMt::InitGL(bool vbo) {
  use_vbo = vbo;

  GLenum ret = glewInit();
  if(ret != GLEW_OK) return false;
  if(vbo && !GLEW_ARB_vertex_buffer_object) {
    cerr << "No vbo available!" << endl;
    use_vbo = false;
  }
  return true;
}

void NexusMt::Render(DrawContest contest) {
  Extraction extraction;
  extraction.frustum.GetView();
  extraction.metric->GetView();
  extraction.Extract(this);
  Render(extraction, contest);
}

void NexusMt::Render(Extraction &extraction, DrawContest &contest,
		     Stats *stats) {
  if(stats) stats->Init();

  preload.post(extraction.selected);

  glEnableClientState(GL_VERTEX_ARRAY);
  if((signature & NXS_COLORS) && (contest.attrs & DrawContest::COLOR))
    glEnableClientState(GL_COLOR_ARRAY);
  if((signature & NXS_NORMALS_SHORT) && (contest.attrs & DrawContest::NORMAL))
    glEnableClientState(GL_NORMAL_ARRAY);

  vector<unsigned int> skipped;
  
  vector<unsigned int>::iterator i;
  for(i = extraction.selected.begin(); i != extraction.selected.end(); i++) {
    Entry &entry = operator[](*i);
    vcg::Sphere3f &sphere = entry.sphere;
    if(extraction.frustum.IsOutside(sphere.Center(), sphere.Radius())) 
      continue;

    if(stats) stats->ktri += entry.nface;

    if(!entry.patch) {
      skipped.push_back(*i);
      continue;
    }

    Draw(*i, contest);
  }

  preload.lock.enter();  
  for(i = skipped.begin(); i != skipped.end(); i++) {
    GetPatch(*i);
    Draw(*i, contest);
  }
  //  Flush(false); //not useful now
  preload.lock.leave();

  glDisableClientState(GL_VERTEX_ARRAY);
  glDisableClientState(GL_COLOR_ARRAY);
  glDisableClientState(GL_NORMAL_ARRAY);

  //flushing mem



}



/*void NexusMt::Draw(vector<unsigned int> &cells) {
  stats.start();
  
  todraw.clear();
  vector<QueuePServer::Data> flush;
  for(unsigned int i = 0; i < cells.size(); i++) {
    PatchInfo &entry = index[cells[i]];    
    QueuePServer::Data &data = patches.Lookup(cells[i], entry.nvert,
					      entry.nface, flush);
    todraw.push_back(&data);
  }

  glEnableClientState(GL_VERTEX_ARRAY);
  if(draw_cnt.use_colors)
    glEnableClientState(GL_COLOR_ARRAY);
  if(draw_cnt.use_normals)
    glEnableClientState(GL_NORMAL_ARRAY);
  //TODO textures and data.

  
  //  unsigned int count = draw.size();
  unsigned int count = todraw.size();
  while(count > 0 || prefetch.draw.get_count()) {
    if(todraw.size()) {
      QueuePServer::Data *data = todraw.back();
      todraw.pop_back();
      Draw((unsigned int)(data->patch), *data);          
    } else {
    // cerr << "Getting message: " << count << endl;
      pt::message *msg = prefetch.draw.getmessage();
      QueuePServer::Data *data = (QueuePServer::Data *)(msg->param);
      if(msg->id == Prefetch::FLUSH) {        
	patches.FlushVbo(*data);
	//        if(data->vbo_element) {
	//glDeleteBuffersARB(1, &(data->vbo_element));
	//glDeleteBuffersARB(1, &(data->vbo_array));
        delete data;
        delete msg;
        continue;
      }

      if(msg->id != Prefetch::DRAW) {
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
  }*/

void NexusMt::Draw(unsigned int cell, DrawContest &contest) {
  Entry &entry = operator[](cell);
  Patch &patch = *(entry.patch);
  char *fstart;
  char *vstart;
  char *cstart;
  char *nstart;
  
  if(use_vbo) {
    if(!entry.vbo_element) 
      LoadVbo(entry);
    
    glBindBufferARB(GL_ARRAY_BUFFER_ARB, entry.vbo_array);
    glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER_ARB, entry.vbo_element);
    
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
  if(contest.attrs & DrawContest::COLOR)
    glColorPointer(4, GL_UNSIGNED_BYTE, 0, cstart);
  if(contest.attrs & DrawContest::NORMAL)
    glNormalPointer(GL_SHORT, 8, nstart);
  
  switch(contest.mode) {
  case DrawContest::POINTS:
    glDrawArrays(GL_POINTS, 0, patch.nv); break;
  case DrawContest::PATCHES:
    glColor3ub((cell * 27)%255, (cell * 37)%255, (cell * 87)%255);
  case DrawContest::SMOOTH:
    if(signature & NXS_FACES)
      glDrawElements(GL_TRIANGLES, patch.nf * 3, 
		     GL_UNSIGNED_SHORT, fstart);
    else if(signature & NXS_STRIP)
      glDrawElements(GL_TRIANGLE_STRIP, patch.nf, 
		     GL_UNSIGNED_SHORT, fstart);
    break;
  case DrawContest::FLAT:
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

void NexusMt::FlushPatch(unsigned int id) {
  Entry &entry = operator[](id);
  if(entry.vbo_element)
    FlushVbo(entry);

  if(entry.patch->start)
    delete [](entry.patch->start);
  delete entry.patch;  
  entry.patch = NULL;    
  ram_used -= entry.ram_size;      
}

void NexusMt::LoadVbo(Entry &entry) { 
  if(entry.vbo_element) return;
  glGenBuffersARB(1, &entry.vbo_element);
  assert(entry.vbo_element);
  glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER_ARB, entry.vbo_element);

  Patch &patch  = *entry.patch;    
  unsigned int size = patch.nf * sizeof(unsigned short);
  if((signature & NXS_FACES) != 0) size *= 3;
  
  glBufferDataARB(GL_ELEMENT_ARRAY_BUFFER_ARB, size, patch.FaceBegin(),
		  GL_STATIC_DRAW_ARB);
  vbo_used += size;
  
  //TODO fix this when we allow data :p
  size = sizeof(float) * patch.dstart;
    
  glGenBuffersARB(1, &entry.vbo_array);
  assert(entry.vbo_array);
  glBindBufferARB(GL_ARRAY_BUFFER_ARB, entry.vbo_array);
    
  glBufferDataARB(GL_ARRAY_BUFFER_ARB, size, patch.VertBegin(), 
		  GL_STATIC_DRAW_ARB);
    
  vbo_used += size;
  delete [](entry.patch->start);
  entry.patch->start = NULL;
}

void NexusMt::FlushVbo(Entry &entry) {
  if(!entry.vbo_element) return;
  glDeleteBuffersARB(1, &entry.vbo_element);
  glDeleteBuffersARB(1, &entry.vbo_array);
  entry.vbo_element = 0;
  entry.vbo_array = 0;

  Patch &patch  = *entry.patch; 
  vbo_used -= patch.nf * sizeof(unsigned short);
  vbo_used -= sizeof(float) * patch.dstart;
}

/*void NexusMt::SetExtractionSize(unsigned int r_size) {    
  extract_cnt.extraction_max = r_size/patches.chunk_size;
}

//void NexusMt::SetMetric(NexusMt::MetricKind kind) {
  //do nothing at the moment.
 
//}

void NexusMt::SetError(float error) {
   extract_cnt.target_error = error;
}

//void NexusMt::SetVboSize(unsigned int _vbo_max) {
//WORKING  patches.vbo_max = _vbo_max;
//}

void NexusMt::SetPrefetchSize(unsigned int size) {
  //TODO do something reasonable with this.
}

bool NexusMt::SetMode(DrawContest::Mode _mode) {
  draw_cnt.mode = _mode;
  return true;
}

bool NexusMt::SetComponent(DrawContest::Component c, bool on) {
  if(c == DrawContest::COLOR && (signature & NXS_COLORS)) 
    draw_cnt.use_colors = on;
  if(c == DrawContest::NORMAL && (signature & NXS_NORMALS_SHORT)) 
    draw_cnt.use_normals = on;
  if(c == DrawContest::TEXTURE && (signature & NXS_TEXTURES_SHORT)) 
    draw_cnt.use_textures = on;
  if(c == DrawContest::DATA && (signature & NXS_DATA32)) 
    draw_cnt.use_data = on;
  
  //  unsigned int components = DrawContest::COLOR * use_colors + 
  //    DrawContest::NORMAL * use_normals +
  //    DrawContest::TEXTURE * use_textures + 
  //    DrawContest::DATA * use_data;
  return true;
}

bool NexusMt::SetComponents(unsigned int mask) {
  SetComponent(DrawContest::COLOR, (mask & DrawContest::COLOR) != 0);
  SetComponent(DrawContest::NORMAL, (mask & DrawContest::NORMAL) != 0);
  SetComponent(DrawContest::TEXTURE, (mask & DrawContest::TEXTURE) != 0);
  SetComponent(DrawContest::DATA, (mask & DrawContest::DATA) != 0);
  
  //  unsigned int components = mask;
  
  if( ((mask & DrawContest::COLOR) && !(signature & NXS_COLORS)) ||
      ((mask & DrawContest::NORMAL) && !(signature & NXS_NORMALS_SHORT)) ||
      ((mask & DrawContest::TEXTURE) && !(signature & NXS_TEXTURES_SHORT)) ||
      ((mask & DrawContest::DATA) && !(signature & NXS_DATA32)) )
    return false;
  return true;
  }*/

//TODO: nodes and fragment sholuld be kept in another way,
// so we can save that structure instead of the history!



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
	      for(cell=(*fragment).begin(); cell != (*fragment).end(); ++cell) 	        selected.push_back(*cell);                   
      }
    }
  }  
} */   

/*void NexusMt::Expand(Node &node, vector<HeapNode> &expand) {
  if(node.visited) return;
  cerr << "Expanding\n";
  for(unsigned int i = 0; i < node.in_size; i++) {
    Link &link = innodes[node.in + i];
    Node &parent = nodes[link.node];
    if(!parent.visited)
      Expand(parent, expand);
  }

  for(unsigned int o = 0; o < node.out_size; o++) {
    Link &link = outnodes[node.out + o];
    Node &outnode = nodes[link.node];
    if(outnode.visited) continue;
    assert(link.error != -1);
    if(link.error > extract_cnt.target_error) {
      HeapNode hnode(&outnode, node.error - extract_cnt.target_error);
      expand.push_back(hnode);
      push_heap(expand.begin(), expand.end());
    }
  }

  unsigned int d_extr = 0;
  unsigned int d_draw = 0;
  unsigned int d_disk = 0;
  Diff(node, d_extr, d_draw, d_disk);
  extract_cnt.extraction_used += d_extr;
  extract_cnt.draw_used += d_draw;
  extract_cnt.disk_used += d_disk;
  node.visited = true;
}

void NexusMt::Contract(Node &node) {
  if(!node.visited) return;
  cerr << "Contracting...\n";
  for(unsigned int i = 0; i < node.out_size; i++) {
    Link &link = outnodes[node.out + i];
    Node &child = nodes[link.node];
    if(child.visited)
      Contract(child);
  }

  unsigned int d_extr = 0;
  unsigned int d_draw = 0;
  unsigned int d_disk = 0;
  Diff(node, d_extr, d_draw, d_disk);
  extract_cnt.extraction_used -= d_extr;
  extract_cnt.draw_used -= d_draw;
  extract_cnt.disk_used -= d_disk;
  node.visited = false;
}

void NexusMt::GetErrorLink(Node &node, Link &link) {
  if(link.error != -1) return;
  for(unsigned int c = link.frag; frags[c].patch != 0xffffffff; c++) {
    Cell &cell = frags[c];
    if(cell.error != -1) continue;
    PatchInfo &info = index[cell.patch];
    cell.error = extract_cnt.GetError(info);
    cerr << "Cell.error: " << cell.error << endl;
    if(link.error < cell.error) link.error = cell.error;
  }
}

void NexusMt::GetErrorNode(Node &node) {
  if(node.error != -1) return;
  for(unsigned int o = 0; o < node.out_size; o++) {
    Link &link = outnodes[node.out + o];
    GetErrorLink(node, link);
    if(node.error < link.error) node.error = link.error;
  }
}

void NexusMt::NodeExplore(Node &node, 
			  vector<HeapNode> &expand,
			  vector<HeapNode> &contract) {
  assert(node.visited);
  GetErrorNode(node);
  
  cerr << "node.error: " << node.error <<
    " target: " << extract_cnt.target_error << endl;
  if(node.error > extract_cnt.target_error) {
    cerr << "Expandible node...\n";
    for(unsigned int o = 0; o < node.out_size; o++) {
      Link &link = outnodes[node.out + o];
      assert(link.error != -1);
      Node &outnode = nodes[link.node];
      if(outnode.visited) continue;

      if(link.error > extract_cnt.target_error) {
	HeapNode hnode(&outnode, node.error - extract_cnt.target_error);
	expand.push_back(hnode);
	push_heap(expand.begin(), expand.end());
      }
    }
  }
  
  float max_err = -1;
  
  bool can_contract = true;
  for(unsigned int o = 0; o < node.out_size; o++) {
    Link &link = outnodes[node.out + o];
    Node &outnode = nodes[link.node];
    if(node.visited) {
      can_contract = false;
      break;
    }
  }
  if(can_contract) {
    for(unsigned int o = 0; o < node.in_size; o++) {
      Link &link = innodes[node.in + o];
      GetErrorLink(node, link);
      if(max_err < link.error) max_err = link.error;
    }
    HeapNode hnode(&node, extract_cnt.target_error - max_err);
    contract.push_back(hnode);
    push_heap(contract.begin(), contract.end());
  }
}

void NexusMt::FrontExplore(vector<HeapNode> &expand,
			   vector<HeapNode> &contract) {


  extract_cnt.extraction_used = 0;
  extract_cnt.draw_used = 0;
  extract_cnt.disk_used = 0;
    
  std::vector<Node>::iterator i;
  for(i = nodes.begin(); i != nodes.end(); i++) {
    Node &node = *i;
    if(!node.visited)       
      continue;                  

    for(int k = 0; k < node.out_size; k++) {
      Link &link = outnodes[node.out + k];
      assert(link.node < nodes.size());
      Node &children = nodes[link.node];
      if(!children.visited) {
	cerr << "Exploting child...\n";
	NodeExplore(node, expand, contract);

	for(unsigned int c = link.frag; frags[c].patch != 0xffffffff; c++) {
	  PServer::Entry &entry = patches.entries[c];
	  if(!entry.patch) extract_cnt.disk_used += entry.disk_size;
	  extract_cnt.extraction_used += entry.ram_size;
	  vcg::Sphere3f &sphere = index[c].sphere;
	  if(!extract_cnt.frustum.IsOutside(sphere.Center(), sphere.Radius()))
	    extract_cnt.draw_used += entry.ram_size;
	}
      }
    }
  }
}
 
void NexusMt::Diff(Node &node, 
		   unsigned int &d_extr, 
		   unsigned int &d_draw,
		   unsigned int &d_disk) {
  d_extr = 0;
  d_draw = 0;
  d_disk = 0;
  for(unsigned int o = 0; o < node.out_size; o++) {
    Link &link = outnodes[node.out + o];
    for(unsigned int c = link.frag; frags[c].patch != 0xffffffff; c++) {
      PServer::Entry &entry = patches.entries[c];
      if(!entry.patch) d_disk += entry.disk_size;
      d_extr += entry.ram_size;
      Sphere3f &sphere = index[c].sphere;
      if(!extract_cnt.frustum.IsOutside(sphere.Center(), sphere.Radius()))
	d_draw += entry.ram_size;
    }
  }
  for(unsigned int o = 0; o < node.in_size; o++) {
    Link &link = innodes[node.in + o];
    for(unsigned int c = link.frag; frags[c].patch != 0xffffffff; c++) {
      PServer::Entry &entry = patches.entries[c];
      if(!entry.patch) d_disk -= entry.disk_size;
      d_extr -= entry.ram_size;
      if(!extract_cnt.frustum.IsOutside(sphere.Center(), sphere.Radius()))
	d_draw -= entry.ram_size;
    }
  }
  }*/

/*void NexusMt::Extract(std::vector<unsigned int> &selected) {
  extract_cnt.frustum.GetView();


  //clear error (no more valid...)
  for(unsigned int i = 0; i < nodes.size(); i++) 
    nodes[i].error = -1;
  for(unsigned int i = 0; i < outnodes.size(); i++) 
    outnodes[i].error = -1;
  for(unsigned int i = 0; i < innodes.size(); i++) 
    innodes[i].error = -1;
  for(unsigned int i = 0; i < frags.size(); i++)
    frags[i].error = -1;


  std::vector<HeapNode> expand;
  std::vector<HeapNode> contract;

  if(!nodes[0].visited)
    nodes[0].visited = true;
  //explore front
  FrontExplore(expand, contract);

  bool can_expand = true;
  bool can_contract = true;
  float current_error = 1e30;
  while(can_expand || can_contract) {
    if(can_expand) {
      if(!expand.size()) {
	can_expand = false;
	continue;
      }
      cerr << "trying to expand..\n";
      pop_heap(expand.begin(), expand.end());
      HeapNode hnode = expand.back();
      
      Node *node = hnode.node;
      //      if(node->visited) {
	//cerr << "Failed cause already visited\n";
//	expand.pop_back();
	//continue; //already expanded...
//	}
      unsigned int d_extr = 0;
      unsigned int d_draw = 0;
      unsigned int d_disk = 0;
      Diff(*node, d_extr, d_draw, d_disk);
      
      if(extract_cnt.disk_used + d_disk > extract_cnt.disk_max) {
	cerr << "Failed cause disk...\n";
	cerr << "d_disl: " << d_disk << " extract_cnt.disk_max "
	     << extract_cnt.disk_max << endl;
	expand.pop_back();
	continue;
      }
      if(extract_cnt.extraction_used + d_extr > extract_cnt.extraction_max ||
	 extract_cnt.draw_used + d_draw > extract_cnt.draw_max) {
	cerr << "Failed...\n";
	can_expand = false;
	push_heap(expand.begin(), expand.end());
	continue;
      } else {
	cerr << "Expanding...\n";
	current_error = hnode.error;
	expand.pop_back();
	Expand(*node, expand);
      }
    } else if(can_contract) {
      if(!contract.size()) {
	can_contract = false;
	continue;
      }
      pop_heap(contract.begin(), contract.end());
      HeapNode hnode = contract.back();

      if(hnode.error < extract_cnt.target_error - current_error) {
	can_contract = false;
	continue;
      }
      Node *node = hnode.node;
      if(!node->visited) {
	contract.pop_back();
	continue;//already contracted
      }
      
      unsigned int d_extr = 0;
      unsigned int d_draw = 0;
      unsigned int d_disk = 0;
      Diff(*node, d_extr, d_draw, d_disk);

      if(extract_cnt.disk_used - d_disk > extract_cnt.disk_max) {
	expand.pop_back();
	continue;
      }
      if(extract_cnt.extraction_used - d_extr < extract_cnt.extraction_max &&
	 extract_cnt.draw_used - d_draw < extract_cnt.draw_max) {
	can_expand = true;
      }
      current_error = hnode.error;
      contract.pop_back();
      Contract(*node);
    }
  }
  Select(selected);   
}*/

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

/*void NexusMt::Select(vector<unsigned int> &selected) {
  selected.clear();
  std::vector<Node>::iterator i;
  for(i = nodes.begin(); i != nodes.end(); i++) {
    Node &node = *i;
    if(!node.visited)       
      continue;                  
    
    for(unsigned int o = 0; o < node.out_size; o++) {
      Link &link = outnodes[node.out + o];
      for(unsigned int c = link.frag; frags[c].patch != 0xffffffff; c++) {
	Cell &cell = frags[c];
	selected.push_back(cell.patch);
      }
    }
    }*/
  //  for(unsigned int i = 0; i < selected.size(); i++) {
  //    cerr << "Sel: " << selected[i] << endl;
  //  }
    
    /*    std::vector<Node *>::iterator n;
    std::vector<>::iterator f;        
    for(n = node.out.begin(), f = node.frags.begin(); 
	    n != node.out.end(); n++, f++) {
      if(!(*n)->visited || (*n)->error == 0) {
	      vector<unsigned int>::iterator c;
	      Frag &frag = (*f);
	      for(c = frag.begin(); c != frag.end(); c++)            
	        selected.push_back(*c);        
      }      
    } 
    }*/
//}
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

/*void NexusMt::VisitNode(Node *node, vector<TNode> &heap) {
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
  }*/

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


/*void NexusMt::LoadHistory() {
  //The last update erases everything.
  assert(history[0].erased.size() == 0);
  
  //maps cell -> node containing it
  map<unsigned int, unsigned int> cell_node;   
  //maps node -> Links
  map<unsigned int, vector<Link> > node_inlinks;
  map<unsigned int, vector<Link> > node_outlinks;
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

      //node.in.push_back(innodes.size());

      unsigned int floor_node = (*e).first;
      Node &oldnode = nodes[floor_node];

      Link inlink;
      inlink.node = floor_node;
      inlink.frag = frags.size();

      Link outlink;
      outlink.node = current_node;
      outlink.frag = frags.size();

      float max_err = -1;
      
      //Fill it with erased cells.
      vector<unsigned int> &cells = (*e).second;
      vector<unsigned int>::iterator k;
      for(k = cells.begin(); k != cells.end(); k++) {
	Cell cell;
	cell.patch = (*k);
        assert(cell.patch < index.size());
	frags.push_back(cell);
	if(index[cell.patch].error > max_err)
	  max_err = index[cell.patch].error;
      }         
      Cell cell;
      cell.patch = 0xffffffff;
      cell.error = -1;
      frags.push_back(cell);

      inlink.error = max_err;
      outlink.error = max_err;

      //Add the new Frag to the node.

      node_outlinks[floor_node].push_back(outlink);
      node_outlinks[current_node].push_back(inlink);
      if(node.error < max_err)
	      node.error = max_err;
      
      //Update in and out of the nodes.
    }
    current_node++;
  }

  map<unsigned int, vector<Link> >::iterator k;
  for(k = node_outlinks.begin(); k != node_outlinks.end(); k++) {
    unsigned int inode = (*k).first;
    vector<Link> &links = (*k).second;
    nodes[inode].out = outnodes.size();
    nodes[inode].out_size = links.size();

    for(unsigned int i = 0; i < links.size(); i++) 
      outnodes.push_back(links[i]);
  }

  for(k = node_inlinks.begin(); k != node_inlinks.end(); k++) {
    unsigned int inode = (*k).first;
    vector<Link> &links = (*k).second;
    nodes[inode].in = innodes.size();
    nodes[inode].in_size = links.size();

    for(unsigned int i = 0; i < links.size(); i++) 
      innodes.push_back(links[i]);
  }
}

void NexusMt::ClearHistory() {
  nodes.clear();
  innodes.clear();
  outnodes.clear();
  frags.clear();
  }*/
