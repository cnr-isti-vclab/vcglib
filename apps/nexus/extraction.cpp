#include "extraction.h"
#include "metric.h"
#include "nexusmt.h"

#include <algorithm>

using namespace std;
using namespace nxs;

Extraction::Extraction(): target_error(4.0f), extr_max(0xffffffff),
			  draw_max(640), disk_max(100) {
  metric = new FrustumMetric;
}

Extraction::~Extraction() {
  if(metric) delete metric;
}

void Extraction::Extract(NexusMt *_mt) {
  mt = _mt;

  //clear statistics
  extr_used = draw_used = disk_used = 0;
  
  //first we clear the visited flags
  visited.clear();
  visited.resize(mt->history.n_nodes(), false);

  heap.clear();

  History::Node *root = mt->history.Root();
  HeapNode hroot(root, 0);
  Diff(hroot);
  extr_used += hroot.extr;
  draw_used += hroot.draw;
  disk_used += hroot.disk;

  Visit(root); 

  while(heap.size()) {
    pop_heap(heap.begin(), heap.end());
    HeapNode hnode = heap.back();
    heap.pop_back();

    History::Node *node = hnode.node;
    unsigned int id = node - root;
    if(visited[id]) continue;

    if(Expand(hnode)) {
      extr_used += hnode.extr;
      draw_used += hnode.draw;
      disk_used += hnode.disk;
      Visit(node);
    }
  }

  Select();
}

void Extraction::Select() {
  selected.clear(); 
  History::Node *root = mt->history.Root();
  
  History::Node *nodes = mt->history.nodes;
  for(unsigned int i = 0; i < visited.size(); i++) {
    if(!visited[i]) continue;
    History::Node &node = nodes[i];

    History::Node::iterator n;
    for(n = node.out_begin(); n != node.out_end(); n++) {
      unsigned int n_out = (*n).node - root;
      if(!visited[n_out]) {
	History::Link &link = *n;
	for(History::Link::iterator k = link.begin(); k != link.end(); k++) {
	  unsigned int patch = (*k).patch;
	  selected.push_back(patch);
	}
      }
    }
  }
}

void Extraction::Visit(History::Node *node) {
  History::Node *root = mt->history.Root();
  unsigned int n_node = node - root;
  if(visited[n_node]) return;

  visited[n_node] = true;

  History::Node::iterator i;
  for(i = node->in_begin(); i != node->in_end(); i++) {
    unsigned int n_in = (*i).node - root;
    if(visited[n_in]) continue;
    HeapNode hin((*i).node, 0);
    Diff(hin);
    extr_used += hin.extr;
    draw_used += hin.draw;
    disk_used += hin.disk;
    Visit((*i).node);
  }

  for(i = node->out_begin(); i != node->out_end(); i++) {
    float maxerror = 0;
    History::Link &link = *i;
    for(History::Link::iterator k = link.begin(); k != link.end(); k++) {
      Entry &entry = (*mt)[(*k).patch];
      float error =  metric->GetError(entry);
      if(error > maxerror) maxerror = error;
    }
    //TODO this check may be dangerous for non saturating things...
    if(maxerror > target_error) {
      HeapNode hnode((*i).node, maxerror);
      Diff(hnode);
      heap.push_back(hnode);
      push_heap(heap.begin(), heap.end());
    }
  }
}

bool Extraction::Expand(HeapNode &node) {
  if(extr_used >= extr_max) return false; 
  if(draw_used >= draw_max) return false;
  //  if(disk_used >= disk_max) return false;
  return node.error > target_error; 
}

void Extraction::Diff(HeapNode &hnode) {
  History::Node &node = *(hnode.node);
  History::Node::iterator i;
  for(i = node.in_begin(); i != node.in_end(); i++) {
    History::Link &link = *i;
    for(History::Link::iterator k = link.begin(); k != link.end(); k++) {
      unsigned int patch = (*k).patch;
      Entry &entry = (*mt)[patch];
      hnode.extr -= entry.ram_size;
      vcg::Sphere3f &sphere = entry.sphere;
      if(!frustum.IsOutside(sphere.Center(), sphere.Radius()))
	hnode.draw -= entry.ram_size;
      if(!entry.patch)
	hnode.disk -= entry.disk_size;
    }
  }
  
  for(i = node.out_begin(); i != node.out_end(); i++) {
    History::Link &link = *i;
    for(History::Link::iterator k = link.begin(); k != link.end(); k++) {
      unsigned int patch = (*k).patch;
      Entry &entry = (*mt)[patch];
      hnode.extr += entry.ram_size;
      vcg::Sphere3f &sphere = entry.sphere;
      if(!frustum.IsOutside(sphere.Center(), sphere.Radius()))
	hnode.draw += entry.ram_size;
      if(!entry.patch)
	hnode.disk += entry.disk_size;
    }
  }
}
