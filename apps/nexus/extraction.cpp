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
  root = mt->history.Root();
  sink = root + (mt->history.n_nodes()-1);

  //clear statistics
  extr_used = draw_used = disk_used = 0;
  
  //first we clear the visited flags
  visited.clear();
  visited.resize(mt->history.n_nodes(), false);

  heap.clear();

  Visit(root); 

  while(heap.size()) {
    pop_heap(heap.begin(), heap.end());
    HeapNode hnode = heap.back();
    heap.pop_back();

    Node *node = hnode.node;
    if(visited[node - root]) continue;

    if(Expand(hnode)) 
      Visit(node);
  }
  Select();
  draw_size = selected.size();
}

void Extraction::Init() { 
  front.clear();
  back.clear();
  
  Cost cost;

  Node *nodes = mt->history.nodes;
  for(unsigned int i = 0; i < visited.size(); i++) {
    if(!visited[i]) continue;
    Node &node = nodes[i];

    bool cancoarse = true;
    Node::iterator n;
    for(n = node.out_begin(); n != node.out_end(); n++) {

      if(!visited[(*n).node - root]) {
	float maxerror = 0;

	Link &link = *n;
	for(Link::iterator k = link.begin(); k != link.end(); k++) {
	  unsigned int patch = (*k).patch;
	  Entry &entry = (*mt)[patch];
	  float error = metric->GetError(entry);
	  if(error > maxerror) maxerror = error;

	  cost.extr += entry.ram_size;
	  vcg::Sphere3f &sphere = entry.sphere;
	  if(!frustum.IsOutside(sphere.Center(), sphere.Radius()))
	    cost.draw += entry.ram_size;
	  if(!entry.patch)
	    cost.disk += entry.disk_size;
	}
	if((*n).node != sink && maxerror > target_error)
	  front.push_back(HeapNode((*n).node, maxerror));
      } else
	cancoarse = false;
    }
    if(cancoarse && &node != root) {
      float error = GetRefineError(&node);
      back.push_back(HeapNode(&node, error));
    }
  }
  make_heap(front.begin(), front.end());
  make_heap(back.begin(), back.end(), greater<HeapNode>());

  extr_used = cost.extr;
  draw_used = cost.draw;
  disk_used = cost.disk;
}

void Extraction::Update(NexusMt *_mt) {
  mt = _mt;
  root = mt->history.Root();
  sink = root + (mt->history.n_nodes()-1);
  //clear statistics

  if(!visited.size()) {
    visited.resize(mt->history.n_nodes(), false);
    visited[0] = true;
  }
  
  Init();

  //first we coarse
  while(back.size()) {
    if(draw_used <= draw_max &&
       extr_used <= extr_max &&
       (front.size() && back.front().error > front.front().error)) 
      break;
    pop_heap(back.begin(), back.end(), greater<HeapNode>());
    HeapNode hnode = back.back();

    if(Visited(hnode.node)) {
      if(!Coarse(hnode)) { //push back on heap the heapnode
	push_heap(back.begin(), back.end(), greater<HeapNode>());
	break;
      }
    }
    back.pop_back();
  }

  while(front.size() && (*front.begin()).error > target_error) {
    pop_heap(front.begin(), front.end());
    HeapNode hnode = front.back();

    if(!Visited(hnode.node)) {
      if(!Refine(hnode.node)) {
	push_heap(front.begin(), front.end());
	break;
      }
    }
    front.pop_back();
  } 
  
  Select();
  draw_size = selected.size();


  //Preloading now
  for(unsigned int i = 0; i < 100; i++) {
    if(!front.size() && !back.size()) break;
    if((i%2) && front.size()) {
      pop_heap(front.begin(), front.end());
      HeapNode hnode = front.back();
      Node *node = hnode.node;
      front.pop_back();
      Node::iterator i;
      for(i = node->out_begin(); i != node->out_end(); i++) {
	Link &link = (*i);
	for(Link::iterator k = link.begin(); k != link.end(); k++) {
	  selected.push_back((*k).patch);
	}
      }
    } else if(back.size()) {
      pop_heap(back.begin(), back.end(), greater<HeapNode>());
      HeapNode hnode = back.back();
      Node *node = hnode.node;
      back.pop_back();
      Node::iterator i;
      for(i = node->in_begin(); i != node->in_end(); i++) {
	Link &link = (*i);
	for(Link::iterator k = link.begin(); k != link.end(); k++) {
	  selected.push_back((*k).patch);
	}
      }
    }
  }
}

float Extraction::GetRefineError(Node *node) {
  float maxerror = 0;
  Node::iterator i;
  for(i = node->in_begin(); i != node->in_end(); i++) {
    Link &link = *i;
    for(Link::iterator k = link.begin(); k != link.end(); k++) {
      Entry &entry = (*mt)[(*k).patch];
      float error =  metric->GetError(entry);
      if(error > maxerror) maxerror = error;
    }
  }
  return maxerror;
}

bool Extraction::Refine(Node *node) {
  //recursively refine parent if applicable.
  Node::iterator i;
  for(i = node->in_begin(); i != node->in_end(); i++) {
    Node *parent = (*i).node;
    if(!Visited(parent)) 
      if(!Refine(parent)) 
	return false;
  }

  Cost cost;
  Diff(node, cost);

  if(disk_used  + cost.disk > disk_max) return false;
  if(extr_used  + cost.extr > extr_max) return false;
  if(draw_used  + cost.draw > draw_max) return false;

  extr_used += cost.extr;
  draw_used += cost.draw;
  disk_used += cost.disk;

  visited[node - root] = true;


  //now add to the front children (unless sink node)

  for(i = node->out_begin(); i != node->out_end(); i++) {
    Link &link = *i;
    if(link.node == sink) continue; 
    float maxerror = GetRefineError(link.node);

    if(maxerror > target_error) 
      front.push_back(HeapNode((*i).node, maxerror));
  }
  push_heap(front.begin(), front.end());
  return true;
}

bool Extraction::Coarse(HeapNode &hnode) {
  Node *node = hnode.node;
  //recursively coarse children if applicable.
  Node::iterator i;
  for(i = node->out_begin(); i != node->out_end(); i++) {
    Node *child = (*i).node;
    float error = GetRefineError(child);
    HeapNode hchild(child, error);
    if(Visited(child)) 
      if(!Coarse(hchild)) return false;
  }


  Cost cost;
  Diff(node, cost);
  extr_used -= cost.extr;
  draw_used -= cost.draw;
  disk_used -= cost.disk;

  if(disk_used > disk_max) return false;

  visited[node - root] = false;

  //now add to the back parents (unless root node)
  for(i = node->in_begin(); i != node->in_end(); i++) {
    Link &link = *i;
    if(link.node == root) continue;
    float maxerror = GetRefineError(link.node);

    back.push_back(HeapNode(link.node, maxerror));
    push_heap(back.begin(), back.end(), greater<HeapNode>());
  }

  front.push_back(hnode);
  push_heap(front.begin(), front.end());
  return true;
}

void Extraction::Select() {
  selected.clear(); 
  Node *root = mt->history.Root();
  
  Node *nodes = mt->history.nodes;
  for(unsigned int i = 0; i < visited.size(); i++) {
    if(!visited[i]) continue;
    Node &node = nodes[i];

    Node::iterator n;
    for(n = node.out_begin(); n != node.out_end(); n++) {
      unsigned int n_out = (*n).node - root;
      if(!visited[n_out]) {
	Link &link = *n;
	for(Link::iterator k = link.begin(); k != link.end(); k++) {
	  unsigned int patch = (*k).patch;
	  selected.push_back(patch);
	}
      }
    }
  }
}

void Extraction::Visit(Node *node) {
  if(visited[node - root]) return;

  visited[node - root] = true;

  Node::iterator i;
  for(i = node->in_begin(); i != node->in_end(); i++) {
    unsigned int n_in = (*i).node - root;
    if(visited[n_in]) continue;
    Visit((*i).node);
  }

  Cost cost;
  Diff(node, cost);
  extr_used += cost.extr;
  draw_used += cost.draw;
  disk_used += cost.disk;
  
  for(i = node->out_begin(); i != node->out_end(); i++) {
    float maxerror = 0;
    Link &link = *i;
    for(Link::iterator k = link.begin(); k != link.end(); k++) {
      Entry &entry = (*mt)[(*k).patch];
      float error =  metric->GetError(entry);
      if(error > maxerror) maxerror = error;
    }
    //TODO this check may be dangerous for non saturating things...
    if(maxerror > target_error) {
      HeapNode hnode((*i).node, maxerror);
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

void Extraction::Diff(Node *node, Cost &cost) {
  Node::iterator i;
  for(i = node->in_begin(); i != node->in_end(); i++) {
    Link &link = *i;
    for(Link::iterator k = link.begin(); k != link.end(); k++) {
      unsigned int patch = (*k).patch;
      Entry &entry = (*mt)[patch];
      cost.extr -= entry.ram_size;
      vcg::Sphere3f &sphere = entry.sphere;
      if(!frustum.IsOutside(sphere.Center(), sphere.Radius()))
	cost.draw -= entry.ram_size;
      if(!entry.patch)
	cost.disk -= entry.disk_size;
    }
  }
  
  for(i = node->out_begin(); i != node->out_end(); i++) {
    Link &link = *i;
    for(Link::iterator k = link.begin(); k != link.end(); k++) {
      unsigned int patch = (*k).patch;
      Entry &entry = (*mt)[patch];
      cost.extr += entry.ram_size;
      vcg::Sphere3f &sphere = entry.sphere;
      if(!frustum.IsOutside(sphere.Center(), sphere.Radius()))
	cost.draw += entry.ram_size;
      if(!entry.patch)
	cost.disk += entry.disk_size;
    }
  }
}

