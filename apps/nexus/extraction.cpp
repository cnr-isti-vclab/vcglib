/****************************************************************************
* VCGLib                                                            o o     *
* Visual and Computer Graphics Library                            o     o   *
*                                                                _   O  _   *
* Copyright(C) 2004                                                \/)\/    *
* Visual Computing Lab                                            /\/|      *
* ISTI - Italian National Research Council                           |      *
*                                                                    \      *
* All rights reserved.                                                      *
*                                                                           *
* This program is free software; you can redistribute it and/or modify      *   
* it under the terms of the GNU General Public License as published by      *
* the Free Software Foundation; either version 2 of the License, or         *
* (at your option) any later version.                                       *
*                                                                           *
* This program is distributed in the hope that it will be useful,           *
* but WITHOUT ANY WARRANTY; without even the implied warranty of            *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
* GNU General Public License (http://www.gnu.org/licenses/gpl.txt)          *
* for more details.                                                         *
*                                                                           *
****************************************************************************/
/****************************************************************************
  History

$Log: not supported by cvs2svn $
Revision 1.8  2005/02/19 16:22:45  ponchio
Minor changes (visited and Cell)

Revision 1.7  2005/02/10 09:18:20  ponchio
Statistics.

Revision 1.6  2005/02/08 12:43:03  ponchio
Added copyright


****************************************************************************/

#include "extraction.h"
#include "metric.h"
#include "nexusmt.h"

#include <algorithm>

using namespace std;
using namespace nxs;

Extraction::Extraction(): target_error(4.0f), extr_max(10000),
			  draw_max(10000), disk_max(100) {
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
  visited.resize(mt->history.n_nodes(), 0);

  heap.clear();

  Visit(root); 

  while(heap.size()) {
    pop_heap(heap.begin(), heap.end());
    HeapNode hnode = heap.back();
    heap.pop_back();

    Node *node = hnode.node;
    if(Visited(node)) continue;

    if(Expand(hnode)) 
      Visit(node);
  }
  Select();
  draw_size = selected.size();
}

void Extraction::Init() { 
  max_error = -1;
  front.clear();
  back.clear();
  errors.clear();
  
  Cost cost;

  Node *nodes = mt->history.nodes;
  for(unsigned int i = 0; i < visited.size(); i++) {
    if(!visited[i]) continue;
    //    if(!visited[i]) continue;
    Node &node = nodes[i];

    bool cancoarse = true;
    Node::iterator n;
    for(n = node.out_begin(); n != node.out_end(); n++) {
      if(!Visited((*n).node)) {
	//      if(!visited[(*n).node - root]) {
	float maxerror = -1;
	
	Link &link = *n;
	for(Link::iterator k = link.begin(); k != link.end(); k++) {
	  //	  unsigned int patch = *k;
	  unsigned int patch = k;
	  Entry &entry = (*mt)[patch];
	  
	  bool visible;
	  float error = metric->GetError(entry, visible);

	  if(error > maxerror) maxerror = error;
	  
	  cost.extr += entry.ram_size;
	  //TODO if(frustum) else
	  //if(visible)
	  //	    cost.draw += entry.ram_size;
	  
	  vcg::Sphere3f &sphere = entry.sphere;
	  if(!frustum.IsOutside(sphere.Center(), sphere.Radius())) 
	    cost.draw += entry.ram_size;
			
	  if(!entry.patch)
	    cost.disk += entry.disk_size;
	}
	if((*n).node != sink && maxerror > target_error)
	  front.push_back(HeapNode((*n).node, maxerror));
	if(maxerror > max_error) max_error = maxerror;
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
    SetVisited(root, true);
  }
  
  Init();
  
  bool no_draw = false;  
  bool no_disk = false;  

  //TODO big problem: nodes a (error 10) with parent b (error -1) 
  //i try to refine a, refine b (recursive) but fail to refine a
  //next step i coarse b whis cause a cycle.

  while(1) {        
    if(!no_draw &&                       //we have buffer
       front.size() &&                   //we are not at max level
       front[0].error > target_error) {  //we are not already target_error            

      max_error = front[0].error;
      pop_heap(front.begin(), front.end());
      HeapNode hnode = front.back();            
      front.pop_back();
      if(!Visited(hnode.node)) {              
        if(!Refine(hnode)) {
          no_draw = true;
        }        
      }            
      continue;           
    }
 
    if(!back.size()) {
      //cerr << "nothing to coarse!\n";        
      break; //nothing to coarse (happen only on extr_max < root.extr
    }

    if(no_draw) { //suppose i have no more buffer      
      //if i do error damages coarsening better get out   
      if(front.size() && ((back.front().error + 0.001) >= front.front().error)) {
        //cerr << "Balanced cut\n";
        break;
      }
      if(!front.size() && back.front().error >= target_error) {
        //cerr << "Maxed out\n";
        break;
      }
    }

    //nothing to refine, coarse only if error <= target_error    
    if(!no_draw && back.front().error >= target_error) {
      //cerr << "error dominating\n";
      break;
    }    
    
    
    pop_heap(back.begin(), back.end(), greater<HeapNode>());
    HeapNode hnode = back.back();
    back.pop_back();
    if(Visited(hnode.node)) {
      bool recursive = false;
      Node::iterator i;
      for(i = hnode.node->out_begin(); i != hnode.node->out_end(); i++) {
         Node *child = (*i).node;
         if(Visited(child)) recursive = true;
      }
      if(!recursive && !Coarse(hnode)) { //no more disk so.. push back on heap the heapnode          
         back.push_back(hnode);
	       push_heap(back.begin(), back.end(), greater<HeapNode>());          
	       break;          
      }
    }
    
    no_draw = false;          
  }
  
  Select();
  draw_size = selected.size();

  //Preloading now
  for(unsigned int i = 0; i < 1000; i++) {
    if(!front.size() && !back.size()) break;
    if((i%2) && front.size()) {
      pop_heap(front.begin(), front.end());
      HeapNode hnode = front.back();
      Node *node = hnode.node;
      front.pop_back();
      Node::iterator l;
      for(l = node->out_begin(); l != node->out_end(); l++) {
	Link &link = (*l);
	for(Link::iterator k = link.begin(); k != link.end(); k++) {
	  //	  selected.push_back(Item(*k, i));
	  //	  errors[*k] = i;
	  selected.push_back(Item(k, i));
	  errors[k] = i;
	}
      }
    } else if(back.size()) {
      pop_heap(back.begin(), back.end(), greater<HeapNode>());
      HeapNode hnode = back.back();
      Node *node = hnode.node;
      back.pop_back();
      Node::iterator l;
      for(l = node->in_begin(); l != node->in_end(); l++) {
	Link &link = (*l);
	for(Link::iterator k = link.begin(); k != link.end(); k++) {
	  //	  selected.push_back(Item(*k, i));
	  //	  errors[*k] = i;
	  selected.push_back(Item(k, i));
	  errors[k] = i;
	}
      }
    }
  }
}

float Extraction::GetRefineError(Node *node) {
  float maxerror = -1;
  Node::iterator i;
  for(i = node->in_begin(); i != node->in_end(); i++) {
    Link &link = *i;
    for(Link::iterator p = link.begin(); p != link.end(); p++) {
      Entry &entry = (*mt)[p];
      bool visible;
      float error =  metric->GetError(entry, visible);
      if(error > maxerror) maxerror = error;
    }
  }
  return maxerror;
}

bool Extraction::Refine(HeapNode hnode) {
  
  Node *node = hnode.node;

  //recursively refine parent if applicable.
  Node::iterator i;
  for(i = node->in_begin(); i != node->in_end(); i++) {
    Node *parent = (*i).node;
    if(!Visited(parent)) {      
      //Here i use parent refine error!!!
      if(!Refine(HeapNode(parent, hnode.error))) return false;      
    }
  }

  Cost cost;
  Diff(node, cost);

  bool failed = false;  
  if(disk_used  + cost.disk > disk_max) {
    //cerr << "Disk failed\n";
    failed = true;    
  }
  if(extr_used  + cost.extr > extr_max) {
    //cerr << "Extr failed\n";
    failed = true;        
  }
  if(draw_used  + cost.draw > draw_max) {
    //cerr << "Draw failed\n";
    failed = true;        
  }

  if(failed) {
    front.push_back(hnode);
    push_heap(front.begin(), front.end());
    return false;
  }
  extr_used += cost.extr;
  draw_used += cost.draw;
  disk_used += cost.disk;

  SetVisited(node, true);


  //now add to the front children (unless sink node)

  for(i = node->out_begin(); i != node->out_end(); i++) {
    Link &link = *i;
    if(link.node == sink) continue; 
    float maxerror = GetRefineError(link.node);

    if(maxerror > target_error) 
      front.push_back(HeapNode((*i).node, maxerror));
  }
  push_heap(front.begin(), front.end());

  back.push_back(hnode);
  push_heap(back.begin(), back.end(), greater<HeapNode>());
  return true;
}

bool Extraction::Coarse(HeapNode hnode) {
  //cerr << "Coarse node: " << (void *)hnode.node << " err: " << hnode.error << endl;    
  Node *node = hnode.node;
  
  //recursively coarse children if applicable.
  Node::iterator i;
  for(i = node->out_begin(); i != node->out_end(); i++) {
    Node *child = (*i).node;
    float error = GetRefineError(child);
    HeapNode hchild(child, error);
    if(Visited(child)) {
      if(!Coarse(hchild)) return false;
    }
  }


  Cost cost;
  Diff(node, cost);
  extr_used -= cost.extr;
  draw_used -= cost.draw;
  disk_used -= cost.disk;

  if(disk_used > disk_max) return false;

  SetVisited(node, false);

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
	/*	for(Link::iterator k = link.begin(); k != link.end(); k++) {
	  unsigned int patch = *k;
	  selected.push_back(Item(patch,0));
	  errors[patch] = 0.0f;
	  }*/
	for(Link::iterator p= link.begin(); p != link.end(); p++) {
	  selected.push_back(Item(p, 0));
	  errors[p] = 0.0f;
	}
      }
    }
  }
}

void Extraction::Visit(Node *node) {
  assert(!Visited(node));

  SetVisited(node, true);

  Node::iterator i;
  for(i = node->in_begin(); i != node->in_end(); i++) {
    if(Visited((*i).node)) continue;
    Visit((*i).node);
  }

  Cost cost;
  Diff(node, cost);
  extr_used += cost.extr;
  draw_used += cost.draw;
  disk_used += cost.disk;
  
  for(i = node->out_begin(); i != node->out_end(); i++) {
    float maxerror = -1;
    Link &link = *i;
    for(Link::iterator p = link.begin(); p != link.end(); p++) {
      Entry &entry = (*mt)[p];
      bool visible;
      float error =  metric->GetError(entry, visible);
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
    for(Link::iterator p = link.begin(); p != link.end(); p++) {
      Entry &entry = (*mt)[p];
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
    for(Link::iterator p = link.begin(); p != link.end(); p++) {
      Entry &entry = (*mt)[p];
      cost.extr += entry.ram_size;
      vcg::Sphere3f &sphere = entry.sphere;
      if(!frustum.IsOutside(sphere.Center(), sphere.Radius()))
	cost.draw += entry.ram_size;
      if(!entry.patch)
	cost.disk += entry.disk_size;
    }
  }
}

