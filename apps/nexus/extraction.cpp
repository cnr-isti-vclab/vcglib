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
Revision 1.12  2005/03/02 10:40:17  ponchio
Extraction rewrittten (to fix recusive problems).

Revision 1.11  2005/02/20 19:49:44  ponchio
cleaning (a bit more).

Revision 1.10  2005/02/20 18:07:00  ponchio
cleaning.

Revision 1.9  2005/02/20 00:43:23  ponchio
Less memory x extraction.  (removed frags)

Revision 1.8  2005/02/19 16:22:45  ponchio
Minor changes (visited and Cell)

Revision 1.7  2005/02/10 09:18:20  ponchio
Statistics.

Revision 1.6  2005/02/08 12:43:03  ponchio
Added copyright


****************************************************************************/

#include "extraction.h"
#include "metric.h"
#include "nexus.h"

#include <algorithm>

using namespace std;
using namespace nxs;

  /* Updateing strategy:
     if i can refine (not at leaves, 
                      have draw and extr buffer, 
                      not past target_error)
        i try to refine BUT
             i can fail because i finish some buffer 
                  (then i put the operation back on the stack)
	     if i have finished disk i should just quit
     if i cannot refine i consider coarsening:
         i need 1) not be at root (eheh)
                2) have finished draw and extr buffer 
                        (unless i am at error < target so i want to coarse
                3) do not make global error worse 
		   (unless it is error < target_error...)
		   a stability term is added (1.1) so that we do not flip
		   nodes in and out quickly
	 i try to coarse BUT
              i can fail because i need disk
	      (then i put the operation back on the stack)
      if i cannt coarse i just quit
  */



Extraction::Extraction(): target_error(4.0f), extr_max(10000),
			  draw_max(10000), disk_max(100) {
  metric = new FrustumMetric;
}

Extraction::~Extraction() {
  if(metric) delete metric;
}

void Extraction::SetMetric(Metric *m) {
  if(metric) 
    delete metric;
  metric = m;
}

void Extraction::Extract(Nexus *_mt) {
  mt = _mt;
  root = mt->history.Root();
  sink = root + (mt->history.n_nodes()-1);

  //clear statistics
  extr_used = draw_used = disk_used = 0;
  
  //first we clear the visited flags
  visited.clear();
  visited.resize(mt->history.n_nodes(), false);
  visible.clear();
  visible.resize(mt->size(), true);
  node_errors.clear();
  node_errors.resize(mt->history.n_nodes(), -1);

  front.clear();

  Visit(root); 

  while(front.size()) {
    pop_heap(front.begin(), front.end());
    HeapNode hnode = front.back();
    front.pop_back();

    Node *node = hnode.node;
    if(Visited(node)) continue;

    if(Expand(hnode)) 
      Visit(node);
  }
  Select();
  draw_size = selected.size();
}

void Extraction::Init() { 
  //I want to add all coarsable nodes
  //and all refinable node (being careful about recursive dependencies)
  for(Node *node = root; node != sink; node++) {
    if(!Visited(node)) continue;
    if(node != root && CanCoarse(node)) 
      back.push_back(HeapNode(node, GetNodeError(node)));

    for(Node::iterator n = node->out_begin; n != node->out_end; n++) {
      Node *child = n->node;
      if(Visited(child)) continue;
      if(node_errors[child - root] != -1) continue; //already visited

      float *error = GetNodeError(child);
      if(CanRefine(node)) // TODO? && error > target_error
	front.push_back(HeapNode(child, error));

      if(*error > max_error) max_error = *error;
    }
  }

  //recursively fix error 
  for(Node *node = root; node != sink; node++) 
    if(node_errors[node - root] != -1)
      SetError(node, node_errors[node-root]);       


  //estimate cost of all the cut arcs (i need the visible info)
  Cost cost;
  for(Node *node =  root; node != sink; node++) {
    if(!Visited(node)) continue;

    for(Node::iterator n = node->out_begin; n != node->out_end; n++) {
      Link &link = *n;
      if(Visited((*n).node)) continue;
      for(unsigned int patch = link.begin; patch != link.end; patch++) {
	Entry &entry = (*mt)[patch];
	                   cost.extr += entry.ram_size;
	if(Visible(patch)) cost.draw += entry.ram_size;
	if(!entry.patch)   cost.disk += entry.disk_size;
      }
    }
  }

  make_heap(front.begin(), front.end());
  make_heap(back.begin(), back.end(), greater<HeapNode>());

  extr_used = cost.extr;
  draw_used = cost.draw;
  disk_used = cost.disk;
}

void Extraction::Update(Nexus *_mt) {
  mt = _mt;
  root = mt->history.Root();
  sink = mt->history.Sink();

  if(!visited.size()) {
    visited.resize(mt->history.n_nodes(), false);
    SetVisited(root, true);
  } 
  visible.clear();
  visible.resize(mt->size(), true);
  node_errors.clear();
  node_errors.resize(mt->history.n_nodes(), -1);
  
  front.clear();
  back.clear();

  max_error = -1;
  
  Init();
  
  bool can_refine = true;

  while(1) {     
    while(can_refine) {                          //we have available budget
      if(!front.size()) break;                   //we are at max level
      if(*front[0].error < target_error) break;  //already at target_error

      max_error = *front[0].error;
      pop_heap(front.begin(), front.end());
      HeapNode hnode = front.back();            
      front.pop_back();

      if(!Visited(hnode.node) && CanRefine(hnode.node)) {
	if(!Refine(hnode)) {
	  can_refine = false;
	  front.push_back(hnode);
	  push_heap(front.begin(), front.end());
	}
      }
    }
 
    if(!back.size()) //nothing to coarse (happen only on extr_max < root.extr)
      break; 
    
    if(*back.front().error >= target_error &&
       (!front.size() ||
	(*back.front().error * 1.4) >= *front.front().error))
      break;
    
    pop_heap(back.begin(), back.end(), greater<HeapNode>());
    HeapNode hnode = back.back();
    back.pop_back();

    if(Visited(hnode.node) &&    //not already coarsed
       CanCoarse(hnode.node) &&  //all children !visited
       !Coarse(hnode)) {         //no more disk
      back.push_back(hnode);
      push_heap(back.begin(), back.end(), greater<HeapNode>());          
      break;          
    }
    can_refine = true;
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
      for(l = node->out_begin; l != node->out_end; l++) {
	Link &link = (*l);
	for(unsigned int k = link.begin; k != link.end; k++) {
	  selected.push_back(Item(k, i));
	}
      }
    } else if(back.size()) {
      pop_heap(back.begin(), back.end(), greater<HeapNode>());
      HeapNode hnode = back.back();
      Node *node = hnode.node;
      back.pop_back();

      for(Node::iterator l = node->in_begin; l != node->in_end; l++) {
	Link &link = (*l);
	for(unsigned int k = link.begin; k != link.end; k++) {
	  selected.push_back(Item(k, i));
	}
      }
    }
  }
}


float *Extraction::GetNodeError(Node *node) {
  float &maxerror = node_errors[node-root];
  for(Node::iterator i = node->in_begin; i != node->in_end; i++) {
    Link &link = *i;
    for(unsigned int p = link.begin; p != link.end; p++) {
      Entry &entry = (*mt)[p];
      bool visible;
      float error =  metric->GetError(entry, visible);
      //      cerr << "Error for patch: " << p << " -> " << error << endl;
      if(error > maxerror) maxerror = error;
      SetVisible(p, visible);
    }
  }
  return &maxerror;
}

bool Extraction::Refine(HeapNode hnode) {
  
  Node *node = hnode.node;
  
  Cost cost;
  Diff(node, cost);
  
  if(disk_used  + cost.disk > disk_max ||
     extr_used  + cost.extr > extr_max ||
     draw_used  + cost.draw > draw_max)
    return false;

  extr_used += cost.extr;
  draw_used += cost.draw;
  disk_used += cost.disk;

  SetVisited(node, true);

  //now add to the front children (unless sink node)

  for(Node::iterator i = node->out_begin; i != node->out_end; i++) {
    Link &link = *i;
    if(link.node == sink) continue; 

    float *error = &node_errors[link.node - root];
    if(*error == -1)
      error = GetNodeError(link.node);
    if(*hnode.error < *error) *hnode.error = *error;
    //TODO    if(maxerror > target_error) 
    if(CanRefine((*i).node)) {
      front.push_back(HeapNode((*i).node, error));
      push_heap(front.begin(), front.end());
    }
  }

  back.push_back(hnode);
  push_heap(back.begin(), back.end(), greater<HeapNode>());
  return true;
}

bool Extraction::Coarse(HeapNode hnode) {
  Node *node = hnode.node;
  
  Cost cost;
  Diff(node, cost);

  extr_used -= cost.extr;
  draw_used -= cost.draw;
  disk_used -= cost.disk;

  if(disk_used > disk_max) return false;

  SetVisited(node, false);

  //now add to the back parents (unless root node)
  for(Node::iterator i = node->in_begin; i != node->in_end; i++) {
    Link &link = *i;
    if(link.node == root) continue;

    float *error = &node_errors[link.node - root];
    if(*error == -1)
      error = GetNodeError(link.node);
    if(*error < *hnode.error) *error = *hnode.error;

    if(CanCoarse(link.node)) {
      back.push_back(HeapNode(link.node, error));
      push_heap(back.begin(), back.end(), greater<HeapNode>());
    }
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
    for(n = node.out_begin; n != node.out_end; n++) {
      unsigned int n_out = (*n).node - root;
      if(!visited[n_out]) {
	Link &link = *n;
	for(unsigned int p = link.begin; p != link.end; p++) {
	  selected.push_back(Item(p, 0));
	}
      }
    }
  }
}

void Extraction::Visit(Node *node) {
  assert(!Visited(node));

  SetVisited(node, true);

  for(Node::iterator i = node->in_begin; i != node->in_end; i++) {
    if(Visited((*i).node)) continue;
    Visit((*i).node);
  }

  Cost cost;
  Diff(node, cost);
  extr_used += cost.extr;
  draw_used += cost.draw;
  disk_used += cost.disk;
  
  for(Node::iterator i = node->out_begin; i != node->out_end; i++) {
    float maxerror = -1;
    Link &link = *i;
    for(unsigned int p = link.begin; p != link.end; p++) {
      Entry &entry = (*mt)[p];
      bool visible;
      float error =  metric->GetError(entry, visible);
      if(error > maxerror) maxerror = error;
      SetVisible(p, visible);
    }
    //TODO this check may be dangerous for non saturating things...
    if(maxerror > target_error) {
      node_errors[(*i).node - root] = maxerror;
      HeapNode hnode((*i).node, &node_errors[(*i).node - root]);
      front.push_back(hnode);
      push_heap(front.begin(), front.end());
    }
  }
}

bool Extraction::Expand(HeapNode &node) {
  if(extr_used >= extr_max) return false; 
  if(draw_used >= draw_max) return false;
  //  if(disk_used >= disk_max) return false;
  return *node.error > target_error; 
}

void Extraction::Diff(Node *node, Cost &cost) {
  Node::iterator i;
  for(i = node->in_begin; i != node->in_end; i++) {
    Link &link = *i;
    for(unsigned int p = link.begin; p != link.end; p++) {
      Entry &entry = (*mt)[p];
                       cost.extr -= entry.ram_size;
      if(Visible(p))   cost.draw -= entry.ram_size;
      if(!entry.patch) cost.disk -= entry.disk_size;
    }
  }
  
  for(i = node->out_begin; i != node->out_end; i++) {
    Link &link = *i;
    for(unsigned int p = link.begin; p != link.end; p++) {
      Entry &entry = (*mt)[p];
                       cost.extr += entry.ram_size;
      if(Visible(p))   cost.draw += entry.ram_size;
      if(!entry.patch) cost.disk += entry.disk_size;
    }
  }
}

 void Extraction::SetError(Node *node, float error) {
   for(Node::iterator i = node->in_begin; i != node->in_end; i++)
     if(node_errors[(*i).node - root] != -1 &&
	node_errors[(*i).node - root] < error) {
       node_errors[(*i).node - root] = error;
       SetError((*i).node, error);
     }
 }

 bool Extraction::CanRefine(Node *node) {
   for(Node::iterator i = node->in_begin; i != node->in_end; i++)
     if(!Visited((*i).node))
       return false;
   return true;
 }

 bool Extraction::CanCoarse(Node *node) {
   for(Node::iterator i = node->out_begin; i != node->out_end; i++)
     if(Visited((*i).node))
       return false;
   return true;
 }
