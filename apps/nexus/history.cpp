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
Revision 1.9  2005/02/20 00:43:23  ponchio
Less memory x extraction.  (removed frags)

Revision 1.8  2005/02/19 17:14:02  ponchio
History quick by default.

Revision 1.7  2005/02/19 16:22:45  ponchio
Minor changes (visited and Cell)

Revision 1.6  2005/02/17 16:40:35  ponchio
Optimized BuildLevels.

Revision 1.5  2005/02/08 12:43:03  ponchio
Added copyright


****************************************************************************/

#include <assert.h>
#include <map>
#include <iostream>

#include "history.h"
#include "nexus.h"


using namespace std;
using namespace nxs;

History::~History() {
  if(buffer) delete []buffer;
  nodes = NULL;
  in_links= NULL;
  out_links = NULL;
}

void History::Clear() {
  if(buffer) delete []buffer;
  buffer = NULL;
  updates.clear();
}

void History::ClearQuick() {
  if(buffer) delete []buffer;
  buffer = NULL;
}

void History::ClearUpdates() {
  updates.clear();
}

bool History::Load(unsigned int _size, char *mem) {
  if(buffer) delete []buffer;
  unsigned int is_quick = *(unsigned int *)mem;
  bool success;
  if(is_quick == 53) {
    success = LoadQuick(_size, mem);
  } else if(is_quick == 32) {
    success = LoadUpdates(_size, mem);
  } else {
    cerr << "Invalid history: " << is_quick << "\n";
    return false;
  }
  return success;
}

bool History::LoadQuick(unsigned int _size, char *mem) {
  buffer = mem;
  nodes = (Node *)(buffer + 4 * sizeof(int));
  in_links = (Link *)(nodes + n_nodes());
  out_links = in_links + n_in_links();

  //check size is ok;
  assert(n_nodes() * sizeof(Node) +
	 (n_in_links() + n_out_links()) * sizeof(Link) +
	 4 * sizeof(int) == _size);
  size = _size;
  return LoadPointers();
}

bool History::LoadUpdates(unsigned int _size, char *mem) {
  unsigned int *tmp = (unsigned int *)mem;
  updates.resize(tmp[1]);

  unsigned int pos = 2;
  for(unsigned int i = 0; i < updates.size(); i++) {
    unsigned int erased = tmp[pos++];
    unsigned int created = tmp[pos++];
    updates[i].erased.resize(erased);
    updates[i].created.resize(created);
    for(unsigned int e = 0; e < erased; e++) 
      updates[i].erased[e] = tmp[pos++];
    for(unsigned int e = 0; e < created; e++) 
      updates[i].created[e] = tmp[pos++];
  }
  delete []mem;
  buffer = 0;
  return true;
}

bool History::LoadPointers() {
  //now convert integer to pointers
  for(unsigned int i = 0; i < n_nodes(); i++) {
    Node &node = nodes[i];
    assert(((unsigned int)node.in_begin) <= n_in_links());
    assert(((unsigned int)node.out_begin) <= n_out_links());
    node.in_begin = in_links + (unsigned int)(node.in_begin);
    node.in_end = in_links + (unsigned int)(node.in_end);
    node.out_begin = out_links + (unsigned int)(node.out_begin);
    node.out_end = out_links + (unsigned int)(node.out_end);
  }
  
  for(unsigned int i = 0; i < n_in_links(); i++) {
    Link &link = in_links[i];
    assert(((unsigned int)link.node) <= n_nodes());
    link.node = nodes + (unsigned int)(link.node);
  } 

  for(unsigned int i = 0; i < n_out_links(); i++) {
    Link &link = out_links[i];
    assert(((unsigned int)link.node) <= n_nodes());
    link.node = nodes + (unsigned int)(link.node);
  }
  return true;
}

char *History::Save(unsigned int &_size) {
  if(buffer) 
    return SaveQuick(_size);
  else 
    return SaveUpdates(_size);
}

char *History::SaveQuick(unsigned int &_size) {
  assert(buffer);
  for(unsigned int i = 0; i < n_nodes(); i++) {
    Node &node = nodes[i];
    node.in_begin = (Link *)(node.in_begin - in_links);
    node.in_end   = (Link *)(node.in_end   - in_links);
    node.out_begin = (Link *)(node.out_begin - out_links);
    node.out_end   = (Link *)(node.out_end   - out_links);
  }
  
  for(unsigned int i = 0; i < n_in_links(); i++) {
    Link &link = in_links[i];
    link.node = (Node *)(link.node - nodes);
  }

  for(unsigned int i = 0; i < n_out_links(); i++) {
    Link &link = out_links[i];
    link.node = (Node *)(link.node - nodes);
  }

  assert(n_nodes() * sizeof(Node) +
	 (n_in_links() + n_out_links()) * sizeof(Link) +
	 4 * sizeof(int) == size);

  _size = size;
  char *tmp = buffer;
  buffer = NULL;
  return tmp;
}

char *History::SaveUpdates(unsigned int &_size) {
  vector<unsigned int> buf;
  buf.push_back(32);
  buf.push_back(updates.size());
  for(unsigned int i = 0; i < updates.size(); i++) {
    Update &update = updates[i];
    buf.push_back(update.erased.size());
    buf.push_back(update.created.size());
    for(unsigned int e = 0; e < update.erased.size(); e++)
      buf.push_back(update.erased[e]);
    for(unsigned int e = 0; e < update.created.size(); e++)
      buf.push_back(update.created[e]);
  }
  
  _size = buf.size() * sizeof(unsigned int);
  char *mem = new char[_size];
  memcpy(mem, &*buf.begin(), _size);
  return mem;
}

bool History::UpdatesToQuick(Nexus &nexus) {
  //maps cell -> node containing it
  map<unsigned int, unsigned int> cell_node;   
  //maps node -> Links
  map<unsigned int, vector<Link> > node_inlinks;
  map<unsigned int, vector<Link> > node_outlinks;

  vector<Node> tmp_nodes;
  tmp_nodes.resize(updates.size());

  vector<Link> tmp_in_links;
  vector<Link> tmp_out_links;
  vector<unsigned int> tmp_frags;

  unsigned int current_node = 0;

  vector<Update>::iterator u;
  for(u = updates.begin(); u != updates.end(); u++) {      
    Node &node = tmp_nodes[current_node];
    
    //created cells belong to this node, 
    for(unsigned int i = 0; i < (*u).created.size(); i++) {      
      unsigned int cell = (*u).created[i];
      cell_node[cell] = current_node;        
    }

    //Every erased cell already belonged to a node   //node -> its cells
    map<unsigned int, vector<unsigned int> > node_erased;      
    
    for(unsigned int i = 0; i < (*u).erased.size(); i++) {
      unsigned int cell = (*u).erased[i];
      assert(cell_node.count(cell));
      node_erased[cell_node[cell]].push_back(cell);               
    }      
    
    //for every node with erased cells we build a fragment and 
    //put the corresponding cells in it.      
    map<unsigned int, vector<unsigned int> >::iterator e;
    for(e = node_erased.begin(); e != node_erased.end(); e++) {

      unsigned int floor_node = (*e).first;
      vector<unsigned int> &cells = (*e).second;

      Node &parent = tmp_nodes[floor_node];

      Link inlink;
      inlink.node = (Node *)floor_node;
      inlink.begin = tmp_frags.size();
      inlink.end = inlink.begin + cells.size();

      Link outlink;
      outlink.node = (Node *)current_node;
      outlink.begin = tmp_frags.size();
      outlink.end = outlink.begin + cells.size();

      //Fill it with erased cells.
      vector<unsigned int>::iterator k;
      for(k = cells.begin(); k != cells.end(); k++) 
	tmp_frags.push_back(*k);

      //Add the new Frag to the nodes (in and out)
      node_outlinks[floor_node].push_back(outlink);
      node_inlinks[current_node].push_back(inlink);
    }
    current_node++; 
  }
  
  map<unsigned int, vector<Link> >::iterator k;
  for(k = node_outlinks.begin(); k != node_outlinks.end(); k++) {
    unsigned int inode = (*k).first;
    vector<Link> &links = (*k).second;
    tmp_nodes[inode].out_begin = (Link *)(tmp_out_links.size());
    tmp_nodes[inode].out_end = (Link *)(tmp_out_links.size() + links.size());
    //    tmp_nodes[inode].out_link_begin = (Link *)(tmp_out_links.size());
    //    tmp_nodes[inode].out_link_size = links.size();
    
    for(unsigned int i = 0; i < links.size(); i++) 
      tmp_out_links.push_back(links[i]);
  }
  
  for(k = node_inlinks.begin(); k != node_inlinks.end(); k++) {
    unsigned int inode = (*k).first;
    vector<Link> &links = (*k).second;
    //    tmp_nodes[inode].in_link_begin = (Link *)(tmp_in_links.size());
    //    tmp_nodes[inode].in_link_size = links.size();
    tmp_nodes[inode].in_begin = (Link *)(tmp_in_links.size());
    tmp_nodes[inode].in_end = (Link *)(tmp_in_links.size() + links.size());
    
    for(unsigned int i = 0; i < links.size(); i++) 
      tmp_in_links.push_back(links[i]);
  }
  
  //Here we reorder entries in nexus...
  nexus.Flush();

  vector<Entry> entries;
  entries.resize(nexus.size());
  for(unsigned int i = 0; i < nexus.size(); i++) {
    assert(!nexus[i].patch);
    entries[i] = nexus[i];
  }
  assert(tmp_frags.size() == nexus.size());
  for(unsigned int i = 0; i < tmp_frags.size(); i++) {
    nexus[i] = entries[tmp_frags[i]];
  }
  //WARNING CRITICAL TODOHey we should do the same on the borders!
  vector<unsigned int> backward;
  backward.resize(tmp_frags.size());
  for(unsigned int i = 0; i < backward.size(); i++) 
    backward[tmp_frags[i]] = i;

  for(unsigned int i = 0; i < nexus.borders.size(); i++) {
    Border &border = nexus.borders.GetBorder(i);
    for(unsigned int k = 0; k < border.Size(); k++)
      border[k].end_patch = backward[border[k].end_patch];
  }
  nexus.borders.Flush();
  vector<Border> borders;
  borders.resize(nexus.borders.size());
  for(unsigned int i = 0; i < nexus.borders.size(); i++) {
    borders[i] = nexus.borders[i];
  }
  assert(tmp_frags.size() == nexus.borders.size());
  for(unsigned int i = 0; i < tmp_frags.size(); i++) {
    nexus.borders[i] = borders[tmp_frags[i]];
  }

  
  size = tmp_nodes.size() * sizeof(Node) +
    tmp_in_links.size() * sizeof(Link) + 
    tmp_out_links.size() * sizeof(Link) + 
    4 * sizeof(int);
  
  if(buffer) delete []buffer;
  buffer = new char[size];

  quick() = 53;
  n_nodes() = tmp_nodes.size();
  n_in_links() = tmp_in_links.size();
  n_out_links() = tmp_out_links.size();

  nodes = (Node *)(buffer + 4 * sizeof(int));
  in_links = (Link *)(nodes + n_nodes());
  out_links = in_links + n_in_links();
  
  memcpy(nodes, &*tmp_nodes.begin(), tmp_nodes.size()*sizeof(Node));
  memcpy(in_links, &*tmp_in_links.begin(), tmp_in_links.size()*sizeof(Link));
  memcpy(out_links, &*tmp_out_links.begin(), 
	 tmp_out_links.size()*sizeof(Link));

  return LoadPointers();
}

void History::BuildLevels(vector<int> &levels) {
  levels.clear();
  if(buffer) {
    //Saved in quick mode:
    for(unsigned int n = 0; n < n_nodes(); n++) {
      Node *node = nodes+n;
      Node::iterator l;
      unsigned int current = 0;
      if(node != nodes) { //not root
	Link *inlink = node->in_begin;
	unsigned int p = inlink->begin;
	assert(p < levels.size());
	assert(p >= 0);
	current = levels[p]+1;
      }
      for(l = node->out_begin; l != node->out_end; l++) {
	Link &link = *l;
	for(unsigned int p = link.begin; p != link.end; p++) {
	  while(p >= levels.size()) levels.push_back(-1);
	  levels[p] = current;
	}
      }
    }
  } else {
    //Saved in updates mode:
    for(unsigned int i = 0; i < updates.size(); i++) {
      Update &u = updates[i];
      unsigned int current = 0;
      if(!u.erased.size()) current = 0;
      else current = levels[u.erased[0]] + 1;
      for(unsigned int i = 0; i < u.created.size(); i++) {
	unsigned int p = u.created[i];
	while(p >= levels.size()) levels.push_back(-1);
	levels[p] = current;
      }
    }
  }
}
