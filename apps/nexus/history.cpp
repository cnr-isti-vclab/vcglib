#include <assert.h>
#include <map>
#include <iostream>

#include "history.h"


using namespace std;
using namespace nxs;

History::~History() {
  if(buffer) delete []buffer;
  nodes = NULL;
  in_links= NULL;
  out_links = NULL;
  frags = NULL;
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
    //    cerr << "Load quick!\n";
    success = LoadQuick(_size, mem);
  } else if(is_quick == 32) {
    //    cerr << "Load updates\n";
    success = LoadUpdates(_size, mem);
  } else {
    cerr << "Invalid history: " << is_quick << "\n";
    return false;
  }
  return success;
}

bool History::LoadQuick(unsigned int _size, char *mem) {
  buffer = mem;
  nodes = (Node *)(buffer + 5 * sizeof(int));
  in_links = (Link *)(nodes + n_nodes());
  out_links = in_links + n_in_links();
  frags = (Cell *)(out_links + n_out_links());

  //check size is ok;
  assert(n_nodes() * sizeof(Node) +
	 (n_in_links() + n_out_links()) * sizeof(Link) +
	 n_frags() * sizeof(Cell) + 
	 5 * sizeof(int) == size);
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
    assert(((unsigned int)node.in_link_begin) <= n_in_links());
    assert(((unsigned int)node.out_link_begin) <= n_out_links());
    node.in_link_begin = in_links + (unsigned int)(node.in_link_begin);
    node.out_link_begin = out_links + (unsigned int)(node.out_link_begin);
  }
  
  for(unsigned int i = 0; i < n_in_links(); i++) {
    Link &link = in_links[i];
    assert(((unsigned int)link.node) <= n_nodes());
    assert(((unsigned int)link.frag_begin) <= n_frags());
    link.node = nodes + (unsigned int)(link.node);
    link.frag_begin = frags + (unsigned int)(link.frag_begin);
  } 

  for(unsigned int i = 0; i < n_out_links(); i++) {
    Link &link = out_links[i];
    assert(((unsigned int)link.node) <= n_nodes());
    assert(((unsigned int)link.frag_begin) <= n_frags());
    link.node = nodes + (unsigned int)(link.node);
    link.frag_begin = frags + (unsigned int)(link.frag_begin);
  }
  return true;
}

char *History::Save(unsigned int &_size) {
  if(buffer) {
    //    cerr << "SaveQuick!\n";
    return SaveQuick(_size);
  } else {
    //    cerr << "Save updates\n";
    return SaveUpdates(_size);
  }
}

char *History::SaveQuick(unsigned int &_size) {
  assert(buffer);
  for(unsigned int i = 0; i < n_nodes(); i++) {
    Node &node = nodes[i];
    node.in_link_begin = (Link *)(node.in_link_begin - in_links);
    node.out_link_begin = (Link *)(node.out_link_begin - out_links);
  }
  
  for(unsigned int i = 0; i < n_in_links(); i++) {
    Link &link = in_links[i];
    link.node = (Node *)(link.node - nodes);
    link.frag_begin = (Cell *)(link.frag_begin - frags);
  }

  for(unsigned int i = 0; i < n_out_links(); i++) {
    Link &link = out_links[i];
    link.node = (Node *)(link.node - nodes);
    link.frag_begin = (Cell *)(link.frag_begin - frags);
  }

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

bool History::UpdatesToQuick() {
  //maps cell -> node containing it
  map<unsigned int, unsigned int> cell_node;   
  //maps node -> Links
  map<unsigned int, vector<Link> > node_inlinks;
  map<unsigned int, vector<Link> > node_outlinks;

  vector<Node> tmp_nodes;
  tmp_nodes.resize(updates.size());

  vector<Link> tmp_in_links;
  vector<Link> tmp_out_links;
  vector<Cell> tmp_frags;

  unsigned int current_node = 0;

  vector<Update>::iterator u;
  for(u = updates.begin(); u != updates.end(); u++) {      
    
    Node &node = tmp_nodes[current_node];
    
    //created cells belong to this node, we look also for max error.      
    for(unsigned int i = 0; i < (*u).created.size(); i++) {      
      unsigned int cell = (*u).created[i];
      cell_node[cell] = current_node;        
    }

    //Every erased cell already belonged to a node.
    //node -> its cells
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

      //node.in.push_back(innodes.size());

      unsigned int floor_node = (*e).first;
      vector<unsigned int> &cells = (*e).second;

      Node &parent = tmp_nodes[floor_node];

      Link inlink;
      inlink.node = (Node *)floor_node;
      inlink.frag_begin = (Cell *)(tmp_frags.size());
      inlink.frag_size = cells.size();

      Link outlink;
      outlink.node = (Node *)current_node;
      outlink.frag_begin = (Cell *)(tmp_frags.size());
      outlink.frag_size = cells.size();

      //Fill it with erased cells.

      vector<unsigned int>::iterator k;
      for(k = cells.begin(); k != cells.end(); k++) {
	      Cell cell;
	      cell.patch = (*k);
	      tmp_frags.push_back(cell);
      }         

      //Add the new Frag to the node.

      node_outlinks[floor_node].push_back(outlink);
      node_inlinks[current_node].push_back(inlink);
      
      //Update in and out of the nodes.
    }
    current_node++; 
  }
  
  map<unsigned int, vector<Link> >::iterator k;
  for(k = node_outlinks.begin(); k != node_outlinks.end(); k++) {
    unsigned int inode = (*k).first;
    vector<Link> &links = (*k).second;
    tmp_nodes[inode].out_link_begin = (Link *)(tmp_out_links.size());
    tmp_nodes[inode].out_link_size = links.size();
    
    for(unsigned int i = 0; i < links.size(); i++) 
      tmp_out_links.push_back(links[i]);
  }
  
  for(k = node_inlinks.begin(); k != node_inlinks.end(); k++) {
    unsigned int inode = (*k).first;
    vector<Link> &links = (*k).second;
    tmp_nodes[inode].in_link_begin = (Link *)(tmp_in_links.size());
    tmp_nodes[inode].in_link_size = links.size();
    
    for(unsigned int i = 0; i < links.size(); i++) 
      tmp_in_links.push_back(links[i]);
  }
  size = tmp_nodes.size() * sizeof(Node) +
    tmp_in_links.size() * sizeof(Link) + 
    tmp_out_links.size() * sizeof(Link) + 
    tmp_frags.size() * sizeof(Cell) +
    5 * sizeof(int);
  
  if(buffer) delete []buffer;
  buffer = new char[size];

  quick() = 53;
  n_nodes() = tmp_nodes.size();
  n_in_links() = tmp_in_links.size();
  n_out_links() = tmp_out_links.size();
  n_frags() = tmp_frags.size();

  nodes = (Node *)(buffer + 5 * sizeof(int));
  in_links = (Link *)(nodes + n_nodes());
  out_links = in_links + n_in_links();
  frags = (Cell *)(out_links + n_out_links());
  
  memcpy(nodes, &*tmp_nodes.begin(), tmp_nodes.size()*sizeof(Node));
  memcpy(in_links, &*tmp_in_links.begin(), tmp_in_links.size()*sizeof(Link));
  memcpy(out_links, &*tmp_out_links.begin(), 
	 tmp_out_links.size()*sizeof(Link));
  memcpy(frags, &*tmp_frags.begin(), tmp_frags.size() * sizeof(Cell));

  return LoadPointers();
}
