#include "nexusmt.h"
#include <map>
#include <queue>

using namespace nxs;
using namespace std;

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
