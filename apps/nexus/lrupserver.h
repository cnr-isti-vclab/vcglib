#ifndef NXS_LRU_PSERVER_H
#define NXS_LRU_PSERVER_H

#include <list>
#include <map>
#include <iostream>
#include "pserver.h"

namespace nxs {
  
class LruPServer: public PServer {
 public:
  //TODO change name to Item
  typedef std::pair<unsigned int, Patch *> Item;
  std::list<Item> items;
  typedef std::list<Item> Items;
  std::map<unsigned int, Items::iterator> index;
  
  ~LruPServer() {
    Flush();
  }

  Patch &Lookup(unsigned int patch) {    
    if(index.count(patch)) {
      Items::iterator &i = index[patch];
      Item item = *i;
      items.erase(i);
      items.push_front(item);
      i = items.begin();
      return *((*i).second);
    } else {          
      while(ram_used > ram_max) {
        Item item = items.back();        
	      index.erase(item.first);
	      FlushPatch(item.first);
	      items.pop_back();
      }
      Item item;
      item.first = patch;
      item.second = LoadPatch(patch);
      items.push_front(item);
      Items::iterator i = items.begin();
      index[patch] = i;
      return *((*i).second);
    }
  }

  bool IsLoaded(unsigned int patch) {
    return index.count(patch);
  }
  void Flush() {    
    std::map<unsigned int, Items::iterator>::iterator i;
    for(i = index.begin(); i != index.end(); i++) {
      Item &item = *((*i).second);
      FlushPatch((*i).first);
    }             
    for(int k = 0; k < size(); k++)
      operator[](k).patch = NULL;
        
    items.clear();
    index.clear();
  }
};

}//namespace

#endif
