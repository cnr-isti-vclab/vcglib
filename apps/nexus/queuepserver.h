#ifndef NXS_QUEUE_PSERVER_H
#define NXS_QUEUE_PSERVER_H

#include <map>
#include <hash_map>
#include <vector>
#include <algorithm>

#include <ptypes/pasync.h>

#include "pserver.h"

namespace nxs {


class QueuePServer: public PServer {
 public:

  enum Action { DRAW = 1, FLUSH = 0 };
  struct Data {
    Patch *patch;
    unsigned int vbo_array;
    unsigned int vbo_element;
    Data(): patch(NULL), vbo_array(0), vbo_element(0) {}
  };

  pt::jobqueue queue;  
  unsigned int vbo_used;

  typedef std::pair<unsigned int, Data> Item;
  std::list<Item> items;
  typedef std::list<Item> Items;
  std::map<unsigned int, Items::iterator> index;
  
  QueuePServer(): queue(64000), vbo_used(0) {}
  ~QueuePServer() {
    Flush();
  }

   Data &Lookup(unsigned int patch, unsigned short nv, unsigned short nf,
                std::vector<Data> &flush) {
    if(index.count(patch)) {
      Items::iterator &i = index[patch];
      Item item = *i;
      items.erase(i);
      items.push_front(item);
      i = items.begin();
      return ((*i).second);
    } else {
      while(ram_used > ram_max) {    
        Data &data = items.back().second;
        //TODO i should not flush current extraction!
	      index.erase(items.back().first);        
	      FlushPatch(patch, data.patch);
        flush.push_back(data);
	      items.pop_back();
      }
      Item item;
      item.first = patch;
      item.second.patch = LoadPatch(patch, nv, nf);
      items.push_front(item);
      Items::iterator i = items.begin();
      index[patch] = i;      
      return ((*i).second); 
    }                        
  }
         
  
  //return flushing too.
  //Data &Lookup(unsigned int patch, unsigned short nv, unsigned short nf, float priority, 
  //    std::vector<QueuePServer::Data> &data);

  //bool IsLoaded(unsigned int patch);
  //float MaxPriority();
                
  
  bool IsLoaded(unsigned int patch) { 
    return index.count(patch);
  }
  void Flush() {
    std::map<unsigned int, Items::iterator>::iterator i;
    for(i = index.begin(); i != index.end(); i++) {
      Item &item = *((*i).second);
      FlushVbo(item.second);
      FlushPatch((*i).first, item.second.patch);
    }             
    for(int k = 0; k < entries.size(); k++)
      entries[k].patch = NULL;
        
    items.clear();
    index.clear();
  }
  void LoadVbo(Data &data);
  void FlushVbo(Data &data);
};

}

#endif
