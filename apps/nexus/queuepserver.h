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
  unsigned int vbo_max;

  std::map<unsigned int, Data> index;
  std::vector<Item> heap;

  QueuePServer(): queue(64000) {}
  
  //Data &Lookup(unsigned int patch, unsigned short nv, unsigned short nf,
//		float priority = 0.0f);  

  //return flushing too.
  Data &Lookup(unsigned int patch, unsigned short nv, unsigned short nf, float priority, 
    std::vector<QueuePServer::Data> &data);

  bool IsLoaded(unsigned int patch);
  float MaxPriority();
  void Flush();

  void LoadVbo(Data &data);
  void FlushVbo(Data &data);
};

}

#endif
