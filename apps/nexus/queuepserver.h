#ifndef NXS_QUEUE_PSERVER_H
#define NXS_QUEUE_PSERVER_H

#include <map>
#include <vector>
#include <algorithm>

#include "pserver.h"

namespace nxs {

class QueuePServer: public PServer {
 public:

  struct Data {
    Patch *patch;
    unsigned int vbo_array;
    unsigned int vbo_element;
    Data(): patch(NULL), vbo_array(0), vbo_element(0) {}
  };

  unsigned int vbo_used;
  unsigned int vbo_max;

  std::map<unsigned int, Data> index;
  std::vector<Item> heap;
  
  Data &Lookup(unsigned int patch, unsigned short nv, unsigned short nf,
		float priority = 0.0f);

  bool IsLoaded(unsigned int patch);
  void Flush();

 protected:

  void LoadVbo(Data &data);
  void FlushVbo(Data &data);
};

}

#endif
