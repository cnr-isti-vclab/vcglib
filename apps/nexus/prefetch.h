#ifndef NXS_PREFETCH_H
#define NXS_PREFETCH_H

#include <map>
#include <vector>
#include <algorithm>

#include <ptypes/pasync.h>

#include "queuepserver.h"

namespace nxs {

class NexusMt;

class Prefetch: public pt::thread{
 public:
  
  pt::mutex safety;
  //unsigned int ram_max;
  //unsigned int ram_used;
  
  NexusMt *mt;
  std::vector<PServer::Item> missing;
  pt::jobqueue draw;
  pt::jobqueue load;

  Prefetch(): thread(false), draw(20000), load(64000) {}
  ~Prefetch() {
    waitfor();
  }
  
  void init(NexusMt *m, 
	    std::vector<unsigned int> &selected,
	    std::vector<PServer::Item> &visited);
  void execute();
  void cleanup() {}
};

}
#endif
