#ifndef NXS_PRELOAD_H
#define NXS_PRELOAD_H

#include <assert.h>

#include <vector>

#include <ptypes/pasync.h>

namespace nxs {

class NexusMt;

class Preload: public pt::thread{
 public:

  NexusMt *mt;
  
  pt::mutex lock;
  
  std::vector<unsigned int> queue;

  Preload(): thread(false) {}
  ~Preload() {
    waitfor();
  }
  
  void execute();
  
  void post(std::vector<unsigned int> &patches) {
    lock.enter();
    queue = patches;
    lock.leave();
  }

  void cleanup() {}
};

}
#endif
