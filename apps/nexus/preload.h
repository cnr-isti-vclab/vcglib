#ifndef NXS_PRELOAD_H
#define NXS_PRELOAD_H

#include <assert.h>

#include <vector>

#include <ptypes/pasync.h>

#include "extraction.h"

namespace nxs {

class NexusMt;

class Preload: public pt::thread{
 public:

  NexusMt *mt;
  
  pt::mutex lock;
  pt::trigger trigger;
  
  std::vector<Item> queue;

  Preload(): thread(false), trigger(false, false) {}
  ~Preload() {
    waitfor();
  }
  
  void execute();
  
  void post(std::vector<Item> &patches) {
    trigger.reset();
    lock.enter();

    queue.reserve(patches.size());
    for(int i = patches.size() -1; i >= 0; i--)
      queue.push_back(patches[i]);

    trigger.post();
    lock.leave();
  }

  void cleanup() {}
};

}
#endif
