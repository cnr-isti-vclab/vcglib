#include "preload.h"
#include "nexusmt.h"

using namespace std;
using namespace nxs;

 void Preload::execute() {
    assert(mt);
    while(!get_signaled()) {
      lock.enter();
      while(!queue.size()) {
        //cerr << "Acc nothing to preload!\n";
	      lock.leave();
	      pt::psleep(10);
	      lock.enter();
      }
      //TODO check we are not loading too much memory!
      assert(queue.size());
      unsigned int patch = queue.back();
      mt->GetPatch(patch, false);
      queue.pop_back();
      lock.leave();
    }
  }
  
