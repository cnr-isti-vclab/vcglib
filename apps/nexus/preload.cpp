#include "preload.h"
#include "nexusmt.h"
#include <iostream>

using namespace std;
using namespace nxs;

 void Preload::execute() {
   assert(mt);
   while(!get_signaled()) {
     trigger.wait();
     lock.enter();
     while(!queue.size()) {
       trigger.reset();
       lock.leave();
       trigger.wait();
       lock.enter();
     }
     //TODO check we are not loading too much memory!
     assert(queue.size());
     Item &item = queue.back();
     if(item.error == 0 || mt->CanAdd(item)) {
       //we cannot flush since we are not in the openGL thread
       //and flushing includes VBO buffer flushing also.
       mt->GetPatch(item.id, item.error, false);
       queue.pop_back();
     } else
       queue.clear();
     lock.leave();
   }
 }
  
