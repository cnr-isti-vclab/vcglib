#include "prefetch.h"
#include "nexusmt.h"

using namespace std;
using namespace nxs;


void Prefetch::init(NexusMt *m, std::vector<unsigned int> &selected,
		    std::vector<PServer::Item> &visited) {

  //  cerr << "Init\n";

  safety.lock();
  mt = m;
  missing.clear();
  
  unsigned int notloaded = 0;
  std::map<unsigned int, float> tmp;
  for(unsigned int i = 0; i < selected.size(); i++) {
    unsigned int patch = selected[i];
    tmp[patch] = 0.0f;
    if(!mt->patches.IsLoaded(patch))
      notloaded++;
    //    if(mt->patches.IsLoaded(patch))
    //      mt->todraw.push_back(make_pair(patch, (Patch *)NULL));
    //    else
    missing.push_back(PServer::Item(patch, 0.0f));
  }
  if(notloaded)
    cerr << "Patches to load: " << notloaded << endl;

  for(unsigned int i = 0; i < visited.size(); i++) {
    PServer::Item &item = visited[i];
    if(tmp.count(item.patch)) continue;
    if(mt->patches.IsLoaded(item.patch))
      tmp[item.patch] = item.priority;
    
    missing.push_back(item);
  }

  QueuePServer &ps = mt->patches;
  for(unsigned int i = 0; i < ps.heap.size(); i++) {
    PServer::Item &item = ps.heap[i];
    if(tmp.count(item.patch)) 
      item.priority = tmp[item.patch];
    else
      item.priority = 1e40;
  }
  make_heap(ps.heap.begin(), ps.heap.end());

  sort(missing.begin(), missing.end()); //CRITICAL reverse pero'!
  reverse(missing.begin(), missing.end());
  safety.unlock();
}

void Prefetch::execute() {
  while(1) {
    if(get_signaled()) return;
    while(1) {
      safety.lock();
      if(missing.size() == 0) {
	safety.unlock();
	break;
      } 
      PServer::Item item = missing.back();
      missing.pop_back();

      if(item.priority > 0 && 
	 mt->patches.ram_used > mt->patches.ram_max &&
	 item.priority >= mt->patches.MaxPriority()) {
	safety.unlock();
      	break;
      }
      PatchInfo &info = mt->index[item.patch];
      
      //      cerr << "prefetching: " << item.patch << endl;
      mt->patches.Lookup(item.patch, info.nvert, info.nface, item.priority);

      safety.unlock();
    }
    relax(5);
  }
}
