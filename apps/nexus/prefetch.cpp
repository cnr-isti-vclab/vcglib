#include "prefetch.h"

using namespace std;
using namespace nxs;


void Prefetch::init(NexusMt *m, std::vector<unsigned int> &selected,
		    std::vector<PServer::Item> &visited) {
  safety.lock();
  mt = m;
  missing.clear();
  
  std::map<unsigned int, float> tmp;
  for(unsigned int i = 0; i < selected.size(); i++) {
    unsigned int patch = selected[i];
    tmp[patch] = 0.0f;
    missing.push_back(PServer::Item(patch, 0.0f));
  }
  for(unsigned int i = 0; i < visited.size(); i++) {
    PServer::Item &item = visited[i];
    if(tmp.count(item.patch)) continue;
    if(mt->patches.IsLoaded(item.patch))
      tmp[item.patch] = item.priority;
    else
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
      if(missing.size() == 0) break;
      PServer::Item item = missing.back();
      missing.pop_back();
      PatchInfo &info = mt->index[item.patch];
      mt->patches.Lookup(item.patch, info.nvert, info.nface, item.priority);
      safety.unlock();
    }
    relax(50);
  }
}
