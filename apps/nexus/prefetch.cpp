#include <set>

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
  mt->todraw.clear();
  
  unsigned int notloaded = 0;
  set<unsigned int> tmp;
  //std::map<unsigned int, float> tmp;
  vector<QueuePServer::Data> flush;
  for(unsigned int i = 0; i < selected.size(); i++) {
    unsigned int patch = selected[i];
    tmp.insert(patch);
    //tmp[patch] = 0.0f;
    if(!mt->patches.entries[patch].patch) {   
      //cerr << "miss: " << patch << endl;
      load.post(patch);
      notloaded++;  
    } else {
       PatchInfo &info = mt->index[patch];             
       QueuePServer::Data &data = mt->patches.Lookup(patch, info.nvert, info.nface, 0.0f, flush);
       if(flush.size() != 0) {
         cerr << "Flushing!\n";
          exit(0);
       }
       mt->todraw.push_back(&data);      
    }    
    //missing.push_back(PServer::Item(patch, 0.0f));
  }
  if(notloaded)
    cerr << "Patches to load: " << notloaded << endl;

  for(unsigned int i = 0; i < visited.size(); i++) {
    PServer::Item &item = visited[i];
    if(tmp.count(item.patch)) continue;
//    if(mt->patches.entries[item.patch].patch)    
//      tmp[item.patch] = item.priority;    
    if(item.priority != 0.0f)
      missing.push_back(item);
  }

  QueuePServer &ps = mt->patches;
  for(unsigned int i = 0; i < ps.heap.size(); i++) {
    PServer::Item &item = ps.heap[i];
    if(tmp.count(item.patch)) 
      item.priority = 0;
    else {
      if(item.priority == 0)
        item.priority = 1;
      item.priority *= 1.1;
    }

    /*if(tmp.count(item.patch)) 
      item.priority = tmp[item.patch];
    else
      item.priority = 1e30;*/
  }
  make_heap(ps.heap.begin(), ps.heap.end());

  sort(missing.begin(), missing.end()); //CRITICAL reverse pero'!
  reverse(missing.begin(), missing.end());    
  load.post(0xffffffff);
 
  safety.unlock();
}

void Prefetch::execute() {      
    while(1) {
      if(get_signaled()) return;               
      vector<QueuePServer::Data> flush;      
      
      if(load.get_count() || missing.size() == 0) {
        pt::message *msg = load.getmessage();
        if(msg->id != 0xffffffff) {          
          safety.lock();
          PatchInfo &info = mt->index[msg->id];                   

          //posting draw message
          QueuePServer::Data &data = mt->patches.Lookup(msg->id, info.nvert, info.nface, 0.0f, flush);
          pt::message *msg = new pt::message(QueuePServer::DRAW, (unsigned int)&data);
          msg->result = msg->id;
          draw.post(msg);

          //p;osting flush messages
          for(unsigned int i = 0; i < flush.size(); i++) {
            QueuePServer::Data *data = new QueuePServer::Data;
            *data = flush[i];
            draw.post(QueuePServer::FLUSH, (unsigned int)data);
          }
          safety.unlock();
        }
        delete msg;        
      } else {             
        safety.lock();
        if(missing.size() != 0) {          
          PServer::Item item = missing.back();
          missing.pop_back();

          if(item.priority > mt->patches.MaxPriority()) {	        
            missing.clear();
            
          } else {        
            PatchInfo &info = mt->index[item.patch];      
            //cerr << "prefetching: " << item.patch << endl;
            mt->patches.Lookup(item.patch, info.nvert, info.nface, item.priority, flush);
            for(unsigned int i = 0; i < flush.size(); i++) {
              QueuePServer::Data  *data = new QueuePServer::Data;
              *data = flush[i];
              draw.post(QueuePServer::FLUSH, (unsigned int)data);
            }
          }
        }
        safety.unlock();
      }      
    }    
}
