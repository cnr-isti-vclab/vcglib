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
  
/*void Preload::init(NexusMt *m, std::vector<unsigned int> &selected,
		    std::vector<PServer::Item> &visited) {

  //  cerr << "Init\n";

  safety.lock();
  mt = m;
  missing.clear();
  mt->todraw.clear();
  float loaded = 0;  
  unsigned int notloaded = 0;
  set<unsigned int> tmp;
  //std::map<unsigned int, float> tmp;
  vector<QueuePServer::Data> flush;
  for(unsigned int i = 0; i < selected.size(); i++) {
    unsigned int patch = selected[i];
    tmp.insert(patch);   
    
    if(!mt->patches.entries[patch].patch) {         
      PServer::Entry &entry = mt->patches.entries[patch];
      load.post(patch);
      loaded += entry.disk_size;
    } else {
       PatchInfo &info = mt->index[patch];             
       //WORKING QueuePServer::Data &data = mt->patches.Lookup(patch, info.nvert, info.nface, 0.0f, flush);
       QueuePServer::Data &data = mt->patches.Lookup(patch, info.nvert, info.nface, flush);
       for(unsigned int i = 0; i < flush.size(); i++) {
          QueuePServer::Data  *data = new QueuePServer::Data;
          *data = flush[i];
          draw.post(FLUSH, (unsigned int)data);
       }
       // WORKING if(flush.size() != 0) {
       //  cerr << "Flushing!\n";
       //   exit(0);
       //} 
       mt->todraw.push_back(&data);      
    }    
    //missing.push_back(PServer::Item(patch, 0.0f));
  }
  loading = 0.2 * loaded + 0.8 * loading;

  for(unsigned int i = 0; i < visited.size(); i++) {
    PServer::Item &item = visited[i];
    if(tmp.count(item.patch)) continue;
//    if(mt->patches.entries[item.patch].patch)    
//      tmp[item.patch] = item.priority;    
    if(item.priority != 0.0f)
      missing.push_back(item);
      }*/

/*  WORKING QueuePServer &ps = mt->patches;
  for(unsigned int i = 0; i < ps.heap.size(); i++) {
    PServer::Item &item = ps.heap[i];
    if(tmp.count(item.patch)) 
      item.priority = 0;
    else {
      if(item.priority == 0)
        item.priority = 1;
      item.priority *= 1.1;
    }    
  }
  make_heap(ps.heap.begin(), ps.heap.end());*/

/*  sort(missing.begin(), missing.end()); //CRITICAL reverse pero'!
  reverse(missing.begin(), missing.end());    
  load.post(0xffffffff);
 
  safety.unlock();
}

void Preload::execute() {      

  float preload;
  prefetching = 0;
  loading = 0;
    while(1) {
      if(get_signaled()) return;               
      vector<QueuePServer::Data> flush;      
      
      if(load.get_count() || missing.size() == 0) {  
        preload = 0;
        pt::message *msg = load.getmessage();
        if(msg->id != 0xffffffff) {          
          safety.lock();
          PatchInfo &info = mt->index[msg->id];                   
          PServer::Entry &entry = mt->patches.entries[msg->id];
          loading += entry.disk_size;
          //posting draw message
          //WORKING QueuePServer::Data &data = mt->patches.Lookup(msg->id, info.nvert, info.nface, 0.0f, flush);
          QueuePServer::Data &data = mt->patches.Lookup(msg->id, info.nvert, info.nface, flush);
          pt::message *msg = new pt::message(DRAW, (unsigned int)&data);
          msg->result = msg->id;
          draw.post(msg);

          //p;osting flush messages
          for(unsigned int i = 0; i < flush.size(); i++) {
            QueuePServer::Data *data = new QueuePServer::Data;
            *data = flush[i];
            draw.post(FLUSH, (unsigned int)data);
          }
          safety.unlock();
        } else {
          prefetching = 0.2 * preload + 0.8 * prefetching;
        }
        delete msg;        
      } else {             
        safety.lock();
        if(missing.size() != 0) {          
          PServer::Item item = missing.back();
          missing.pop_back();

/*Working          if(item.priority > mt->patches.MaxPriority()) {	        
            missing.clear();
            
          } else {        */
/*PatchInfo &info = mt->index[item.patch];      
            //cerr << "prefetching: " << item.patch << endl;
            //WORKING mt->patches.Lookup(item.patch, info.nvert, info.nface, item.priority, flush);
            if(!mt->patches.entries[item.patch].patch) {
              PServer::Entry &entry = mt->patches.entries[item.patch];
              preload += entry.disk_size;
          
              mt->patches.Lookup(item.patch, info.nvert, info.nface, flush);              
              for(unsigned int i = 0; i < flush.size(); i++) {
                QueuePServer::Data  *data = new QueuePServer::Data;
                *data = flush[i];
                draw.post(FLUSH, (unsigned int)data);
              }
            }
        //  }
        }
        safety.unlock();
      }      
    }    
}*/
