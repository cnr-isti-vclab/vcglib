#include <GL/glew.h>

#include "queuepserver.h"

using namespace std;
using namespace nxs;
using namespace pt;

 /*
QueuePServer::Data &QueuePServer::Lookup(unsigned int patch, 
			    unsigned short nv, unsigned short nf,
			    float priority,
          vector<QueuePServer::Data> &flush) {
  flush.clear();
  if(index.count(patch)) {    
    return index[patch];    
  } else {
    while(ram_used > ram_max) {
      if(MaxPriority() == 0) break;
      pop_heap(heap.begin(), heap.end());
      Item item = heap.back();      
      heap.pop_back();
      assert(item.priority != 0);
      Data &data = index[item.patch];
      flush.push_back(data);      
      FlushPatch(item.patch, data.patch);
      index.erase(item.patch);
    }
    Item item(patch, priority);
    heap.push_back(item);
    push_heap(heap.begin(), heap.end());
    Data &data = index[patch];
    //    cerr << "Loading: " << patch << endl;
    data.patch = LoadPatch(patch, nv, nf);

   
    return data;
  }
}



bool QueuePServer::IsLoaded(unsigned int patch) {
  return index.count(patch);
}

float QueuePServer::MaxPriority() {
  if(!heap.size()) return 0;
  return heap.front().priority;
}

void QueuePServer::Flush() {
  std::map<unsigned int, Data>::iterator i;
  for(i = index.begin(); i != index.end(); i++) {
    FlushVbo((*i).second);
    FlushPatch((*i).first, (*i).second.patch);
  }
}       */

void QueuePServer::LoadVbo(Data &data) {
//WORKING  if(!vbo_max) return;
  Patch &patch  = *data.patch;
  glGenBuffersARB(1, &data.vbo_element);
  assert(data.vbo_element);
  glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER_ARB, data.vbo_element);
    
  unsigned int size = patch.nf * sizeof(unsigned short);
  if((signature & NXS_FACES) != 0) size *= 3;
    
  glBufferDataARB(GL_ELEMENT_ARRAY_BUFFER_ARB, size, patch.FaceBegin(),
		  GL_STATIC_DRAW_ARB);
  vbo_used += size;

  //TODO fix this when we allow data :p
  size = sizeof(float) * patch.dstart;
    
  glGenBuffersARB(1, &data.vbo_array);
  assert(data.vbo_array);
  glBindBufferARB(GL_ARRAY_BUFFER_ARB, data.vbo_array);
    
  glBufferDataARB(GL_ARRAY_BUFFER_ARB, size, patch.VertBegin(), 
		  GL_STATIC_DRAW_ARB);
    
  vbo_used += size;
}
void QueuePServer::FlushVbo(Data &data) {
//WORKING  if(!vbo_max) return;
  if(!data.vbo_element) return;
  //glDeleteBuffersARB(1, &data.vbo_element);
  //glDeleteBuffersARB(1, &data.vbo_array);
  //data.vbo_element = 0;
  // data.vbo_array = 0;

  Patch &patch  = *data.patch; 
  vbo_used -= patch.nf * sizeof(unsigned short);
  vbo_used -= sizeof(float) * patch.dstart;
}


