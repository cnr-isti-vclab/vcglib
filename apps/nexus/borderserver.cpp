#include "borderserver.h"
#include <iostream>

using namespace std;
using namespace nxs;

bool BorderServer::Create(const string &file) {
  ram_used = 0;
  return IndexFile<Border>::Create(file, 2 * sizeof(Link));
}

bool BorderServer::Load(const string &file, bool rdonly) {
  ram_used = 0;  
  bool success = IndexFile<Border>::Load(file, rdonly);
  if(!success) return false;
  for(unsigned int i = 0; i < size(); i++) 
    operator[](i).links = NULL;
  return true;
}

void BorderServer::Close() {
  if(!Opened()) return;
  Flush();
  IndexFile<Border>::Close();
}

void BorderServer::Flush() {
  std::map<unsigned int, list<unsigned int>::iterator>::iterator i;
  for(i = index.begin(); i != index.end(); i++) {
    unsigned int patch = (*i).first;
    FlushBorder(patch);
  }
  pqueue.clear();
  index.clear();
}

void BorderServer::AddBorder(unsigned short _size, unsigned int used) {
  Border entry;
  assert((Length() % sizeof(Link)) == 0);

  entry.start = Length()/ sizeof(Link);
  entry.size = _size;
  entry.used = used;
  entry.links = NULL;
  push_back(entry);
  Redim((int64)entry.start * (int64)sizeof(Link) + (int64)_size * (int64)sizeof(Link));
}

Border &BorderServer::GetBorder(unsigned int border, bool flush) { 
  Border &entry = operator[](border);
   //assert(entry.size != 0);
  if(index.count(border)) {
    //assert(entry.links);
    list<unsigned int>::iterator i = index[border];
    pqueue.erase(i);
    pqueue.push_front(border);
    index[border] = pqueue.begin();
  } else {
    while(flush && ram_used > ram_max) { 
      assert(pqueue.size());
      unsigned int to_flush = pqueue.back();
      pqueue.pop_back();
      index.erase(to_flush);        
      FlushBorder(to_flush);
    }        
    entry.links = GetRegion(entry.start, entry.size);
    pqueue.push_front(border);    
    index[border] = pqueue.begin();   
    ram_used += entry.size;
  }                        
  return entry;
}
//TODO Change when remving borderentry class.
void BorderServer::ResizeBorder(unsigned int border, unsigned int used) {  
  assert(border < size());      
  Border &entry = GetBorder(border);
  if(used <= entry.size) {
    entry.used = used;
    return;
  }
  
  int capacity = used;  
  if(capacity < entry.size * 2) 
    capacity = entry.size * 2;
  
  unsigned int newstart = Length()/sizeof(Link);  
  Redim((int64)(newstart + capacity) * (int64)sizeof(Link));
  Link *newlinks = new Link[capacity];
  if(entry.used > 0) {
    assert(entry.links);
    memcpy(newlinks, entry.links, entry.used * sizeof(Link));
    delete []entry.links;
    entry.links = NULL;
  }
  assert(entry.links == NULL);
  entry.links = newlinks;  
  entry.start = newstart;
  entry.size = capacity;
  entry.used = used;  
}

void BorderServer::FlushBorder(unsigned int border) {
  Border &entry = operator[](border);
  //assert(entry.links);
  if(entry.size && !MFile::IsReadOnly()) { //write back patch
    MFile::SetPosition((int64)entry.start * (int64)sizeof(Link));
    MFile::WriteBuffer(entry.links, entry.used * (int64)sizeof(Link));
  }
  if(entry.links)
    delete [](entry.links);
  entry.links = NULL;    
  ram_used -= entry.size;
}

Link *BorderServer::GetRegion(unsigned int start, unsigned int size) {
  if(size == 0) return NULL;  
  SetPosition((int64)start * (int64)sizeof(Link));
  Link *buf = new Link[size];
  assert(buf);
  ReadBuffer(buf, (int64)size * (int64)sizeof(Link));
  return buf;
}

bool BorderServer::LoadHeader() {
  unsigned int magic;
  ReadBuffer(&magic, sizeof(unsigned int));
  if(magic != 0x3042584e) { //NXB0
    cerr << "Invalid magic. Not a nxs file\n";
    return false;
  } 
  ReadBuffer(&offset, sizeof(int64));  
  return true;
}

void BorderServer::SaveHeader() {
  unsigned int magic = 0x3042584e; // NXB0
  WriteBuffer(&magic, sizeof(unsigned int));
  WriteBuffer(&offset, sizeof(int64));
}
