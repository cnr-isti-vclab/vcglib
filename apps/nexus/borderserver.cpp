#include "borderserver.h"
#include <iostream>

using namespace std;
using namespace nxs;

bool BorderServer::Create(const string &file) {
  ram_used = 0;
  return IndexFile<BorderEntry>::Create(file, 2 * sizeof(Link));
}

bool BorderServer::Load(const string &file, bool rdonly) {
  ram_used = 0;
  cerr << "Loading...\n";
  return IndexFile<BorderEntry>::Load(file, rdonly);
}

void BorderServer::Close() {
  if(!Opened()) return;
  Flush();
  IndexFile<BorderEntry>::Close();
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

void BorderServer::AddBorder(unsigned short nbord, unsigned int used) {
  BorderEntry entry;
  assert((Length() % sizeof(Link)) == 0);

  entry.start = Length()/ sizeof(Link);
  entry.size = nbord;
  entry.used = used;
  entry.links = NULL;
  push_back(entry);
  Redim(entry.start * sizeof(Link) + nbord * sizeof(Link));
}

Border BorderServer::GetBorder(unsigned int border, bool flush) { 
  BorderEntry &entry = operator[](border);
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
    assert(!entry.links);
    if(entry.size != 0)
      entry.links = GetRegion(entry.start, entry.size);        
    pqueue.push_front(border);    
    index[border] = pqueue.begin();   
    ram_used += entry.size;
  }                        
  return Border(entry.links, entry.used, entry.size);
}
//TODO Change when remving borderentry class.
bool BorderServer::ResizeBorder(unsigned int border, unsigned int nbord) {
  assert(nbord < 65500);
  assert(border < size());    
  GetBorder(border);
  BorderEntry &entry = operator[](border);
  if(nbord > entry.size) {
    int capacity = nbord;
    if(capacity < entry.size*2) 
      capacity = entry.size * 2;
    if(capacity > 65500) 
      capacity = 65500;
    unsigned int newstart = Length()/sizeof(Link);
    Redim((newstart + capacity) * sizeof(Link));
    Link *dst = GetRegion(newstart, capacity);
    if(entry.used > 0) {
      Link *src = GetRegion(entry.start, entry.size);      
      memcpy(dst, src, entry.used * sizeof(Link));
    }
    entry.links = dst;
    entry.start = newstart;
    entry.size = capacity;
    entry.used = nbord;
    return true;
  }
  entry.used = nbord;
  return false;
}

void BorderServer::FlushBorder(unsigned int border) {
  BorderEntry &entry = operator[](border);
  //assert(entry.links);
  if(entry.size && !MFile::IsReadOnly()) { //write back patch
    MFile::SetPosition((int64)entry.start * sizeof(Link));
    MFile::WriteBuffer(entry.links, entry.used * sizeof(Link));
  }
  if(entry.size)
    delete [](entry.links);
  entry.links = NULL;    
  ram_used -= entry.size;
}

Link *BorderServer::GetRegion(unsigned int start, unsigned int size) {
  assert(size > 0);
  SetPosition(start * sizeof(Link));
  Link *buf = new Link[size];
  assert(buf);
  ReadBuffer(buf, size * sizeof(Link));
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
  cerr << "Offset: " << offset << endl;
  return true;
}

void BorderServer::SaveHeader() {
  unsigned int magic = 0x3042584e; // NXB0
  WriteBuffer(&magic, sizeof(unsigned int));
  WriteBuffer(&offset, sizeof(int64));
}
