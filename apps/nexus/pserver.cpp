#include <iostream>
#include "pserver.h"

using namespace std;
using namespace nxs;

bool PServer::Create(const std::string &_filename, 
		     Signature _signature, 
		     unsigned int _chunk_size) {
  signature = _signature;
  chunk_size = _chunk_size;

  ram_used = 0;

  return MFile::Create(filename);
}

bool PServer::Load(const std::string &filename, Signature _signature, 
		   bool _readonly, unsigned int _chunk_size) {

  signature = _signature;
  chunk_size = _chunk_size;

  ram_used = 0;

  return MFile::Load(filename, _readonly);
}

void PServer::Close() {  
  Flush();
  MFile::Close();
}

//TODO add error checking.
bool PServer::ReadEntries(FILE *fp) {
  unsigned int n;
  fread(&n, 1, sizeof(int), fp);
  resize(n);
  for(unsigned int i = 0; i < n; i++) {
    Entry &entry = operator[](i);
    fread(&entry, sizeof(Entry), 1, fp);
    entry.patch = NULL;
  }
  return true;
}

bool PServer::WriteEntries(FILE *fp) {
  unsigned int n = size();
  fwrite(&n, 1, sizeof(int), fp);
  for(unsigned int i = 0; i < size(); i++) {
    Entry &entry = operator[](i);
    fwrite(&entry, sizeof(Entry), 1, fp);
  }
  return true;
}

void PServer::AddPatch(unsigned short nvert, unsigned short nface) {
  Entry entry;
  entry.patch_start = 0xffffffff;
  entry.ram_size = Patch::ChunkSize(signature, nvert, nface, chunk_size);
  entry.disk_size = 0xffff;
  entry.nvert = nvert;
  entry.nface = nface;
  //sphere and error undefined.
  entry.patch = NULL;
  entry.vbo_array = 0;
  entry.vbo_element = 0;
  
  push_back(entry);
}

Patch *PServer::LoadPatch(unsigned int idx) {
  assert(idx < size());
  Entry &entry = operator[](idx);
  if(entry.patch) return entry.patch;
  
  char *ram = new char[entry.ram_size * chunk_size];
#ifndef NDEBUG
  if(!ram) {
    cerr << "COuld not allocate ram!\n";
    exit(0);
  }
#endif

  Patch *patch = new Patch(signature, ram, entry.nvert, entry.nface);
  
  if(entry.patch_start != 0xffffffff) { //was allocated.
    assert(entry.disk_size != 0xffff);
    
    MFile::SetPosition((int64)entry.patch_start * (int64)chunk_size);
	
    if((signature & NXS_COMPRESSED) == 0) { //not compressed
      MFile::ReadBuffer(ram, entry.disk_size * chunk_size);
    } else {
      unsigned char *disk = new unsigned char[entry.disk_size * chunk_size];
      MFile::ReadBuffer(disk, entry.disk_size * chunk_size);
      
      patch->Decompress(entry.ram_size * chunk_size, 
			disk, entry.disk_size * chunk_size);
      delete []disk;
    } 
  }
  ram_used += entry.ram_size;  
  entry.patch = patch;  
  return patch;
}

void PServer::FlushPatch(unsigned int id) {
  Entry &entry = operator[](id);    
  //TODO move this into an assert!!!!  
  if(!entry.patch) return;


  if(!MFile::IsReadOnly()) { //write back patch
    if((signature & NXS_COMPRESSED)) {
      unsigned int compressed_size;
      char *compressed = entry.patch->Compress(entry.ram_size * chunk_size,
					       compressed_size);
      if(entry.disk_size == 0xffff) {//allocate space 
	assert(entry.patch_start == 0xffffffff);
	entry.disk_size = (unsigned int)((compressed_size-1)/chunk_size) + 1;
	entry.patch_start = (unsigned int)( MFile::Length()/chunk_size);
	MFile::Redim(MFile::Length() + entry.disk_size * chunk_size);
      } else {
	//cerr << "OOOOPSPPPS not supported!" << endl;
	exit(-1);
      }
      MFile::SetPosition((int64)entry.patch_start * (int64)chunk_size);
      MFile::WriteBuffer(compressed, entry.disk_size * chunk_size);
      delete []compressed;
    } else {
      if(entry.disk_size == 0xffff) {
	entry.disk_size = entry.ram_size;
	entry.patch_start = (unsigned int)(MFile::Length()/chunk_size);
	MFile::Redim(MFile::Length() + entry.disk_size * chunk_size);
      }
      MFile::SetPosition((int64)entry.patch_start * (int64)chunk_size);
      MFile::WriteBuffer(entry.patch->start, entry.disk_size * chunk_size);
    }
  }

  delete [](entry.patch->start);
  delete entry.patch;  
  entry.patch = NULL;    
  ram_used -= entry.ram_size;      
}
  
bool PServer::IsLoaded(unsigned int patch) { 
  return operator[](patch).patch != NULL;
}

Entry &PServer::Lookup(unsigned int patch, std::vector<unsigned int> &flush) {
  Entry &entry = operator[](patch);
  if(index.count(patch)) {
    list<unsigned int>::iterator &i = index[patch];
    pqueue.erase(i);
    pqueue.push_front(patch);
    assert(entry.patch);
  } else {
    while(ram_used > ram_max) {    
      unsigned int to_flush = pqueue.back();
      pqueue.pop_back();
      //TODO i should not flush current extraction!
      //CRITICAL actually works just if ram_max is big enough.
      index.erase(to_flush);        
      FlushPatch(to_flush);
      flush.push_back(to_flush);
    }
    assert(!entry.patch);
    entry.patch = LoadPatch(patch);
    pqueue.push_front(patch);
    list<unsigned int>::iterator i = pqueue.begin();
    index[patch] = i;      
  }                        
  return entry;
}

void PServer::Flush() {
  std::map<unsigned int, list<unsigned int>::iterator>::iterator i;
  for(i = index.begin(); i != index.end(); i++) {
    unsigned int patch = *((*i).second);
    Item &item = *((*i).second);
    FlushVbo(patch);
    FlushPatch(patch);
  }
  
  pqueue.clear();
  index.clear();
}

void PServer::LoadVbo(unsigned int npatch) {
//WORKING  if(!vbo_max) return;
  Entryg &entry = operator[](npatch);
  Patch &patch  = *entry.patch;
  glGenBuffersARB(1, &entry.vbo_element);
  assert(entry.vbo_element);
  glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER_ARB, entry.vbo_element);
    
  unsigned int size = patch.nf * sizeof(unsigned short);
  if((signature & NXS_FACES) != 0) size *= 3;
    
  glBufferDataARB(GL_ELEMENT_ARRAY_BUFFER_ARB, size, patch.FaceBegin(),
		  GL_STATIC_DRAW_ARB);
  vbo_used += size;

  //TODO fix this when we allow data :p
  size = sizeof(float) * patch.dstart;
    
  glGenBuffersARB(1, &entry.vbo_array);
  assert(entry.vbo_array);
  glBindBufferARB(GL_ARRAY_BUFFER_ARB, entry.vbo_array);
    
  glBufferDataARB(GL_ARRAY_BUFFER_ARB, size, patch.VertBegin(), 
		  GL_STATIC_DRAW_ARB);
    
  vbo_used += size;
}
void PServer::FlushVbo(unsigned int npatch) {
  Entryg &entry = operator[](npatch);
  if(!entry.vbo_element) return;
  glDeleteBuffersARB(1, &entry.vbo_element);
  glDeleteBuffersARB(1, &entry.vbo_array);
  entry.vbo_element = 0;
  entry.vbo_array = 0;

  Patch &patch  = *entry.patch; 
  vbo_used -= patch.nf * sizeof(unsigned short);
  vbo_used -= sizeof(float) * patch.dstart;
}


