#include "patchserver.h"

#include <iostream>
#include <algorithm>

using namespace std;
using namespace nxs;

//TODO support compression!


bool PatchServer::Create(const std::string &filename, 
			 Signature sig, 
			 unsigned int csize,
			 unsigned int rsize) {
  signature = sig;
  chunk_size = csize;
  frame = 0;
  ram_size = rsize;
  ram_used = 0;
  lru.clear();
  return File::Create(filename);
}

bool PatchServer::Load(const std::string &filename, Signature sig, 
		       unsigned int csize, bool readonly, 
		       unsigned int rsize) {
  signature = sig;
  chunk_size = csize;
  ram_size = rsize;
  frame = 0;
  ram_used = 0;
  lru.clear();
  return File::Load(filename, readonly);
}

void PatchServer::Close() {  
  FlushAll();
  File::Close();
}

//TODO add error checking.
bool PatchServer::ReadEntries(FILE *fp) {
  unsigned int n;
  fread(&n, 1, sizeof(int), fp);
  patches.resize(n);
  for(unsigned int i = 0; i < n; i++) {
    patches[i].patch = NULL;
    fread(&(patches[i].patch_start), 1, sizeof(unsigned int), fp);
    fread(&(patches[i].patch_size),  1, sizeof(unsigned short), fp);
    fread(&(patches[i].ram_used),    1, sizeof(unsigned short), fp);
    patches[i].lru_pos = 0xffffffff;
  }
  return true;
}

bool PatchServer::WriteEntries(FILE *fp) {
  unsigned int n = patches.size();
  fwrite(&n, 1, sizeof(int), fp);
  for(unsigned int i = 0; i < patches.size(); i++) {
    fwrite(&(patches[i].patch_start), 1, sizeof(unsigned int), fp);
    fwrite(&(patches[i].patch_size),  1, sizeof(unsigned short), fp);
    fwrite(&(patches[i].ram_used),    1, sizeof(unsigned short), fp);
  }
  return true;
}

void PatchServer::AddPatch(unsigned short nvert, unsigned short nface) {
  PatchEntry entry;
  entry.patch = NULL;
  entry.patch_start = Length()/chunk_size;
  entry.patch_size = Patch::ChunkSize(signature, nvert, nface, chunk_size);
  entry.ram_used = entry.patch_size;
  entry.lru_pos = 0xffffffff;
  patches.push_back(entry);

  Redim(Length() + entry.patch_size * chunk_size);
}

Patch &PatchServer::GetPatch(unsigned int idx, 
			     unsigned short nvert, unsigned short nface,
			     bool flush) {

  assert(idx < patches.size());
  PatchEntry &entry = patches[idx];

  if(entry.patch) { 
    assert(entry.lru_pos < lru.size());
    assert(lru[entry.lru_pos].patch == idx);
    lru[entry.lru_pos].frame = frame++;
  } else {
    SetPosition(entry.patch_start * chunk_size);
    
    assert(entry.patch_size != 0);
    char *start = new char[entry.patch_size * chunk_size];
    ReadBuffer(start, entry.patch_size * chunk_size);
   
    entry.patch = new Patch(signature, start, nvert, nface);
    entry.lru_pos = lru.size();
    lru.push_back(PTime(idx, frame++));
    ram_used += entry.ram_used;
  }
  
  //avoid frame overflow!
  if(frame > (1<<30)) {
    cerr << "oVERFLOW! (nothing dangerous... just warning." << endl;;
    for(unsigned int i = 0; i < lru.size(); i++) {
      if(lru[i].frame < (1<<29)) lru[i].frame = 0;
      else lru[i].frame -= (1<<29);
    }
    make_heap(lru.begin(), lru.end());
    for(unsigned int i = 0; i < lru.size(); i++)
      patches[lru[i].patch].lru_pos = i;
  }

  if(flush && ram_used > ram_size * 1.1)
    Flush();

  return *(entry.patch); 
}

void PatchServer::Flush() {
  cerr << "FLUSHING\n\n\n n";
  cerr << "ram_size: " << ram_size << endl;
  cerr << "ram_used: " << ram_used << endl;
  make_heap(lru.begin(), lru.end());
  while(ram_used > ram_size) {
    pop_heap(lru.begin(), lru.end());
    PTime &ptime = lru.back();

    Flush(ptime.patch);

    lru.pop_back();
  }
  make_heap(lru.begin(), lru.end());
  for(unsigned int i = 0; i < lru.size(); i++)
    patches[lru[i].patch].lru_pos = i;
}

void PatchServer::FlushAll() {
  for(unsigned int i = 0; i < lru.size(); i++) {
    PTime &ptime = lru[i];
    Flush(ptime.patch);
  }
  assert(ram_used == 0);
  lru.clear();
}


void PatchServer::Flush(unsigned int patch) {

  PatchEntry &entry = patches[patch];
  assert(entry.patch);

  if(!readonly) { //write back patch
    SetPosition(entry.patch_start * chunk_size);
    WriteBuffer(entry.patch->start, entry.patch_size * chunk_size);
  }
  delete [](entry.patch->start);
  delete entry.patch;
  entry.patch = NULL;
  entry.lru_pos = 0xffffffff;
  ram_used -= entry.ram_used;
}
