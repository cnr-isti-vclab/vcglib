#include <iostream>
#include "pserver.h"

using namespace std;
using namespace nxs;


bool PServer::Create(const std::string &filename, 
		     Signature sig, 
		     unsigned int csize,
		     unsigned int rsize) {
  signature = sig;
  chunk_size = csize;

  ram_max = rsize/chunk_size + 1;
  ram_used = 0;

  return MFile::Create(filename);
}

bool PServer::Load(const std::string &filename, Signature sig, 
		       unsigned int csize, bool readonly, 
		       unsigned int rsize) {

  signature = sig;
  chunk_size = csize;

  ram_max = rsize/chunk_size + 1;
  ram_used = 0;

  return MFile::Load(filename, readonly);
}

void PServer::Close() {  
  Flush();
  MFile::Close();
}

//TODO add error checking.
bool PServer::ReadEntries(FILE *fp) {
  unsigned int n;
  fread(&n, 1, sizeof(int), fp);
  entries.resize(n);
  for(unsigned int i = 0; i < n; i++) {
    fread(&(entries[i].patch_start), 1, sizeof(unsigned int), fp);
    fread(&(entries[i].ram_size),  1, sizeof(unsigned short), fp);
    fread(&(entries[i].disk_size),    1, sizeof(unsigned short), fp);
    entries[i].patch = NULL;
  }
  return true;
}

bool PServer::WriteEntries(FILE *fp) {
  unsigned int n = entries.size();
  fwrite(&n, 1, sizeof(int), fp);
  for(unsigned int i = 0; i < entries.size(); i++) {
    fwrite(&(entries[i].patch_start), 1, sizeof(unsigned int), fp);
    fwrite(&(entries[i].ram_size),  1, sizeof(unsigned short), fp);
    fwrite(&(entries[i].disk_size),    1, sizeof(unsigned short), fp);
  }
  return true;
}

void PServer::AddPatch(unsigned short nvert, unsigned short nface) {
  Entry entry;
  entry.patch_start = 0xffffffff;
  entry.ram_size = Patch::ChunkSize(signature, nvert, nface, chunk_size);
  entry.disk_size = 0xffff;
  entry.patch = NULL;
  entries.push_back(entry);
}

Patch *PServer::LoadPatch(unsigned int idx, 
			  unsigned short nvert, unsigned short nface) {
  
  //  ramlock.rdlock();
  
  assert(idx < entries.size());
  Entry &entry = entries[idx];  
  if(entry.patch) return entry.patch;
  
  char *ram = new char[entry.ram_size * chunk_size];
#ifndef NDEBUG
  if(!ram) {
    cerr << "COuld not allocate ram!\n";
    exit(0);
  }
#endif

  Patch *patch = new Patch(signature, ram, nvert, nface);
      
  if(entry.patch_start != 0xffffffff) { //was allocated.
    assert(entry.disk_size != 0xffff);
	
    SetPosition((int64)entry.patch_start * (int64)chunk_size);
	
    if((signature & NXS_COMPRESSED) == 0) { //not compressed
      ReadBuffer(ram, entry.disk_size * chunk_size);
    } else {
      unsigned char *disk = new unsigned char[entry.disk_size * chunk_size];
      ReadBuffer(disk, entry.disk_size * chunk_size);
      
      patch->Decompress(entry.ram_size * chunk_size, 
			disk, entry.disk_size * chunk_size);
      delete []disk;
    } 
  }
  ram_used += entry.ram_size;  
  entry.patch = patch;  
  return patch;
}

void PServer::FlushPatch(unsigned int id, Patch *patch) {
  //TODO move this into an assert!!!!  
  if(!patch) return;  
  Entry &entry = entries[id];    
//  cerr << "entry: " << (void *)(entry.patch) << " patch: " << (void *)patch << endl;  
  entry.patch = NULL;    

  if(!readonly) { //write back patch
    if((signature & NXS_COMPRESSED)) {
      unsigned int compressed_size;
      char *compressed = patch->Compress(entry.ram_size * chunk_size,
					       compressed_size);
      if(entry.disk_size == 0xffff) {//allocate space 
	assert(entry.patch_start == 0xffffffff);
	entry.disk_size = (unsigned int)((compressed_size-1)/chunk_size) + 1;
	entry.patch_start = (unsigned int)(Length()/chunk_size);
	Redim(Length() + entry.disk_size * chunk_size);
      } else {
	//cerr << "OOOOPSPPPS not supported!" << endl;
	exit(-1);
      }
      SetPosition((int64)entry.patch_start * (int64)chunk_size);
      WriteBuffer(compressed, entry.disk_size * chunk_size);
      delete []compressed;
    } else {
      if(entry.disk_size == 0xffff) {
	entry.disk_size = entry.ram_size;
	entry.patch_start = (unsigned int)(Length()/chunk_size);
	Redim(Length() + entry.disk_size * chunk_size);
      }
      SetPosition((int64)entry.patch_start * (int64)chunk_size);
      WriteBuffer(patch->start, entry.disk_size * chunk_size);
    }
  }

  delete [](patch->start);
  delete patch;  
  ram_used -= entry.ram_size;      
}

void PServer::MaxRamBuffer(unsigned int r_buffer) {  
  ram_max = (unsigned int)(r_buffer/chunk_size) + 1;  
}
