#include "patchserver.h"
#include <iostream>
#include <algorithm>

#include <GL/glew.h>

using namespace std;
using namespace nxs;


bool PatchServer::Create(const std::string &filename, 
			 Signature sig, 
			 unsigned int csize,
			 unsigned int rsize) {
  signature = sig;
  chunk_size = csize;

  ram_size = rsize/chunk_size + 1;
  ram_used = 0;
  vbo_size = 0;
  vbo_used = 0;

  frame = 0;
  ram_readed = 0;
  ram_flushed = 0;

  lru.clear();
  return MFile::Create(filename);
}

bool PatchServer::Load(const std::string &filename, Signature sig, 
		       unsigned int csize, bool readonly, 
		       unsigned int rsize) {

  signature = sig;
  chunk_size = csize;

  ram_size = rsize/chunk_size + 1;
  ram_used = 0;
  vbo_size = 0;
  vbo_used = 0;

  frame = 0;
  ram_readed = 0;
  ram_flushed = 0;

  lru.clear();
  return MFile::Load(filename, readonly);
}

void PatchServer::Close() {  
  FlushAll();
  MFile::Close();
}

//TODO add error checking.
bool PatchServer::ReadEntries(FILE *fp) {
  unsigned int n;
  fread(&n, 1, sizeof(int), fp);
  patches.resize(n);
  for(unsigned int i = 0; i < n; i++) {
    fread(&(patches[i].patch_start), 1, sizeof(unsigned int), fp);
    fread(&(patches[i].ram_size),  1, sizeof(unsigned short), fp);
    fread(&(patches[i].disk_size),    1, sizeof(unsigned short), fp);
    patches[i].lru_pos = 0xffffffff;
  }
  return true;
}

bool PatchServer::WriteEntries(FILE *fp) {
  unsigned int n = patches.size();
  fwrite(&n, 1, sizeof(int), fp);
  for(unsigned int i = 0; i < patches.size(); i++) {
    fwrite(&(patches[i].patch_start), 1, sizeof(unsigned int), fp);
    fwrite(&(patches[i].ram_size),  1, sizeof(unsigned short), fp);
    fwrite(&(patches[i].disk_size),    1, sizeof(unsigned short), fp);
  }
  return true;
}

void PatchServer::AddPatch(unsigned short nvert, unsigned short nface) {
  PatchEntry entry;

  entry.patch_start = 0xffffffff;
  entry.ram_size = Patch::ChunkSize(signature, nvert, nface, chunk_size);
  entry.disk_size = 0xffff;
  entry.lru_pos = 0xffffffff;
  patches.push_back(entry);
}

Patch &PatchServer::GetPatch(unsigned int idx, 
			     unsigned short nvert, unsigned short nface,
			     bool flush) {
  assert(idx < patches.size());
  PatchEntry &entry = patches[idx];

  if(entry.lru_pos == 0xffffffff) {  //not on buffer
    if(flush) Flush();
    PTime nptime(idx);

    char *ram = new char[entry.ram_size * chunk_size];
    nptime.patch = new Patch(signature, ram, nvert, nface);
    
    if(entry.patch_start != 0xffffffff) { //was allocated.
      assert(entry.disk_size != 0xffff);

      SetPosition(entry.patch_start * chunk_size);
    
      if((signature & NXS_COMPRESSED) == 0) { //not compressed
	ReadBuffer(ram, entry.disk_size * chunk_size);
      } else {

	unsigned char *disk = new unsigned char[entry.disk_size * chunk_size];
	ReadBuffer(disk, entry.disk_size * chunk_size);

	nptime.patch->Decompress(entry.ram_size * chunk_size, 
				 disk, entry.disk_size * chunk_size);
	delete []disk;
      } 
    }
    
    entry.lru_pos = lru.size();
    lru.push_back(nptime);
    ram_used += entry.ram_size;
    ram_readed += entry.ram_size;
  }

  PTime &ptime = lru[entry.lru_pos];
  ptime.frame = frame++;
  
  //avoid frame overflow!
  if(frame > (1<<30)) {
    cerr << "oVERFLOW! (nothing dangerous... just warning." << endl;;
    for(unsigned int i = 0; i < lru.size(); i++) {
      if(lru[i].frame < (1<<29)) lru[i].frame = 0;
      else lru[i].frame -= (1<<29);
    }
    make_heap(lru.begin(), lru.end());
    for(unsigned int i = 0; i < lru.size(); i++)
      patches[lru[i].npatch].lru_pos = i;
  }
  return *(ptime.patch); 
}


void PatchServer::GetVbo(unsigned int p, 
			 unsigned int &element, unsigned int &array) {
  PatchEntry &entry = patches[p];
  assert(entry.lru_pos != 0xffffffff);
  PTime &ptime = lru[entry.lru_pos];
  if(!ptime.vbo_element) {
    //TODO  cerr << "Adding vbo: " << p << endl;
    assert(ptime.patch);
    
    Patch &patch = *ptime.patch;
    
    glGenBuffersARB(1, &ptime.vbo_element);
    assert(ptime.vbo_element);
    glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER_ARB, ptime.vbo_element);
    
    unsigned int size = patch.nf * sizeof(unsigned short);
    if((signature & NXS_FACES) != 0) size *= 3;
    
    glBufferDataARB(GL_ELEMENT_ARRAY_BUFFER_ARB, size, patch.FaceBegin(),
		    GL_STATIC_DRAW_ARB);
    
    //TODO fix this when we allow data :p
    size = sizeof(float) * patch.dstart;
    
    glGenBuffersARB(1, &ptime.vbo_array);
    glBindBufferARB(GL_ARRAY_BUFFER_ARB, ptime.vbo_array);
    
    glBufferDataARB(GL_ARRAY_BUFFER_ARB, size, patch.VertBegin(), 
		    GL_STATIC_DRAW_ARB);
    
    vbo_used += patches[p].ram_size;
  }

  element = ptime.vbo_element;
  array = ptime.vbo_array;
}


void PatchServer::Flush() {

  if(ram_used < ram_size * 1.1) return;

  make_heap(lru.begin(), lru.end());
  for(unsigned int i = 0; i < lru.size(); i++)
    patches[lru[i].npatch].lru_pos = i;

  while(ram_used > ram_size) {
    pop_heap(lru.begin(), lru.end());
    PTime &ptime = lru.back();
    Flush(ptime);        
    lru.pop_back();
  }
  for(unsigned int i = 0; i < lru.size(); i++)
    patches[lru[i].npatch].lru_pos = i;
}

void PatchServer::FlushAll() {
  for(unsigned int i = 0; i < lru.size(); i++) {
    PTime &ptime = lru[i];
    Flush(ptime);
  }
  assert(ram_used == 0);
  lru.clear();
}

void PatchServer::Flush(PTime &ptime) {
  PatchEntry &entry = patches[ptime.npatch];  
  assert(ptime.patch);
  if(!readonly) { //write back patch


    if((signature & NXS_COMPRESSED)) {
      unsigned int compressed_size;
      char *compressed = ptime.patch->Compress(entry.ram_size * chunk_size,
					       compressed_size);
      if(entry.disk_size == 0xffff) {//allocate space 
	assert(entry.patch_start == 0xffffffff);
	entry.disk_size = (unsigned int)((compressed_size-1)/chunk_size) + 1;
	entry.patch_start = (unsigned int)(Length()/chunk_size);
	Redim(Length() + entry.disk_size * chunk_size);
      } else {
	cerr << "OOOOPSPPPS not supported!" << endl;
	exit(-1);
      }
      SetPosition(entry.patch_start * chunk_size);
      WriteBuffer(compressed, entry.disk_size * chunk_size);
      delete []compressed;
    } else {
      if(entry.disk_size == 0xffff) {
	entry.disk_size = entry.ram_size;
	entry.patch_start = (unsigned int)(Length()/chunk_size);
	Redim(Length() + entry.disk_size * chunk_size);
      }
      SetPosition(entry.patch_start * chunk_size);
      WriteBuffer(ptime.patch->start, entry.disk_size * chunk_size);
    }
    /*    FILE *fo = fopen("tmp", "wb+");
    fwrite(entry.patch->start, 1, entry.disk_size * chunk_size, fo);
    fclose(fo);
    exit(0);*/
  }

  if(FlushVbo(ptime))
    vbo_used -= entry.ram_size;

  delete [](ptime.patch->start);
  delete ptime.patch;
  ptime.patch = NULL;

  

  entry.lru_pos = 0xffffffff;
  ram_used -= entry.ram_size;
  ram_flushed += entry.ram_size;
}

bool PatchServer::FlushVbo(PTime &ptime) {
  //TODO  
  //cerr << "Flushing vbo: " << patch << endl;
  if(!ptime.vbo_element) return false;

  glDeleteBuffersARB(1, &ptime.vbo_element);
  glDeleteBuffersARB(1, &ptime.vbo_array);
  ptime.vbo_element = 0;
  ptime.vbo_array = 0;
  return true;
}

void PatchServer::SetRamBufferSize(unsigned int r_buffer) {
  cerr << "Chunk_size: " << chunk_size << endl;
  ram_size = (unsigned int)(r_buffer/chunk_size) + 1;
}
