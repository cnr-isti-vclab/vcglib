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
  frame = 0;
  ram_size = rsize;
  ram_used = 0;

  vbo_size = 0;
  vbo_used = 0;

  ram_readed = 0;
  ram_flushed = 0;
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

  ram_readed = 0;
  ram_flushed = 0;
  vbo_size = 0;
  vbo_used = 0;
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
  entry.patch = NULL;
  entry.ram_size = Patch::ChunkSize(signature, nvert, nface, chunk_size);

  //if compressed we do not allocate space now (how much anyway?)
  //  if((signature & NXS_COMPRESSED) != 0) {
    entry.disk_size = 0xffff;
    entry.patch_start = 0xffffffff;
    /*  }  else {
    entry.disk_size = entry.ram_size;
    entry.patch_start = Length()/chunk_size;
    Redim(Length() + entry.disk_size * chunk_size);
    }*/
  entry.lru_pos = 0xffffffff;
  patches.push_back(entry);
  

}

Patch &PatchServer::GetPatch(unsigned int idx, 
			     unsigned short nvert, unsigned short nface,
			     bool flush) {
  assert(idx < patches.size());
  PatchEntry &entry = patches[idx];

  if(entry.patch) {  //already on buffer
    assert(entry.lru_pos < lru.size());
    assert(lru[entry.lru_pos].patch == idx);
    lru[entry.lru_pos].frame = frame++;

  } else {

    assert(entry.lru_pos == 0xffffffff);
    if(flush) Flush();


    char *ram = new char[entry.ram_size * chunk_size];
    entry.patch = new Patch(signature, ram, nvert, nface);
    
    if(entry.patch_start != 0xffffffff) { //was allocated.
      assert(entry.disk_size != 0xffff);

      SetPosition(entry.patch_start * chunk_size);
    
      if((signature & NXS_COMPRESSED) == 0) { //not compressed
	ReadBuffer(ram, entry.disk_size * chunk_size);
      } else {

	unsigned char *disk = new unsigned char[entry.disk_size * chunk_size];
	ReadBuffer(disk, entry.disk_size * chunk_size);

	entry.patch->Decompress(entry.ram_size * chunk_size, 
				disk, entry.disk_size * chunk_size);
	delete []disk;
      } 
    }

    entry.lru_pos = lru.size();
    lru.push_back(PTime(idx, frame++));
    ram_used += entry.ram_size;
    ram_readed += entry.ram_size;
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
  return *(entry.patch); 
}

VboBuffer &PatchServer::GetVbo(unsigned int p) {
  VboBuffer &buffer = vbos[p];
  if(buffer.index) return buffer;

  //TODO  cerr << "Adding vbo: " << p << endl;
  assert(patches[p].patch);
  Patch &patch = *patches[p].patch;

  glGenBuffersARB(1, &buffer.index);
  assert(buffer.index);
  glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER_ARB, buffer.index);

  unsigned int size = patch.nf * sizeof(unsigned short);
  if((signature & NXS_FACES) != 0) size *= 3;

  glBufferDataARB(GL_ELEMENT_ARRAY_BUFFER_ARB, size, patch.FaceBegin(),
  GL_STATIC_DRAW_ARB);

  //TODO fix this when we allow data :p
  size = sizeof(float) * patch.dstart;
  
  glGenBuffersARB(1, &buffer.vertex);
  glBindBufferARB(GL_ARRAY_BUFFER_ARB, buffer.vertex);

  glBufferDataARB(GL_ARRAY_BUFFER_ARB, size, patch.VertBegin(), 
		  GL_STATIC_DRAW_ARB);

  vbo_used += patches[p].ram_size;
  return buffer;
}


void PatchServer::Flush() {

  if(ram_used < ram_size * 1.1) return;

  make_heap(lru.begin(), lru.end());
  for(unsigned int i = 0; i < lru.size(); i++)
    patches[lru[i].patch].lru_pos = i;

  while(ram_used > ram_size) {
    pop_heap(lru.begin(), lru.end());
    PTime &ptime = lru.back();
    Flush(ptime.patch);
    lru.pop_back();
  }
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


    if((signature & NXS_COMPRESSED)) {
      unsigned int compressed_size;
      char *compressed = entry.patch->Compress(entry.ram_size * chunk_size,
					       compressed_size);
      if(entry.disk_size == 0xffff) {//allocate space 
	assert(entry.patch_start == 0xffffffff);
	entry.disk_size = (unsigned int)((compressed_size-1)/chunk_size) + 1;
	entry.patch_start = Length()/chunk_size;
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
	entry.patch_start = Length()/chunk_size;
	Redim(Length() + entry.disk_size * chunk_size);
      }
      SetPosition(entry.patch_start * chunk_size);
      WriteBuffer(entry.patch->start, entry.disk_size * chunk_size);
    }
    /*    FILE *fo = fopen("tmp", "wb+");
    fwrite(entry.patch->start, 1, entry.disk_size * chunk_size, fo);
    fclose(fo);
    exit(0);*/
  }

  delete [](entry.patch->start);
  delete entry.patch;
  entry.patch = NULL;
  entry.lru_pos = 0xffffffff;
  if(FlushVbo(patch))
    vbo_used -= entry.ram_size;
  ram_used -= entry.ram_size;
  ram_flushed += entry.ram_size;
}

bool PatchServer::FlushVbo(unsigned int patch) {
  //TODO  
  //cerr << "Flushing vbo: " << patch << endl;
  VboBuffer &buffer = vbos[patch];
  if(!buffer.index) return false;
  glDeleteBuffersARB(1, &buffer.index);
  glDeleteBuffersARB(1, &buffer.vertex);
  return true;
}

void PatchServer::SetRamBufferSize(unsigned int r_buffer) {
  ram_size = (unsigned int)(r_buffer/chunk_size) + 1;
}
