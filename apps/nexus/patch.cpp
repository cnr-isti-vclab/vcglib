#include "patch.h"
#include <lzo1x.h>
#include <iostream>
using namespace std;
using namespace nxs;

static double wrkmem[LZO1X_999_MEM_COMPRESS/sizeof(double) +1];

void pad(unsigned int &size) {
  while(size&0x3) size++;
}

Patch::Patch(Signature signature, char *s, 
	     unsigned short nvert, unsigned short nface):
             start(s) {
  Init(signature, nvert, nface);
}

void Patch::Init(Signature signature, 
		 unsigned short nvert, unsigned short nface) {
  nv = nvert;
  nf = nface;
  
  if(signature & NXS_FACES)
    vstart = (float *)(((char *)start) + nf * sizeof(unsigned short) * 3);
  else if(signature & NXS_STRIP)
    vstart = (float *)(((char *)start) + nf * sizeof(unsigned short));
  else
    vstart = (float *)start;

  //align memory
  if(((int)vstart) & 0x2) vstart = (float *)(((char *)vstart) + 2);

  cstart = nv * 3;
  if(signature & NXS_COLORS)
    nstart = cstart + nv;
  else 
    nstart = cstart;

  if(signature & NXS_NORMALS_SHORT)
    tstart = nstart + nv * 2;
  else if(signature & NXS_NORMALS_FLOAT)
    tstart = nstart + nv * 3;
  else 
    tstart = nstart;

  if(signature & NXS_TEXTURES_SHORT)
    dstart = tstart + nv;
  else if(signature & NXS_TEXTURES_FLOAT)
    dstart = tstart + nv;
  else 
    dstart = tstart;
}

unsigned int Patch::ChunkSize(Signature signature, 
			      unsigned short nvert, 
			      unsigned short nface,
			      unsigned int chunk_size) {
  unsigned int size = ByteSize(signature, nvert, nface);
  size = (size/chunk_size + 1);
  return size;
}

unsigned int Patch::ByteSize(Signature signature, 
			     unsigned short nvert, 
			     unsigned short nface) {
  unsigned int size = 0;
  if(signature & NXS_STRIP)
    size += nface * sizeof(unsigned short);
  else if(signature & NXS_FACES)
    size += nface * 3 * sizeof(unsigned short);
  
  //memory alignment
  pad(size);
  
  size += nvert * sizeof(vcg::Point3f);
  
  if(signature & NXS_COLORS)
    size += nvert * sizeof(unsigned int);
  
  if(signature & NXS_NORMALS_SHORT)
    size += nvert * 4 * sizeof(short);
  
  if(signature & NXS_NORMALS_FLOAT)
    size += nvert * 3 * sizeof(float);
  
  if(signature & NXS_TEXTURES_SHORT)
    size += nvert * 2 * sizeof(short);
  
  if(signature & NXS_TEXTURES_FLOAT)
    size += nvert * 2 * sizeof(float);
  
  if(signature & NXS_DATA8)
    size += nvert * sizeof(char);
  pad(size);
  if(signature & NXS_DATA16)
    size += nvert * 2 * sizeof(char);
  pad(size);
  if(signature & NXS_DATA32)
    size += nvert * 4 * sizeof(char);
  if(signature & NXS_DATA64)
    size += nvert * 8 * sizeof(char);
  
  
  //this condition should really rarely happen but helps save space
  //during construction
  if(size < nface * 3 * sizeof(unsigned int))
    size = nface * 3 * sizeof(unsigned int);
  
  return size;
}


char *Patch::Compress(unsigned int ram_size, unsigned int &size) {
  //TODO use OVERLAP and test speed
  //TODO fill chunk padding with zeroes?
  //TODO compress only used memory!
  size = ram_size + ram_size/64 + 23;
  char *buffer = new char[size];
  lzo1x_1_compress(((unsigned char *)start), ram_size,
		     (unsigned char *)buffer + sizeof(int), &size,
		   (char *)wrkmem);
  
  *(int *)buffer = size;
  size += sizeof(int);
  

  //  memcpy(buffer, start, ram_size);
  //  size = ram_size;
  //TODO optimize!
  //      lzo1x_optimize((unsigned char *)entry.patch->start,
  //			 entry.ram_size * chunk_size,
  //			 compressed, &compressed_size, NULL);


  return buffer;
}

void Patch::Decompress(unsigned int ram_size, void *src, unsigned int src_sz) {
		     
  unsigned int size = *(int *)src;
  assert(size < src_sz + sizeof(int));
  unsigned int dst_size = ram_size;
  
  //  memcpy(start, src, ram_size);
  int ret = lzo1x_decompress_safe(((unsigned char *)src) + sizeof(int), size,
				  (unsigned char *)start, &dst_size, 0);
  if(ret != 0) {
    cerr << "Ret from decompress: " << ret << endl;
    exit(-1);
  }
  assert(dst_size == ram_size);
  //TODO add 3 to start... so we can use asm_fast decompressor  
}

