#include "patch.h"
#ifdef WIN32
#include "minilzo.108/minilzo.h"
#else
#include <lzo1x.h>
#endif
#include <iostream>
using namespace std;
using namespace nxs;

#ifdef WIN32
static double wrkmem[LZO1X_1_MEM_COMPRESS/sizeof(double) +1];
#else
static double wrkmem[LZO1X_999_MEM_COMPRESS/sizeof(double) +1];
#endif


void pad(unsigned int &size) {
  while(size&0x3) size++;
}

void shuffle(float *buffer, unsigned int size, unsigned int stride) {
  float *tmp = new float[size];
  unsigned int count = 0;

  unsigned int nelem = size/stride;
  for(unsigned int s = 0; s < stride; s++) {
    float *ptr = buffer + s;
    for(unsigned int i = 0; i < nelem; i++) {
      tmp[count++] = *ptr;
      ptr += stride;
    }
  }
  memcpy(buffer, tmp, size * sizeof(float));
  delete []tmp;
}

void unshuffle(float *buffer, unsigned int size, unsigned int stride) {
  float *tmp = new float[size];

  unsigned int count = 0;
  unsigned int nelem = size/stride;
  for(unsigned int s = 0; s < stride; s++) {
    float *ptr = tmp + s;
    for(unsigned int i = 0; i < nelem; i++) {
      *ptr = buffer[count++];
      ptr += stride;
    }
  }
  memcpy(buffer, tmp, size * sizeof(float));
  delete []tmp;
}

void subtract(float *buffer, unsigned int size) {
  float p = buffer[0];
  float q;
  for(unsigned int i = 1; i < size; i++) {
    q = buffer[i];
    buffer[i] -= p;
    p = q;
  }
}

void unsubtract(float *buffer, unsigned int size) {
  for(unsigned int i = 1; i < size; i++) 
    buffer[i] += buffer[i-1];
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

  //lets use differences
  //  shuffle((float *)VertBegin(), nv * 3, 3);
  //  subtract((float *)VertBegin(), nv * 3);

  //TODO use OVERLAP and test speed
  //TODO fill chunk padding with zeroes?
  //TODO compress only used memory!
  size = ram_size + ram_size/64 + 23;
  char *buffer = new char[size];
#ifdef WIN32
    lzo1x_1_compress(((unsigned char *)start), ram_size,
		     (unsigned char *)buffer + sizeof(int), &size,
		   (char *)wrkmem);
#else
    lzo1x_999_compress(((unsigned char *)start), ram_size,
		       (unsigned char *)buffer + sizeof(int), &size,
		   (char *)wrkmem);

  lzo1x_optimize((unsigned char *)buffer + sizeof(int), size,
		 ((unsigned char *)start), &ram_size,		 
		 NULL);
#endif

  *(int *)buffer = size;
  size += sizeof(int);

  //  memcpy(buffer, start, ram_size);
  //  size = ram_size;

  return buffer;

}

void Patch::Decompress(unsigned int ram_size, void *src, unsigned int src_sz) {
  
  unsigned int size = *(int *)src;
  assert(size < src_sz + sizeof(int));
  unsigned int dst_size = ram_size;
    
  int ret = lzo1x_decompress_safe(((unsigned char *)src) + sizeof(int), size,
				  (unsigned char *)start, &dst_size, 0);
  if(ret != 0) {
    cerr << "Ret from decompress: " << ret << endl;
    exit(-1);
  }
  assert(dst_size == ram_size);
  //TODO add 3 to start... so we can use asm_fast decompressor  

  //  unsubtract((float *)VertBegin(), nv * 3);
  //  unshuffle((float *)VertBegin(), nv * 3, 3);
}

