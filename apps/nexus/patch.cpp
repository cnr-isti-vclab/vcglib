#include "patch.h"

using namespace nxs;


Patch::Patch(Signature signature, Chunk *s, 
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

  cstart = nv * sizeof(float) * 3;
  if(signature & NXS_COLORS)
    nstart = cstart + nv * sizeof(unsigned int);
  else 
    nstart = cstart;

  if(signature & NXS_NORMALS_SHORT)
    tstart = nstart + nv * sizeof(short) * 4;
  else if(signature & NXS_NORMALS_FLOAT)
    tstart = nstart + nv * sizeof(float) * 3;
  else 
    tstart = nstart;

  if(signature & NXS_TEXTURES_SHORT)
    dstart = tstart + nv * sizeof(short) * 2;
  else if(signature & NXS_TEXTURES_FLOAT)
    dstart = tstart + nv * sizeof(float) * 2;
  else 
    dstart = tstart;
}

unsigned int Patch::ChunkSize(Signature signature, 
			      unsigned short nvert, 
			      unsigned short nface) {
  unsigned int size = ByteSize(signature, nvert, nface);
  size = (size/sizeof(Chunk) + 1);
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
  if(size & 0x2) size += 2;
  
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
  if(signature & NXS_DATA16)
    size += nvert * 2 * sizeof(char);
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
