/****************************************************************************
* VCGLib                                                            o o     *
* Visual and Computer Graphics Library                            o     o   *
*                                                                _   O  _   *
* Copyright(C) 2004                                                \/)\/    *
* Visual Computing Lab                                            /\/|      *
* ISTI - Italian National Research Council                           |      *
*                                                                    \      *
* All rights reserved.                                                      *
*                                                                           *
* This program is free software; you can redistribute it and/or modify      *   
* it under the terms of the GNU General Public License as published by      *
* the Free Software Foundation; either version 2 of the License, or         *
* (at your option) any later version.                                       *
*                                                                           *
* This program is distributed in the hope that it will be useful,           *
* but WITHOUT ANY WARRANTY; without even the implied warranty of            *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
* GNU General Public License (http://www.gnu.org/licenses/gpl.txt)          *
* for more details.                                                         *
*                                                                           *
****************************************************************************/
/****************************************************************************
  History

$Log: not supported by cvs2svn $
Revision 1.13  2005/02/22 10:38:15  ponchio
Debug, cleaning and optimization.

Revision 1.12  2005/02/21 19:05:58  ponchio
i already fixed this bug. I hate you cvs.

Revision 1.11  2005/02/19 12:06:55  ponchio
Debug...

Revision 1.10  2005/02/19 10:45:05  ponchio
Patch generalized and small fixes.

Revision 1.9  2005/02/08 12:43:03  ponchio
Added copyright


****************************************************************************/

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


Encodings Patch::encodings;

Encodings::Encodings() {
  for(unsigned int i = 0; i < 17; i++) {
    e[i].bytes = 0;
    e[i].comps = 0;
    e[i].pack = NULL;
    e[i].unpack = NULL;
  }
  e[1].bytes = 1;
  e[2].bytes = 2;
  e[3].bytes = 4;
  e[4].bytes = 8;
  e[1].comps = e[2].comps = e[3].comps = e[4].comps = 1;
  e[5].bytes = 1;
  e[6].bytes = 2;
  e[7].bytes = 4;
  e[8].bytes = 8;
  e[5].comps = e[6].comps = e[7].comps = e[8].comps = 2;
  e[9].bytes = 1;
  e[10].bytes = 2;
  e[11].bytes = 4;
  e[12].bytes = 8;
  e[9].comps = e[10].comps = e[11].comps = e[12].comps = 3;
  e[13].bytes = 1;
  e[14].bytes = 2;
  e[15].bytes = 4;
  e[16].bytes = 8;
  e[13].comps = e[14].comps = e[15].comps = e[16].comps = 4;
}


void pad8(unsigned int &s) {
  if((s & 0x00000007) != 0) {
    s>>=3; s++; s<<=3;
  }
}
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


Patch::Patch(Signature &signature, char *s, 
	     unsigned short nvert, unsigned short nface):
  fstart(s) {
  Init(signature, nvert, nface);
}

void Patch::Init(Signature &signature, 
		 unsigned short nvert, unsigned short nface) {
  nv = nvert;
  nf = nface;
  
  unsigned int offset = 0;

  if(signature.face == Signature::TRIANGLES)
    offset += nf * 3 * sizeof(unsigned short);
  else if (signature.face == Signature::STRIPS)
    offset += nf * sizeof(unsigned short);
  else if (signature.face == Signature::TETRAS)
    offset += nf * 4 * sizeof(unsigned short);
  else if (signature.face == Signature::SLICE) {
    assert(0);
    //non lo so...
  }
  pad8(offset);

  fstartc = fstart + offset;
  offset += encodings[signature.fcolor].size(nf);
  fstartn = fstart + offset;
  offset += encodings[signature.fnorm].size(nf);
  fstartt = fstart + offset;
  offset += encodings[signature.ftext].size(nf);
  fstartd = fstart + offset;
  offset += encodings[signature.fdata].size(nf);

  vstart = fstart + offset;
  if(signature.vert == Signature::POINT3F)
    offset += nv * sizeof(float) * 3;
  else if(signature.vert == Signature::POINT4F)
    offset += nv * sizeof(float) * 4;
  else 
    assert(0);
  pad8(offset);

  vstartc = fstart + offset;
  offset += encodings[signature.vcolor].size(nv);
  vstartn = fstart + offset;
  offset += encodings[signature.vnorm].size(nv);
  vstartt = fstart + offset;
  offset += encodings[signature.vtext].size(nv);
  vstartd = fstart + offset;
  offset += encodings[signature.vdata].size(nv);
  

  /*  if(signature & NXS_FACES)
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
  dstart = tstart;*/
}

unsigned int Patch::ChunkSize(Signature &signature, 
			      unsigned short nvert, 
			      unsigned short nface,
			      unsigned int chunk_size) {
  unsigned int size = ByteSize(signature, nvert, nface);
  size = (size/chunk_size + 1);
  return size;
}

unsigned int Patch::ByteSize(Signature &signature, 
			     unsigned short nvert, 
			     unsigned short nface) {

  unsigned int size = 0;
  if(signature.face == Signature::TRIANGLES)
    size += nface * 3 * sizeof(unsigned short);
  else if (signature.face == Signature::STRIPS)
    size += nface * sizeof(unsigned short);
  else if (signature.face == Signature::TETRAS)
    size += nface * 4 * sizeof(unsigned short);
  else if (signature.face == Signature::SLICE) {
    assert(0);
    //non lo so...
  }
  pad8(size);

  size += encodings[signature.fcolor].size(nface);
  size += encodings[signature.fnorm].size(nface);
  size += encodings[signature.ftext].size(nface);
  size += encodings[signature.fdata].size(nface);

  if(signature.vert == Signature::POINT3F)
    size += nvert * sizeof(float) * 3;
  else if(signature.vert == Signature::POINT4F)
    size += nvert * sizeof(float) * 4;
  else 
    assert(0);
  pad8(size);

  size += encodings[signature.vcolor].size(nvert);
  size += encodings[signature.vnorm].size(nvert);
  size += encodings[signature.vtext].size(nvert);
  size += encodings[signature.vdata].size(nvert);

  //this condition should really rarely happen but helps save space
  //during construction
  if(size < nface * 3 * sizeof(unsigned int))
    size = nface * 3 * sizeof(unsigned int);
  
    return size;


  /*  unsigned int size = 0;
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
  
    return size;*/
}


char *Patch::Compress(unsigned int ram_size, unsigned int &size) {

  //lets use differences
  //  shuffle((float *)VertBegin(), nv * 3, 3);
  //  subtract((float *)VertBegin(), nv * 3);

  //TODO use OVERLAP and test speed
  //TODO fill chunk padding with zeroes?
  size = ram_size + ram_size/64 + 23;
  char *buffer = new char[size];
#ifdef WIN32
    lzo1x_1_compress(((unsigned char *)fstart), ram_size,
		     (unsigned char *)buffer + sizeof(int), &size,
		   (char *)wrkmem);
#else
    lzo1x_999_compress(((unsigned char *)fstart), ram_size,
		       (unsigned char *)buffer + sizeof(int), &size,
		   (char *)wrkmem);

  lzo1x_optimize((unsigned char *)buffer + sizeof(int), size,
		 ((unsigned char *)fstart), &ram_size,		 
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
				  (unsigned char *)fstart, &dst_size, 0);
  if(ret != 0) {
    cerr << "Ret from decompress: " << ret << endl;
    exit(-1);
  }
  assert(dst_size == ram_size);
  //TODO add 3 to start... so we can use asm_fast decompressor  

  //  unsubtract((float *)VertBegin(), nv * 3);
  //  unshuffle((float *)VertBegin(), nv * 3, 3);
}

