#ifndef NXS_PATCH_H
#define NXS_PATCH_H

#include <vcg/space/point3.h>
#include <iostream>
namespace nxs {

struct Chunk {
  unsigned char p[4096];
};

class Patch {
 public:
  
  Patch(Chunk *s = NULL, unsigned short nv = 0, unsigned short nf = 0): 
    start(s) {
    Resize(nv, nf);
  }

  void Resize(unsigned short nv, unsigned short nf) {
    nvert = nv;
    nface = nf;
    fstart = (unsigned short *)(((char *)start) + 
				VertSize() * sizeof(vcg::Point3f));
  }
  unsigned short VertSize()    { return nvert; }

  vcg::Point3f *VertBegin()     { return (vcg::Point3f *)(start); }

  unsigned short FaceSize()    { return nface; }

  unsigned short *FaceBegin()   { return fstart; }

  vcg::Point3f &Vert(unsigned int v)   { return VertBegin()[v]; }
  
  unsigned short *Face(unsigned int f) { return FaceBegin() + f * 3; }
 
  unsigned int ChunkSize() { return ChunkSize(VertSize(), FaceSize()); }

  unsigned int ByteSize() { return ByteSize(VertSize(), FaceSize()); }

  static unsigned int ChunkSize(unsigned short nvert, unsigned short nface) {
    unsigned int size = ByteSize(nvert, nface);
    size = (size/sizeof(Chunk) + 1);
    return size;
  }

  static unsigned int ByteSize(unsigned short nvert, unsigned short nface) {
    unsigned int size = nvert * sizeof(vcg::Point3f);
    size += nface * 3 * sizeof(unsigned short);

    //this condition should really rarely happen but helps save space
    //during construction
    if(size < nface * 3 * sizeof(unsigned int))
      size = nface * 3 * sizeof(unsigned int);
    
    return size;
  }
  // private:
  Chunk *start;
  unsigned short *fstart;
  unsigned short nvert;
  unsigned short nface;
};

};

#endif
