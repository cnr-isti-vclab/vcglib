#ifndef NXS_PATCH_H
#define NXS_PATCH_H

#include <vcg/space/point3.h>

namespace nxs {

struct Chunk {
  unsigned char p[4096];
};

class Patch {
 public:
  
  Patch(Chunk *s = NULL): start(s) {}

  unsigned short &VertSize() { return *(unsigned short *)start; }

  vcg::Point3f *VertBegin() { 
    return (vcg::Point3f *)(start + 2*sizeof(short)); }

  unsigned short &FaceSize() { return *(((unsigned short *)start) + 1); }

  unsigned short *FaceBegin() {
    return (unsigned short *)(start + 2*sizeof(short) + 
			      VertSize() * sizeof(vcg::Point3f)); }

  unsigned int ChunkSize() {
    return ChunkSize(VertSize(), FaceSize());
  }

  unsigned int ByteSize() {
    return ByteSize(VertSize(), FaceSize());
  }

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
    
    size = 2 * sizeof(unsigned short);    
    return size;
  }
 private:
  Chunk *start;
};

};

#endif
