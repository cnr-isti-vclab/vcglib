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
 
  enum Signature { DEFAULT = 0,
              HAS_STRIP          = 0x002, //if true faces saved as strip
	      HAS_COLORS         = 0x004, 
	      HAS_NORMALS_SHORT  = 0x008, 
	      HAS_NORMALS_FLOAT  = 0x008, 
	      HAS_TEXTURES_SHORT = 0x010,  
	      HAS_TEXTURES_FLOAT = 0x020, 
	      HAS_DATA8          = 0x040,  
	      HAS_DATA16         = 0x080, 
	      HAS_DATA32         = 0x100, 
	      HAS_DATA64         = 0x200 };
   
  Patch(Signature signature, Chunk *s = NULL, 
	unsigned short nv = 0, unsigned short nf = 0): 
    start(s) {
    Resize(signature, nv, nf);
  }

  void Resize(Signature signature, unsigned short nv, unsigned short nf) {
    nvert = nv;
    nface = nf;
    fstart = (unsigned short *)(((char *)start) + 
				VertSize() * sizeof(vcg::Point3f));

    unsigned int size = nf * sizeof(unsigned short);
    if(!(signature & HAS_STRIP)) size *= 3;    
    cstart = (void *)(((char *)fstart) + size);

    if(signature & HAS_COLORS) size = nv * sizeof(unsigned int);
    else size = 0;

    nstart = (void *)(((char *)cstart) + size);
    
    size = 0;
    if(signature & HAS_NORMALS_SHORT) size = nv * 4*sizeof(short);
    if(signature & HAS_NORMALS_FLOAT) size = nv * 4*sizeof(float);

    tstart = (void *)(((char *)cstart) + size);

    size = 0;
    if(signature & HAS_TEXTURES_SHORT) size = nv * 2*sizeof(short);
    if(signature & HAS_TEXTURES_FLOAT) size = nv * 2*sizeof(float);
    
    dstart = (void *)(((char *)tstart) + size);
  }

  unsigned short VertSize()    { return nvert; }

  vcg::Point3f *VertBegin()     { return (vcg::Point3f *)(start); }

  unsigned short FaceSize()    { return nface; }

  unsigned short *FaceBegin()   { return fstart; }

  vcg::Point3f &Vert(unsigned int v)   { return VertBegin()[v]; }
  
  unsigned short *Face(unsigned int f) { return FaceBegin() + f * 3; }
 
  //  unsigned int ChunkSize() { return ChunkSize(VertSize(), FaceSize()); }

  //  unsigned int ByteSize() { return ByteSize(VertSize(), FaceSize()); }

  static unsigned int ChunkSize(Signature signature, 
				unsigned short nvert, unsigned short nface) {
    unsigned int size = ByteSize(signature, nvert, nface);
    size = (size/sizeof(Chunk) + 1);
    return size;
  }

  static unsigned int ByteSize(Signature signature, 
			       unsigned short nvert, unsigned short nface) {
    unsigned int size = nvert * sizeof(vcg::Point3f);
    if(signature & HAS_STRIP)
      size += nface * sizeof(unsigned short);
    else
      size += nface * 3 * sizeof(unsigned short);
    
    if(signature & HAS_COLORS)
      size += nvert * sizeof(unsigned int);

    if(signature & HAS_NORMALS_SHORT)
      size += nvert * 4 * sizeof(short);

    if(signature & HAS_NORMALS_FLOAT)
      size += nvert * 3 * sizeof(float);

    if(signature & HAS_TEXTURES_SHORT)
      size += nvert * 2 * sizeof(short);

    if(signature & HAS_TEXTURES_FLOAT)
      size += nvert * 2 * sizeof(float);

    if(signature & HAS_DATA8)
      size += nvert * sizeof(char);
    if(signature & HAS_DATA16)
      size += nvert * 2 * sizeof(char);
    if(signature & HAS_DATA32)
      size += nvert * 4 * sizeof(char);
    if(signature & HAS_DATA64)
      size += nvert * 8 * sizeof(char);


    //this condition should really rarely happen but helps save space
    //during construction
    if(size < nface * 3 * sizeof(unsigned int))
      size = nface * 3 * sizeof(unsigned int);
    
    return size;
  }
  // private:
  unsigned short nvert;
  unsigned short nface;

  Chunk *start;
  unsigned short *fstart;
  void *cstart;
  void *nstart;
  void *tstart;
  void *dstart;
};

};

#endif
