#ifndef NXS_PATCH_H
#define NXS_PATCH_H

#include <vcg/space/point3.h>
#include <iostream>
namespace nxs {

enum Signature { 
		 HAS_FACES          = 0x00000001,
		 HAS_STRIP          = 0x00000002, 
		 HAS_COLORS         = 0x00000010, 
		 HAS_NORMALS_SHORT  = 0x00000100, 
		 HAS_NORMALS_FLOAT  = 0x00000200, 
		 HAS_TEXTURES_SHORT = 0x00001000,  
		 HAS_TEXTURES_FLOAT = 0x00002000, 
		 HAS_DATA8          = 0x00010000,  
		 HAS_DATA16         = 0x00020000, 
		 HAS_DATA32         = 0x00040000, 
		 HAS_DATA64         = 0x00080000 };

struct Chunk {
  unsigned char p[4096];
};


class Patch {
 public:
 
  Patch(Signature signature, Chunk *s, 
	unsigned short nv, unsigned short nf);

  void Init(Signature signature, unsigned short nv, unsigned short nf);

  inline vcg::Point3f *VertBegin();
  inline unsigned short *FaceBegin();

  inline vcg::Point3f &Vert(unsigned short v);
  inline unsigned short *Face(unsigned short f);
 
  static unsigned int ChunkSize(Signature signature, 
				unsigned short nvert, 
				unsigned short nface);

  static unsigned int ByteSize(Signature signature, 
			       unsigned short nvert, 
			       unsigned short nface);


  Chunk *start;

  unsigned short nv;
  unsigned short nf;

  float *vstart;
  unsigned short cstart;
  unsigned short nstart;
  unsigned short tstart;
  unsigned short dstart;
};

inline vcg::Point3f *Patch::VertBegin() { 
  return (vcg::Point3f *)vstart; 
}
inline unsigned short *Patch::FaceBegin() {
  return (unsigned short *)start; 
}

inline vcg::Point3f &Patch::Vert(unsigned short v) {
  return VertBegin()[v]; 
}
inline unsigned short *Patch::Face(unsigned short f) {
  return FaceBegin() + f * 3; 
}

} //namespace

#endif
