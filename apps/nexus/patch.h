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

****************************************************************************/

#ifndef NXS_PATCH_H
#define NXS_PATCH_H

#include <vcg/space/point3.h>
#include <vcg/space/sphere3.h>

namespace nxs {

enum Signature { NXS_FACES          = 0x00000001,
		 NXS_STRIP          = 0x00000002, 
		 NXS_COLORS         = 0x00000010, 
		 NXS_NORMALS_SHORT  = 0x00000100, 
		 NXS_NORMALS_FLOAT  = 0x00000200, 
		 NXS_TEXTURES_SHORT = 0x00001000,  
		 NXS_TEXTURES_FLOAT = 0x00002000, 
		 NXS_DATA8          = 0x00010000,  
		 NXS_DATA16         = 0x00020000, 
		 NXS_DATA32         = 0x00040000, 
		 NXS_DATA64         = 0x00080000,
                 NXS_COMPRESSED      = 0x10000000};



class Patch {
 public:
 
  Patch(Signature signature, char *s, 
	unsigned short nv, unsigned short nf);

  void Init(Signature signature, unsigned short nv, unsigned short nf);

  inline vcg::Point3f *VertBegin();
  inline unsigned short *FaceBegin();

  inline vcg::Point3f &Vert(unsigned short v);
  inline unsigned short *Face(unsigned short f);

  inline unsigned int *ColorBegin();
  inline short *Norm16Begin();
  inline short *Norm16(unsigned short v);
  inline vcg::Point3f *Norm32Begin();
  inline vcg::Point3f &Norm32(unsigned short v);
 
  static unsigned int ChunkSize(Signature signature, 
				unsigned short nvert, 
				unsigned short nface,
				unsigned int chunk_size);

  static unsigned int ByteSize(Signature signature, 
			       unsigned short nvert, 
			       unsigned short nface);
  
  char *Compress(unsigned int ram_size, unsigned int &size);
  void Decompress(unsigned int ram_size, void *src, unsigned int src_sz);
		  

  unsigned short nv;
  unsigned short nf;

  char *start;
  float *vstart;
  //these offset are from vstart!
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

inline unsigned int *Patch::ColorBegin() { 
  return (unsigned int *)(vstart + cstart); 
}

inline short *Patch::Norm16Begin() { 
  return (short *)(vstart + nstart); 
}

inline short *Patch::Norm16(unsigned short v) { 
  return Norm16Begin() + 4 * v;
}

 inline vcg::Point3f *Patch::Norm32Begin() {
   return (vcg::Point3f *)(vstart + nstart); 
 }
 inline vcg::Point3f &Patch::Norm32(unsigned short v) {
   return Norm32Begin()[v];
 }

} //namespace

#endif
