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
Revision 1.11  2005/02/08 12:43:03  ponchio
Added copyright


****************************************************************************/

#ifndef NXS_PATCH_H
#define NXS_PATCH_H

#include <vcg/space/point3.h>
#include <vcg/space/sphere3.h>

namespace nxs {

 struct Signature {

   enum Face { TRIANGLES = 1, STRIPS = 2, TETRAS = 3, SLICE = 4 };

   enum Vert { POINT2F = 1, POINT2D = 2, 
	       POINT3F = 2, POINT3D = 3,
	       POINT4F = 4, POINT4D = 5 };

   enum Compr { LZO = 1 };

   unsigned char face;   
   unsigned char vert;   
   unsigned char compr;  
   unsigned char future; //who knows...

   unsigned char fcolor;
   unsigned char fnorm;
   unsigned char ftext;
   unsigned char fdata;

   unsigned char vcolor;
   unsigned char vnorm;
   unsigned char vtext;
   unsigned char vdata;

   Signature(): face(1), vert(2), compr(0), future(0),
	fcolor(0),  fnorm(0), ftext(0), fdata(0),
	vcolor(0),  vnorm(0), vtext(0), vdata(0) {}
 };


 struct Encoding {
   
   unsigned char bytes; //size per element
   unsigned char comps; //number of components
   void (*pack)(char *start, unsigned int nelem);
   void (*unpack)(char *start, unsigned int nelem);
   
   unsigned int size(unsigned short n) { 
     unsigned int s = (int)n * (int)bytes * (int)comps; 
     //padding a 64 bytes
     if((s & 0x0000003f) != 0) {
       s>>=6; s++; s<<=6;
     }
     return s;
   }
 };
 
 struct Encodings {
   enum Name { EMPTY = 0, 
	       BYTE1 = 1,  SHORT1 = 2,  FLOAT1 = 3,  DOUBLE1 = 4, 
	       BYTE2 = 5,  SHORT2 = 6,  FLOAT2 = 7,  DOUBLE2 = 8, 
	       BYTE3 = 9,  SHORT3 = 10, FLOAT3 = 11, DOUBLE3 = 12, 
	       BYTE4 = 13, SHORT4 = 14, FLOAT4 = 15, DOUBLE4 = 16 };
   Encodings();
   Encoding &operator[](int n) { return e[n]; }
 protected:
   Encoding e[17];
   
   
 };

class Patch {
 public:
 
  static Encodings encodings;

  Patch(Signature &signature, char *s, 
	unsigned short nv, unsigned short nf);

  void Init(Signature &signature, unsigned short nv, unsigned short nf);

  vcg::Point3f *Vert3fBegin() { return (vcg::Point3f *)vstart; }
  vcg::Point3f &Vert3f(int n) { return Vert3fBegin()[n]; }
  unsigned short *FaceBegin() { return (unsigned short *)fstart; }
  unsigned short *Face(int n) { return FaceBegin() + 3 * n; }

  //vcg::Point3f &Vert(unsigned short v)   { return VertBegin()[v]; }
  //  unsigned short *Face(unsigned short f) { return FaceBegin() + f * 3; }

  char *VColorBegin() { return vstart + 64*vstartc; }
  char *VNormBegin()  { return vstart + 64*vstartn; }
  char *VTextBegin()  { return vstart + 64*vstartt; }
  char *VDataBegin()  { return vstart + 64*vstartd; }

  char *FColorBegin() { return fstart + 64*fstartc; }
  char *FNormBegin()  { return fstart + 64*fstartn; }
  char *FTextBegin()  { return fstart + 64*fstartt; }
  char *FDataBegin()  { return fstart + 64*fstartd; }

  static unsigned int ChunkSize(Signature &signature, 
				unsigned short nvert, 
				unsigned short nface,
				unsigned int chunk_size);

  static unsigned int ByteSize(Signature &signature, 
			       unsigned short nvert, 
			       unsigned short nface);
  
  char *Compress(unsigned int ram_size, unsigned int &size);
  void Decompress(unsigned int ram_size, void *src, unsigned int src_sz);
		  

  char *fstart;
  char *vstart;

  unsigned short nf;
  unsigned short nv;

  //these offset are from fstart in 64 bytes
  unsigned short fstartc;
  unsigned short fstartn;
  unsigned short fstartt;
  unsigned short fstartd;

  //these offset are from vstart in 64 bytes
  unsigned short vstartc;
  unsigned short vstartn;
  unsigned short vstartt;
  unsigned short vstartd;
};

} //namespace

#endif
