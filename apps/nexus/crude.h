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
Revision 1.5  2004/11/30 22:50:30  ponchio
Level 0.

Revision 1.4  2004/07/05 15:49:38  ponchio
Windows (DevCpp, mingw) port.

Revision 1.3  2004/07/04 15:25:33  ponchio
Backup (work in progress)

Revision 1.2  2004/06/25 16:47:13  ponchio
Various debug

Revision 1.1  2004/06/24 14:32:45  ponchio
Moved from wrap/nexus

Revision 1.1  2004/06/22 15:31:40  ponchio
Created


****************************************************************************/

#ifndef NXS_CRUDE_H
#define NXS_CRUDE_H


#include <string>
#include <vcg/space/point3.h>
#include <vcg/space/box3.h>

#include "vfile.h"

namespace nxs {

class Crude {
public:

  struct Face {
    Face() {}
    Face(unsigned int a, unsigned int b, unsigned int c) {
      v[0] = a; v[1] = b; v[2] = c;
    }
    unsigned int v[3];
    unsigned int &operator[](int k) { return v[k]; }
    unsigned int *ptr() { return v; }
  };

  VFile<vcg::Point3f> vert;
  VFile<Face> face;

  Crude(): fp(NULL), nvert(0), nface(0) {}
  ~Crude();
  
  bool Create(const std::string &file, unsigned int nvert = 0, 
	      unsigned int nface = 0);
  bool Load(const std::string &file, bool rdonly = false);
  void Close();
  void Resize(unsigned int nvert, unsigned int nface);

  unsigned int Vertices();		
  unsigned int Faces();

  vcg::Point3f GetVertex(unsigned int i);
  void SetVertex(unsigned int i, vcg::Point3f &p);
  void SetVertex(unsigned int i, float *p);

  Face GetFace(unsigned int i);
  void SetFace(unsigned int i, Face &f);
  void SetFace(unsigned int i, unsigned int *f);

  vcg::Point3f GetBari(unsigned int i);	

  vcg::Box3f &GetBox();

  /*template <class MESH> void Export(MESH &mesh) {
    MESH::VertexType v;
    v.ClearFlags();
    for(unsigned int i = 0; i < <Vertices(); i++) {
      v.P().Import(GetVertex(i));
      mesh.vert.push_back(v);
    }
    
    MESH::FaceType f;
    f.ClearFlags();
    for(unsigned int i = 0; i < Faces(); i++) {
      for(int k = 0; k < 3; k++) {
	f.
	mesh.vert.push_back(GetVertex(i));
      }
    }
  }*/


 protected:

  FILE *fp;
  unsigned int nvert;
  unsigned int nface;	
  vcg::Box3f box;
};

} //namespace

#endif
