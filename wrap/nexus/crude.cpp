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

#include "crude.h"

using namespace std;
using namespace vcg;
using namespace nxs;


Crude::~Crude() {
  if(fp)
    Close();
}

bool Crude::Create(const std::string &file, unsigned int nv, unsigned int nf) {
  if(!vert.Create(file + ".vrt")) return false;
  if(!face.Create(file + ".frt")) return false;

  fp = fopen(file.c_str(), "wb+");
  if(!fp) return false;
  Resize(nv, nf);
  return true;
}
bool Crude::Load(const std::string &file) {
  if(!vert.Load(file + ".vrt")) return false;
  if(!face.Load(file + ".frt")) return false;

  fp = fopen(file.c_str(), "rb+");
  if(!fp) return false;
  fread(&nvert, sizeof(unsigned int), 1, fp);
  fread(&nface, sizeof(unsigned int), 1, fp);
  fread(&box, sizeof(Box3f), 1, fp);
  return true;                
}
void Crude::Close() {
  vert.Close();
  face.Close();
  rewind(fp);
  fwrite(&nvert, sizeof(unsigned int), 1, fp);
  fwrite(&nface, sizeof(unsigned int), 1, fp);
  fwrite(&box, sizeof(Box3f), 1, fp);
  fclose(fp);
  fp = NULL;
}

void Crude::Resize(unsigned int nv, unsigned int nf) {
  nvert = nv;
  nface = nf;
  vert.Resize(nv);
  face.Resize(nf);
}

unsigned int Crude::Vertices() {
  return nvert;
}
unsigned int Crude::Faces() {
  return nface;
}

void Crude::SetVertex(unsigned int i, float *f) {
  Point3f &p = vert[i];
  p[0] = f[0];
  p[1] = f[1];
  p[2] = f[2];
}

Point3f &Crude::GetVertex(unsigned int i) {
  return vert[i];
}
Crude::Face &Crude::GetFace(unsigned int i) {
  return face[i];
}
 void Crude::SetFace(unsigned int i, unsigned int *f) { 
   Face &ff = face[i];
   ff[0] = f[0];
   ff[1] = f[1];
   ff[2] = f[2];
 }
  
vcg::Point3f Crude::GetBari(unsigned int i) {
  Point3f bari(0, 0, 0);
  Face &f = face[i];
  for(int k = 0; k < 3; k++)
    bari += vert[f[k]];
  bari /= 3;
  return bari;
}

vcg::Box3f &Crude::GetBox() {
  return box;
}

  
