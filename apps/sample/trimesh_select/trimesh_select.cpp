/****************************************************************************
* VCGLib                                                            o o     *
* Visual and Computer Graphics Library                            o     o   *
*                                                                _   O  _   *
* Copyright(C) 2004-2016                                           \/)\/    *
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
/*! \file trimesh_select.cpp
\ingroup code_sample

\brief the minimal example of using the selection functionality
*/

#include<vcg/complex/complex.h>
#include<vcg/complex/algorithms/create/platonic.h>

using namespace vcg;
using namespace std;


class MyFace;
class MyVertex;

struct MyUsedTypes : public UsedTypes<	Use<MyVertex>::AsVertexType, Use<MyFace>::AsFaceType>{};
class MyVertex  : public Vertex< MyUsedTypes, vertex::VFAdj, vertex::Coord3f, vertex::Normal3f, vertex::BitFlags  >{};
class MyFace    : public Face  < MyUsedTypes, face::VFAdj, face::Normal3f, face::VertexRef, face::BitFlags > {};
class MyMesh    : public vcg::tri::TriMesh<vector<MyVertex>, vector<MyFace> > {};

int main(int ,char ** )
{
  MyMesh m; 
  // build a 10x10 grid of 100 vertices and 200 faces. 
  tri::Grid(m,10,10,1.0f,1.0f);
  tri::UpdateFlags<MyMesh>::VertexBorderFromNone(m);
  int cornerNum = tri::UpdateSelection<MyMesh>::VertexCornerBorder(m,math::ToRad(100.0f));
  printf("Mesh has %i vertexes  %i faces and %i corners\n",m.VN(),m.FN(),cornerNum);
  
  int borderNum = tri::UpdateSelection<MyMesh>::VertexFromBorderFlag(m);
  int internalNum = tri::UpdateSelection<MyMesh>::VertexInvert(m);
  int faceNotOnBoundaryNum = tri::UpdateSelection<MyMesh>::FaceFromVertexStrict(m);
  printf("%i vertexes are on the boundary and %i internal and internal faces are %i \n",borderNum,internalNum,faceNotOnBoundaryNum);
  int innerFacesNum = tri::UpdateSelection<MyMesh>::FaceErode(m);
  printf("Eroding the face selection of on ring of faces bring it to %i faces\n",innerFacesNum);
  
  
  return 0;
}

