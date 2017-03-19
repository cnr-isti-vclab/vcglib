/****************************************************************************
* VCGLib                                                            o o     *
* Visual and Computer Graphics Library                            o     o   *
*                                                                _   O  _   *
* Copyright(C) 2004-2012                                           \/)\/    *
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
/*! \file trimesh_normal.cpp
\ingroup code_sample

\brief An example of all the methods for computing normals over a mesh.

*/
#include<vcg/complex/complex.h>

#include<wrap/io_trimesh/import_ply.h>

#include<vcg/complex/algorithms/clean.h>
#include<vcg/complex/algorithms/update/bounding.h>

using namespace vcg;
using namespace std;

class MyEdge;
class MyFace;
class MyVertex;
struct MyUsedTypes : public UsedTypes<	Use<MyVertex>   ::AsVertexType,
                                        Use<MyEdge>     ::AsEdgeType,
                                        Use<MyFace>     ::AsFaceType>{};

class MyVertex  : public Vertex<MyUsedTypes, vertex::Coord3f, vertex::Normal3f, vertex::BitFlags  >{};
class MyFace    : public Face< MyUsedTypes, face::Mark,  face::VertexRef, face::Normal3f, face::BitFlags > {};
class MyEdge    : public Edge<MyUsedTypes>{};
class MyMesh    : public tri::TriMesh< vector<MyVertex>, vector<MyFace> , vector<MyEdge>  > {};

int main( int argc, char **argv )
{
  MyMesh m;

  if(tri::io::ImporterPLY<MyMesh>::Open(m,"../../meshes/SelfIntersection2.ply")!=0)
  {
    printf("Error reading file  %s\n",argv[1]);
    exit(0);
  }
  std::vector<MyFace *> retVec;
  tri::UpdateNormal<MyMesh>::PerVertexNormalizedPerFaceNormalized(m);
  tri::UpdateBounding<MyMesh>::Box(m);
  tri::Clean<MyMesh>::SelfIntersections(m, retVec);
  printf("Mesh has %i selfintersecting faces\n",retVec.size());
  return 0;
}
