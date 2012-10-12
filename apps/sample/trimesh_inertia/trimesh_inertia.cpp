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
#include<vcg/complex/complex.h>

#include<wrap/io_trimesh/import_off.h>

#include<vcg/complex/algorithms/inertia.h>
#include<vcg/complex/algorithms/create/platonic.h>

using namespace vcg;
using namespace std;

class MyEdge;
class MyFace;
class MyVertex;
struct MyUsedTypes : public UsedTypes<	Use<MyVertex>   ::AsVertexType,
                                        Use<MyEdge>     ::AsEdgeType,
                                        Use<MyFace>     ::AsFaceType>{};

class MyVertex  : public Vertex<MyUsedTypes, vertex::Coord3f, vertex::Normal3f, vertex::BitFlags  >{};
class MyFace    : public Face< MyUsedTypes, face::FFAdj, face::Normal3f, face::VertexRef, face::BitFlags > {};
class MyEdge    : public Edge<MyUsedTypes>{};
class MyMesh    : public tri::TriMesh< vector<MyVertex>, vector<MyFace> , vector<MyEdge>  > {};

int main( int argc, char **argv )
{
  MyMesh tet,oct,hex,dod,ico;

  tri::Hexahedron(hex);
  tri::Tetrahedron(tet);
  tri::Octahedron(oct);
  tri::Dodecahedron(dod);
  tri::Icosahedron(ico);
  Matrix44f ScaleM,TransM;
  ScaleM.SetScale(1,2,1);
  TransM.SetTranslate(1,1,1);
//  tri::UpdatePosition<MyMesh>::Matrix(hex,ScaleM);
  tri::UpdatePosition<MyMesh>::Matrix(hex,TransM);

  tri::Inertia<MyMesh> I;
  I.Compute(hex);
  Point3f cc = I.CenterOfMass();
  printf("Mass %f \n",I.Mass());
  printf("CenterOfMass %f %f %f\n",cc[0],cc[1],cc[2]);
  Matrix33f IT;
  Point3f ITv;
  I.InertiaTensorEigen(IT,ITv);
  printf("InertiaTensor  %f %f %f\n\n",ITv[0],ITv[1],ITv[2]);

  printf("InertiaTensor  %f %f %f\n",IT[0][0],IT[0][1],IT[0][2]);
  printf("InertiaTensor  %f %f %f\n",IT[1][0],IT[1][1],IT[1][2]);
  printf("InertiaTensor  %f %f %f\n",IT[2][0],IT[2][1],IT[2][2]);

  return 0;
}
