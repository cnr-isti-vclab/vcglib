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
/*! \file trimesh_create.cpp
\ingroup code_sample

\brief Some examples of building meshes out of nothing

*/
#include<vcg/complex/complex.h>
#include<vcg/complex/algorithms/create/platonic.h>

// input output
#include <wrap/io_trimesh/export_off.h>

class MyFace;
class MyVertex;

struct MyUsedTypes : public vcg::UsedTypes<	vcg::Use<MyVertex>::AsVertexType,    vcg::Use<MyFace>::AsFaceType>{};

class MyVertex  : public vcg::Vertex< MyUsedTypes, vcg::vertex::Coord3f, vcg::vertex::Normal3f, vcg::vertex::Color4b, vcg::vertex::BitFlags  >{};
class MyFace    : public vcg::Face < MyUsedTypes, vcg::face::VertexRef, vcg::face::Normal3f, vcg::face::FFAdj, vcg::face::BitFlags > {};
class MyMesh    : public vcg::tri::TriMesh< std::vector<MyVertex>, std::vector<MyFace> > {};

using namespace std;
using namespace vcg;

int  main()
{
  MyMesh diskMesh;

  // Create a simple triangle mesh using just a vector of coords and a vector of indexes
  vector<Point3f> coordVec;
  vector<Point3i> indexVec;
  coordVec.push_back(Point3f(0,0,0));
  for(int i=0;i<36;++i) {
    float angleRad = float(i)*M_PI/18.0;
    coordVec.push_back(Point3f(sin(angleRad),cos(angleRad),0));
    indexVec.push_back(Point3i(0,i+1,1+(i+1)%36));
  }

  tri::BuildMeshFromCoordVectorIndexVector(diskMesh,coordVec,indexVec);
  tri::io::ExporterOFF<MyMesh>::Save(diskMesh,"disc.off");

  // Create the platonic solids  
  MyMesh platonicMesh;
  tri::Tetrahedron(platonicMesh);
  tri::io::ExporterOFF<MyMesh>::Save(platonicMesh,"tetrahedron.off");
  tri::Octahedron(platonicMesh);
  tri::io::ExporterOFF<MyMesh>::Save(platonicMesh,"octahedron.off");
  tri::Hexahedron(platonicMesh);
  tri::io::ExporterOFF<MyMesh>::Save(platonicMesh,"hexahedron.off");
  tri::Dodecahedron(platonicMesh);
  tri::io::ExporterOFF<MyMesh>::Save(platonicMesh,"dodecahedron.off");
  tri::Icosahedron(platonicMesh);
  tri::io::ExporterOFF<MyMesh>::Save(platonicMesh,"icosahedron.off");
  
  // Procedurally transform a mesh into a solid collection of triangular prisms
  MyMesh facePrismMesh;
  tri::BuildPrismFaceShell(platonicMesh, facePrismMesh, 0.1f, 0.1f);
  tri::io::ExporterOFF<MyMesh>::Save(facePrismMesh,"facePrism.off");
  
  return 0;
}
