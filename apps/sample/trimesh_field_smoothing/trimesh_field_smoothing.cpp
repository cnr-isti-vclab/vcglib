/****************************************************************************
* VCGLib                                                            o o     *
* Visual and Computer Graphics Library                            o     o   *
*                                                                _   O  _   *
* Copyright(C) 2004-2009                                           \/)\/    *
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
#include<vcg/complex/algorithms/create/platonic.h>

#include<wrap/io_trimesh/import_ply.h>
#include<wrap/io_trimesh/export_ply.h>
#include<vcg/complex/algorithms/create/platonic.h>
#include<vcg/complex/algorithms/smooth_field.h>

using namespace vcg;
using namespace std;

class MyEdge;
class MyFace;
class MyVertex;
struct MyUsedTypes : public UsedTypes<	Use<MyVertex>   ::AsVertexType,
                                        Use<MyEdge>     ::AsEdgeType,
                                        Use<MyFace>     ::AsFaceType>{};

class MyVertex  : public Vertex<MyUsedTypes, vertex::Coord3f, vertex::Normal3f, vertex::VFAdj, vertex::Qualityf, vertex::Color4b, vertex::BitFlags, vertex::CurvatureDirf  >{};
class MyFace    : public Face< MyUsedTypes, face::VertexRef, face::Normal3f, face::BitFlags, face::FFAdj, face::VFAdj, face::CurvatureDirf, face::Qualityf, face::Color4b> {};
class MyEdge    : public Edge< MyUsedTypes>{};
class MyMesh    : public tri::TriMesh< std::vector<MyVertex>, std::vector<MyFace> , vector<MyEdge>  > {};


int main( int argc, char **argv )
{
  MyMesh m;
  tri::Hexahedron(m);
  tri::UpdateTopology<MyMesh>::FaceFace(m);
  tri::UpdateFlags<MyMesh>::FaceBorderFromFF(m);
  for(int i=0;i<3;++i)
    tri::Refine<MyMesh, tri::MidPoint<MyMesh> >(m,tri::MidPoint<MyMesh>(&m));

  tri::UpdateNormal<MyMesh>::PerVertexNormalized(m);
  tri::Smooth<MyMesh>::VertexCoordLaplacian(m,3);
  tri::FieldSmoother<MyMesh>::InitByCurvature(m);
  tri::FieldSmoother<MyMesh>::SmoothParam par;

  tri::FieldSmoother<MyMesh>::SmoothDirections(m,par);
  tri::io::ExporterPLY<MyMesh>::Save(m,"Full.ply",tri::io::Mask::IOM_VERTCOLOR|tri::io::Mask::IOM_WEDGTEXCOORD );
  return 0;
}
