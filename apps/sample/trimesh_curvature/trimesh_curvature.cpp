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

#include<wrap/io_trimesh/export_off.h>
#include <vcg/complex/algorithms/create/platonic.h>

#include<vcg/complex/algorithms/update/curvature.h>
#include<vcg/complex/algorithms/update/normal.h>

class MyEdge;
class MyFace;
class MyVertex;
struct MyUsedTypes : public UsedTypes<	Use<MyVertex>   ::AsVertexType,
                                        Use<MyEdge>     ::AsEdgeType,
                                        Use<MyFace>     ::AsFaceType>{};

class MyVertex  : public Vertex<MyUsedTypes, vertex::Coord3f, vertex::Normal3f, vertex::BitFlags  >{};
class MyFace    : public Face< MyUsedTypes, face::FFAdj,  face::VertexRef, face::BitFlags > {};
class MyEdge    : public Edge<MyUsedTypes>{};
class MyMesh    : public vcg::tri::TriMesh< std::vector<MyVertex>, std::vector<MyFace> , std::vector<MyEdge>  > {};

int main( int argc, char **argv )
{
  MyMesh m;
  vcg::tri::Torus(m,30,10);
  vcg::tri::io::ExporterOFF<MyMesh>::Save(m,"torus.off");

  vcg::tri::UpdateTopology<MyMesh>::FaceFace(m);
  vcg::tri::UpdateCurvature<MyMesh>::PerVertex(m);
  vcg::tri::UpdateCurvature<MyMesh>::MeanAndGaussian(m);
  vcg::tri::UpdateCurvature<MyMesh>::PrincipalDirections(m);
  vcg::tri::UpdateCurvature<MyMesh>::PrincipalDirectionsNormalCycle(m);
  vcg::tri::UpdateCurvature<MyMesh>::PrincipalDirectionsPCA(m,m.bbox.Diag()/100);

  vcg::tri::UpdateNormal<MyMesh>::PerVertexNormalized(m);
  printf("Input mesh  vn:%i fn:%i\n",m.VN(),m.FN());
  printf( "Mesh has %i vert and %i faces\n", m.VN(), m.FN() );

  return 0;
}
