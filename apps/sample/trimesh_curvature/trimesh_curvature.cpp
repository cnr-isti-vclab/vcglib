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
/*! \file trimesh_curvature.cpp
\ingroup code_sample

\brief an example showing the various techniques for computing curvatures

*/
#include <vcg/complex/complex.h>

#include <vcg/complex/algorithms/create/platonic.h>
#include <vcg/complex/algorithms/update/curvature.h>

#include<wrap/io_trimesh/export_off.h>

class MyEdge;
class MyFace;
class MyVertex;
struct MyUsedTypes : public vcg::UsedTypes<	vcg::Use<MyVertex>   ::AsVertexType,
                                            vcg::Use<MyEdge>     ::AsEdgeType,
                                            vcg::Use<MyFace>     ::AsFaceType>{};

class MyVertex  : public vcg::Vertex<MyUsedTypes, vcg::vertex::Coord3f, vcg::vertex::Normal3f, vcg::vertex::VFAdj, vcg::vertex::CurvatureDirf, vcg::vertex::Curvaturef, vcg::vertex::BitFlags  >{};
class MyFace    : public vcg::Face< MyUsedTypes, vcg::face::FFAdj, vcg::face::VFAdj, vcg::face::VertexRef, vcg::face::BitFlags > {};
class MyEdge    : public vcg::Edge<MyUsedTypes>{};
class MyMesh    : public vcg::tri::TriMesh< std::vector<MyVertex>, std::vector<MyFace> , std::vector<MyEdge>  > {};

int main( int /*argc*/, char **/*argv*/ )
{
  MyMesh m;
  vcg::tri::Torus(m,30,10);
  
  vcg::tri::UpdateTopology<MyMesh>::FaceFace(m);
  vcg::tri::UpdateTopology<MyMesh>::VertexFace(m);
  
  // Two different techniques for computing Discrete Gaussian and Mean Curvature
  // they require the presence of the vertex::Curvature component
  vcg::tri::UpdateCurvature<MyMesh>::PerVertex(m);
  vcg::tri::UpdateCurvature<MyMesh>::MeanAndGaussian(m);
  
  // Two different techniques for computing Principal Curvature Directions 
  // they require the presence of the vertex::CurvatureDir component
  vcg::tri::UpdateCurvature<MyMesh>::PrincipalDirections(m);
  vcg::tri::UpdateCurvature<MyMesh>::PrincipalDirectionsNormalCycle(m);
  printf("Input mesh  vn:%i fn:%i\n",m.VN(),m.FN());

  return 0;
}
