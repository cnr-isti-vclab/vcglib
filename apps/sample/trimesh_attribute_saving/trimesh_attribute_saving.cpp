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
/*! \file trimesh_attribute.cpp
\ingroup code_sample

\brief the minimal example of using the attributes

Attributes are a simple mechanism to associate user-defined 'attributes' to the simplicies and to the mesh.
See the page '\ref attributes' for more details.
*/

#include<vcg/complex/complex.h>
#include<vcg/complex/algorithms/create/platonic.h>
#include<wrap/io_trimesh/export_ply.h>
class MyEdge;
class MyFace;
class MyVertex;
struct MyUsedTypes : public vcg::UsedTypes<	vcg::Use<MyVertex>::AsVertexType,
                                            vcg::Use<MyFace>  ::AsFaceType>{};

class MyVertex  : public vcg::Vertex< MyUsedTypes, vcg::vertex::Coord3f,vcg::vertex::Normal3f, vcg::vertex::BitFlags>{};
class MyFace    : public vcg::Face< MyUsedTypes, vcg::face::VertexRef, vcg::face::Normal3f, vcg::face::BitFlags> {};

class MyMesh : public vcg::tri::TriMesh< std::vector<MyVertex>, std::vector<MyFace> > {};
int main()
{
  MyMesh m;
  Torus<MyMesh>(m, 3.0f, 1.0f);
  //! [Adding a few attributes]
  // add a per-vertex attribute with type float named "GaussianCurvature"
  MyMesh::PerVertexAttributeHandle<float> 
      hvf = vcg::tri::Allocator<MyMesh>:: GetPerVertexAttribute<float> (m,std::string("GaussianCurvature"));
  
  MyMesh::PerVertexAttributeHandle<vcg::Point3f> 
      hv3f = vcg::tri::Allocator<MyMesh>:: GetPerVertexAttribute<vcg::Point3f> (m,std::string("InvertedNormal"));
  
  // add a per-face attribute with type float named "FaceArea"
  MyMesh::PerFaceAttributeHandle<float> 
      hff = vcg::tri::Allocator<MyMesh>:: GetPerFaceAttribute<float> (m,std::string("FaceArea"));
  //! [filling the attribute]
  vcg::tri::Allocator<MyMesh>::ClearPerVertexAttribute<float>(m,hvf, float(M_PI*2));
  vcg::tri::Allocator<MyMesh>::ClearPerVertexAttribute<vcg::Point3f>(m,hv3f, vcg::Point3f(0,0,0));
  
  ForEachFace(m, [&](MyFace &f){
    hff[&f]=vcg::DoubleArea(f)*0.5f;
    for(int i=0;i<3;++i){
      hvf[f.V(i)] -= vcg::Angle(f.P1(i)-f.P0(i),f.P2(i)-f.P0(i)); 
      hv3f[f.V(i)] -= vcg::NormalizedTriangleNormal(f);
    }
  }); 
  
  //! [Saving 3 attributes in ply, one of the 3 disguised as quality]
  vcg::tri::io::PlyInfo pi;
  pi.AddPerVertexFloatAttribute("GaussianCurvature","quality");
  pi.AddPerFaceFloatAttribute("FaceArea");
  pi.AddPerVertexPoint3fAttribute(m,"InvertedNormal");
  vcg::tri::io::ExporterPLY<MyMesh>::Save(m,"MeshWithCurvature.ply",false,pi);     
}
