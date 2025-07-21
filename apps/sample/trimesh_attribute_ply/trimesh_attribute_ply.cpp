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

\brief the minimal example of saving and loading attributes into PLY

Attributes are a simple mechanism to associate user-defined 'attributes' to the simplicies and to the mesh.
See the page '\ref attributes' for more details.
*/

#include<vcg/complex/complex.h>
#include<vcg/complex/algorithms/create/platonic.h>
#include<wrap/io_trimesh/import_ply.h>
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
  vcg::tri::Torus<MyMesh>(m, 3.0f, 1.0f);
  //! [Adding a few attributes]
  // add a per-vertex attribute with type float named "GaussianCurvature"
  MyMesh::PerVertexAttributeHandle<float> 
      hvf = vcg::tri::Allocator<MyMesh>:: GetPerVertexAttribute<float> (m,std::string("GaussianCurvature"));
  
  MyMesh::PerVertexAttributeHandle<vcg::Point3f> 
      hv3f = vcg::tri::Allocator<MyMesh>:: GetPerVertexAttribute<vcg::Point3f> (m,std::string("InvertedNormal"));
  
  // add a per-face attribute with type float named "FaceArea"
  MyMesh::PerFaceAttributeHandle<float> 
      hff = vcg::tri::Allocator<MyMesh>:: GetPerFaceAttribute<float> (m,std::string("FaceArea"));
  //! [filling the attributes]
  vcg::tri::Allocator<MyMesh>::ClearPerVertexAttribute<float>(m,hvf, float(M_PI*2));
  vcg::tri::Allocator<MyMesh>::ClearPerVertexAttribute<vcg::Point3f>(m,hv3f, vcg::Point3f(0,0,0));
  
  ForEachFace(m, [&](MyFace &f){
    hff[&f]=vcg::DoubleArea(f)*0.5f;
    for(int i=0;i<3;++i){
      hvf[f.V(i)] -= vcg::Angle(f.P1(i)-f.P0(i),f.P2(i)-f.P0(i)); 
      hv3f[f.V(i)] -= vcg::NormalizedTriangleNormal(f);
    }
  }); 
  
  //! [Saving 3 attributes in ply, one of the 3 is saved as the standard ply property quality]
  vcg::tri::io::PlyInfo pi;
  
  // // The first attribute is a per-vertex float attribute named "GaussianCurvature" and is saved as quality even if this mesh has not quality
  // pi.AddPerVertexFloatAttribute("GaussianCurvature","quality");
  pi.AddPerVertexFloatAttribute("GaussianCurvature");
  pi.AddPerFaceFloatAttribute("FaceArea");
  pi.AddPerVertexPoint3fAttribute(m,"InvertedNormal");
  vcg::tri::io::ExporterPLY<MyMesh>::Save(m,"MeshWithCurvature.ply",false,pi);
  
  printf("Face 0 has area %f\n",hff[&m.face[0]]);
  printf("Vertex 1 has GaussianCurvature %f\n",hvf[1]);
  
  /******************************************************************************/
  
  //! Now Trying to reload it in a new mesh
  MyMesh mi;
  //! use an attribute with a slightly different name for testing
  MyMesh::PerFaceAttributeHandle<float> 
      ihff = vcg::tri::Allocator<MyMesh>:: GetPerFaceAttribute<float> (mi,std::string("MyFaceArea"));
  
  MyMesh::PerVertexAttributeHandle<float> 
      ihvf = vcg::tri::Allocator<MyMesh>:: GetPerVertexAttribute<float> (mi,std::string("GaussianCurvature"));
  
  vcg::tri::io::PlyInfo piOpen0;
  // Use load Mask to get info about the file
  int mask=0;
  vcg::tri::io::ImporterPLY<MyMesh>::LoadMask("MeshWithCurvature.ply",mask,piOpen0);
  
  
  
  vcg::tri::io::PlyInfo piOpen;
  piOpen.AddPerFaceFloatAttribute("MyFaceArea", "FaceArea");
  piOpen.AddPerVertexFloatAttribute("GaussianCurvature");
  piOpen.AddPerVertexPoint3fAttribute(mi,"InvertedNormal");
  
  
  int ret = vcg::tri::io::ImporterPLY<MyMesh>::Open(mi,"MeshWithCurvature.ply",piOpen);
  if(ret!=0)
  {
    printf("Unable to open MeshWithCurvature.ply for '%s'\n",vcg::tri::io::ImporterPLY<MyMesh>::ErrorMsg(ret));
    return -1;
  }
  
  printf("Loaded %i vertices and %i faces\n",mi.VN(),mi.FN());
  printf("Face 0 has area %f\n",ihff[&mi.face[0]]);
  printf("Vertex 1 has GaussianCurvature %f\n",ihvf[1]);
  
  
}
