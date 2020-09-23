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
#include<vcg/complex/complex.h>

// input output
#include <wrap/io_trimesh/import_ply.h>
#include <wrap/io_trimesh/export_ply.h>

#include<vcg/complex/algorithms/outline_support.h>

#include <vcg/space/outline2_packer.h>
#include <wrap/qt/outline2_rasterizer.h>

#include <vcg/space/rasterized_outline2_packer.h>
#include <wrap/qt/Outline2ToQImage.h>
/*! \file trimesh_texture.cpp
\ingroup code_sample

\brief a small example about packing texture regions

It takes a sample mesh as input, computes the connected
components of its parametrization and pack them again using two different strategies
*/

using namespace vcg;

class MyEdge;
class MyFace;
class MyVertex;
struct MyUsedTypes : public UsedTypes<	Use<MyVertex>::AsVertexType, Use<MyFace>::AsFaceType>{};

class MyVertex  : public Vertex< MyUsedTypes, vertex::Coord3f, vertex::BitFlags  >{};
class MyFace    : public Face  < MyUsedTypes, face::VertexRef, face::FFAdj, face::WedgeTexCoord2f, face::Mark, face::BitFlags > {};
class MyMesh : public tri::TriMesh< std::vector<MyVertex>, std::vector<MyFace > >{};


int main(int ,char ** )
{
  MyMesh m,tm;
  tri::io::ImporterPLY<MyMesh>::Open(m,"../../meshes/bunny10k_textured.ply");

  // 1) Build a mesh with tex coords as coords. 
  
  for( auto &&f : m.face)
  {
   tri::Allocator<MyMesh>::AddFace(tm,
                                   Point3f(f.WT(0).U(),f.WT(0).V(),0),
                                   Point3f(f.WT(1).U(),f.WT(1).V(),0),
                                   Point3f(f.WT(2).U(),f.WT(2).V(),0));
  }
  
  // 2) compute connected components (e.g. texture regions in an atlas)
  tri::Clean<MyMesh>::RemoveDuplicateVertex(tm);
  tri::UpdateTopology<MyMesh>::FaceFace(tm);
  std::vector<std::pair<int,MyMesh::FacePointer> > fpVec;
  tri::Clean<MyMesh>::ConnectedComponents(tm,fpVec);
  printf("Mesh has %lu texture components\n",fpVec.size());
  tri::io::ExporterPLY<MyMesh>::Save(tm,"out.ply");
  std::vector< std::vector<Point2f> > outline2Vec;
  
  // build the 2D outlines of each regions
  for(size_t i=0; i<fpVec.size();++i)
  {
    tri::UpdateSelection<MyMesh>::FaceClear(tm);
    fpVec[i].second->SetS();
    tri::UpdateSelection<MyMesh>::FaceConnectedFF(tm);
    tri::UpdateSelection<MyMesh>::VertexClear(tm);
    tri::UpdateSelection<MyMesh>::VertexFromFaceLoose(tm);

    MyMesh comp;
    tri::Append<MyMesh,MyMesh>::Mesh(comp, tm, true);

    std::vector< std::vector<Point3f> > outline3Vec;
    tri::OutlineUtil<float>::ConvertMeshBoundaryToOutline3Vec(comp, outline3Vec);
    std::vector< std::vector<Point2f> > compOutline2Vec;
    tri::OutlineUtil<float>::ConvertOutline3VecToOutline2Vec(outline3Vec,compOutline2Vec);
    int largestInd=tri::OutlineUtil<float>::LargestOutline2(compOutline2Vec);
    if(tri::OutlineUtil<float>::Outline2Area(compOutline2Vec[largestInd])<0)
      tri::OutlineUtil<float>::ReverseOutline2(compOutline2Vec[largestInd]);

    outline2Vec.push_back(compOutline2Vec[largestInd]);
  }

  printf("Mesh has %lu texture components\n",outline2Vec.size());

  Outline2Dumper::Param pp;
  Similarity2f sim;
  sim.sca=1024.0f;
  std::vector<Similarity2f> trVec(outline2Vec.size(),sim);
  printf("Mesh has %lu texture components\n",outline2Vec.size());
  
  // Dump the original parametrization as a png
  Outline2Dumper::dumpOutline2VecPNG("PrePack.png",outline2Vec,trVec,pp);

  // Pack using Axis Aligned Rect
  const Point2i containerSize(1024,1024);
  Point2f finalSize(1024,1024);
  PolyPacker<float>::PackAsAxisAlignedRect(outline2Vec,containerSize,trVec,finalSize);

  Outline2Dumper::dumpOutline2VecPNG("PostPack.png",outline2Vec,trVec,pp);

   // Pack using Oriented Rect
  RasterizedOutline2Packer<float, QtOutline2Rasterizer>::Parameters packingParam;
  packingParam.costFunction  = RasterizedOutline2Packer<float, QtOutline2Rasterizer>::Parameters::LowestHorizon;
  packingParam.doubleHorizon = true;
  packingParam.innerHorizon = true;
  packingParam.permutations = false;
  packingParam.rotationNum = 16; //number of rasterizations in 90°
  packingParam.gutterWidth = 2;

  RasterizedOutline2Packer<float, QtOutline2Rasterizer>::Pack(outline2Vec,containerSize,trVec,packingParam);
  Outline2Dumper::dumpOutline2VecPNG("PostPackRR.png",outline2Vec,trVec,pp);

  return 0;
}

