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

// input output
#include <wrap/io_trimesh/import_ply.h>
#include <wrap/io_trimesh/export_ply.h>

#include<vcg/complex/algorithms/update/topology.h>
#include<vcg/complex/algorithms/outline_support.h>

#include<vcg/simplex/face/pos.h>
#include <vcg/space/poly_packer.h>
#include <wrap/qt/PolyToQImage.h>
//#include<vcg/../../sandbox/bernabeia/polypacker/new_polypacker.h>
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

  tri::Allocator<MyMesh>::AddVertices(tm,m.fn*3);
  tri::Allocator<MyMesh>::AddFaces(tm,m.fn);

  for(size_t i=0;i<m.fn;++i)
  {
    for(int j=0;j<3;++j)
    {
      tm.face[i].V(j)=&tm.vert[i*3+j];
      tm.vert[i*3+j].P()[0] = m.face[i].WT(j).U();
      tm.vert[i*3+j].P()[1] = m.face[i].WT(j).V();
      tm.vert[i*3+j].P()[2] = 0;
    }
  }
  tri::Clean<MyMesh>::RemoveDuplicateVertex(tm);
  std::vector<std::pair<int,MyMesh::FacePointer> > fpVec;
  tri::UpdateTopology<MyMesh>::FaceFace(tm);
  tri::Clean<MyMesh>::ConnectedComponents(tm,fpVec);
  printf("Mesh has %i texture components\n",fpVec.size());
  tri::io::ExporterPLY<MyMesh>::Save(tm,"out.ply");

  std::vector< std::vector<Point3f> > outline3Vec;
  std::vector< std::vector<Point2f> > outline2Vec;
  tri::OutlineUtil<MyMesh>::ConvertMeshBoundaryToOutline3Vec(tm, outline3Vec);
  printf("Mesh has %i texture components\n",outline3Vec.size());
  tri::OutlineUtil<MyMesh>::ConvertOutline3VeToOutline2Vec(outline3Vec,outline2Vec);
  PolyDumper::Param pp;
  Similarity2f sim;
  sim.sca=1024.0f;
  std::vector<Similarity2f> trVec(outline2Vec.size(),sim);
  printf("Mesh has %i texture components\n",outline2Vec.size());
  PolyDumper::dumpOutline2VecPNG("PrePack.png",outline2Vec,trVec,pp);

  const Point2f containerSize(1024,1024);
  Point2f finalSize(1024,1024);
  PolyPacker<float>::PackAsAxisAlignedRect(outline2Vec,containerSize,trVec,finalSize);

  PolyDumper::dumpOutline2VecPNG("PostPack.png",outline2Vec,trVec,pp);

  return 0;
}

