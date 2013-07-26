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
#include<vcg/complex/algorithms/point_sampling.h>
#include<vcg/complex/algorithms/voronoi_clustering.h>


using namespace vcg;
using namespace std;

class MyEdge;
class MyFace;
class MyVertex;
struct MyUsedTypes : public UsedTypes<	Use<MyVertex>   ::AsVertexType,
                                        Use<MyEdge>     ::AsEdgeType,
                                        Use<MyFace>     ::AsFaceType>{};

class MyVertex  : public Vertex<MyUsedTypes,  vertex::Coord3f, vertex::Normal3f, vertex::VFAdj , vertex::Qualityf, vertex::Color4b, vertex::BitFlags  >{};
class MyFace    : public Face< MyUsedTypes,   face::VertexRef, face::BitFlags, face::VFAdj > {};
//class MyEdge    : public Edge< MyUsedTypes> {};
class MyEdge    : public Edge< MyUsedTypes, edge::VertexRef, edge::BitFlags>{};
class MyMesh    : public tri::TriMesh< vector<MyVertex>, vector<MyEdge>, vector<MyFace>   > {};

class EmEdge;
class EmFace;
class EmVertex;
struct EmUsedTypes : public UsedTypes<	Use<EmVertex>   ::AsVertexType,
                                        Use<EmEdge>     ::AsEdgeType,
                                        Use<EmFace>     ::AsFaceType>{};

class EmVertex  : public Vertex<EmUsedTypes,  vertex::Coord3f, vertex::Normal3f, vertex::VFAdj , vertex::Qualityf, vertex::Color4b, vertex::BitFlags  >{};
class EmFace    : public Face< EmUsedTypes,   face::VertexRef, face::BitFlags, face::VFAdj > {};
class EmEdge    : public Edge< EmUsedTypes, edge::VertexRef> {};
//class EmEdge    : public Edge< EmUsedTypes, edge::VertexRef, edge::BitFlags>{};
class EmMesh    : public tri::TriMesh< vector<EmVertex>, vector<EmEdge>, vector<EmFace>   > {};


int main( int argc, char **argv )
{
  MyMesh baseMesh, outMesh, polyMesh;
  if(argc < 4 )
  {
    printf("Usage trimesh_voronoisampling mesh sampleNum variance iterNum\n");
     return -1;
  }
  int sampleNum = atoi(argv[2]);
  int iterNum   = atoi(argv[3]);
  float radiusVariance = atof(argv[4]);

  printf("Reading %s and sampling %i points with %i iteration and using %f variance\n",argv[1],sampleNum,iterNum,radiusVariance);
  int ret= tri::io::ImporterPLY<MyMesh>::Open(baseMesh,argv[1]);
  if(ret!=0)
  {
    printf("Unable to open %s for '%s'\n",argv[1],tri::io::ImporterPLY<MyMesh>::ErrorMsg(ret));
    return -1;
  }

  tri::UpdateTopology<MyMesh>::VertexFace(baseMesh);
  std::vector<MyVertex *> seedVec;
  vector<Point3f> pointVec;
  float radius=0;
  tri::VoronoiProcessingParameter vpp;
  vpp.deleteUnreachedRegionFlag=true;


  tri::PoissonSampling<MyMesh>(baseMesh,pointVec,sampleNum,radius);
  tri::VoronoiProcessing<MyMesh>::SeedToVertexConversion(baseMesh,pointVec,seedVec);

  tri::EuclideanDistance<MyMesh> dd;
  tri::VoronoiProcessing<MyMesh, tri::EuclideanDistance<MyMesh> >::VoronoiRelaxing(baseMesh, seedVec, iterNum, dd, vpp);
  tri::VoronoiProcessing<MyMesh, tri::EuclideanDistance<MyMesh> >::ConvertVoronoiDiagramToMesh(baseMesh,outMesh,polyMesh,seedVec, dd, vpp);
  tri::io::ExporterPLY<MyMesh>::Save(outMesh,"out.ply",tri::io::Mask::IOM_VERTCOLOR );
  tri::io::ExporterPLY<MyMesh>::Save(polyMesh,"poly.ply",tri::io::Mask::IOM_VERTCOLOR| tri::io::Mask::IOM_EDGEINDEX ,false);

  tri::io::ImporterPLY<MyMesh>::Open(baseMesh,argv[1]);
  tri::UpdateTopology<MyMesh>::VertexFace(baseMesh);
  tri::PoissonSampling<MyMesh>(baseMesh,pointVec,sampleNum,radius,radiusVariance);
  tri::VoronoiProcessing<MyMesh>::SeedToVertexConversion(baseMesh,pointVec,seedVec);
  tri::IsotropicDistance<MyMesh> id(baseMesh,radiusVariance);
  tri::VoronoiProcessing<MyMesh, tri::IsotropicDistance<MyMesh> >::VoronoiRelaxing(baseMesh, seedVec, iterNum,id,vpp);
  tri::VoronoiProcessing<MyMesh, tri::IsotropicDistance<MyMesh> >::ConvertVoronoiDiagramToMesh(baseMesh,outMesh,polyMesh,seedVec, id, vpp);

  tri::io::ExporterPLY<MyMesh>::Save(outMesh,"outW.ply",tri::io::Mask::IOM_VERTCOLOR );
  tri::io::ExporterPLY<MyMesh>::Save(polyMesh,"polyW.ply",tri::io::Mask::IOM_VERTCOLOR | tri::io::Mask::IOM_EDGEINDEX,false);
  return 0;
}
