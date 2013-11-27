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
#include<wrap/io_trimesh/export_off.h>
#include<wrap/io_trimesh/export_ply.h>
#include<wrap/io_trimesh/export_dxf.h>
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

class MyVertex  : public Vertex<MyUsedTypes,  vertex::Coord3f, vertex::Normal3f, vertex::VFAdj, vertex::Qualityf, vertex::Color4b, vertex::BitFlags  >{};
class MyFace    : public Face< MyUsedTypes,   face::VertexRef, face::Normal3f, face::BitFlags, face::VFAdj, face::FFAdj > {};
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
    printf("Usage trimesh_voronoisampling mesh sampleNum iterNum edgeCollapsePerc \n");
     return -1;
  }
  int sampleNum = atoi(argv[2]);
  int iterNum   = atoi(argv[3]);
  float collapseShortEdgePerc = atof(argv[4]);

  printf("Reading %s and sampling %i points with %i iteration and using %f variance\n",argv[1],sampleNum,iterNum,collapseShortEdgePerc);
  int ret= tri::io::ImporterPLY<MyMesh>::Open(baseMesh,argv[1]);
  if(ret!=0)
  {
    printf("Unable to open %s for '%s'\n",argv[1],tri::io::ImporterPLY<MyMesh>::ErrorMsg(ret));
    return -1;
  }

  tri::UpdateTopology<MyMesh>::VertexFace(baseMesh);
  tri::UpdateFlags<MyMesh>::FaceBorderFromVF(baseMesh);

  // -- Build the mesh with corners
  MyMesh cornerMesh;
  std::vector<Point3f> sampleVec;
  tri::TrivialSampler<MyMesh> mps(sampleVec);
  tri::SurfaceSampling<MyMesh,tri::TrivialSampler<MyMesh> >::VertexBorderCorner(baseMesh,mps,math::ToRad(150.f));
  tri::Build(cornerMesh,sampleVec);

  // -- Build the montercarlo sampling of the surface
  MyMesh MontecarloSurfaceMesh;
  sampleVec.clear();
  tri::SurfaceSampling<MyMesh,tri::TrivialSampler<MyMesh> >::Montecarlo(baseMesh,mps,50000);
  tri::Build(MontecarloSurfaceMesh,sampleVec);
  tri::io::ExporterPLY<MyMesh>::Save(MontecarloSurfaceMesh,"MontecarloSurfaceMesh.ply");

  // -- Prune the montecarlo sampling with poisson strategy using the precomputed corner vertexes.
  tri::SurfaceSampling<MyMesh,tri::TrivialSampler<MyMesh> >::PoissonDiskParam pp;
  pp.preGenMesh = &cornerMesh;
  pp.preGenFlag=true;
  sampleVec.clear();
  float radius = tri::SurfaceSampling<MyMesh,tri::TrivialSampler<MyMesh> >::ComputePoissonDiskRadius(baseMesh,sampleNum);
  tri::SurfaceSampling<MyMesh,tri::TrivialSampler<MyMesh> >::PoissonDiskPruning(mps, MontecarloSurfaceMesh, radius, pp);
  MyMesh PoissonMesh;
  tri::Build(PoissonMesh,sampleVec);
  tri::io::ExporterPLY<MyMesh>::Save(PoissonMesh,"PoissonMesh.ply");

  std::vector<MyVertex *> seedVec;
  tri::VoronoiProcessing<MyMesh>::SeedToVertexConversion(baseMesh,sampleVec,seedVec);
  float eps = baseMesh.bbox.Diag()/10000.0f;
  for(size_t i=0;i<cornerMesh.vert.size();++i)
  {
    for(size_t j=0;j<seedVec.size();++j)
      if(Distance(cornerMesh.vert[i].P(),seedVec[j]->P()) < eps)
        seedVec[j]->SetS();
  }

  tri::VoronoiProcessingParameter vpp;
  vpp.deleteUnreachedRegionFlag=true;
  vpp.fixSelectedSeed=true;
  vpp.collapseShortEdge=true;
  vpp.collapseShortEdgePerc=collapseShortEdgePerc;
  vpp.triangulateRegion = true;
  vpp.unbiasedSeedFlag =true;

  tri::EuclideanDistance<MyMesh> dd;
  int t0=clock();
  tri::VoronoiProcessing<MyMesh, tri::EuclideanDistance<MyMesh> >::VoronoiRelaxing(baseMesh, seedVec, iterNum, dd, vpp);
  int t1=clock();
  tri::VoronoiProcessing<MyMesh, tri::EuclideanDistance<MyMesh> >::ConvertVoronoiDiagramToMesh(baseMesh,outMesh,polyMesh, seedVec, dd, vpp);
  tri::io::ExporterPLY<MyMesh>::Save(baseMesh,"base.ply",tri::io::Mask::IOM_VERTCOLOR + tri::io::Mask::IOM_VERTQUALITY );
  tri::io::ExporterPLY<MyMesh>::Save(outMesh,"out.ply",tri::io::Mask::IOM_VERTCOLOR + tri::io::Mask::IOM_FLAGS );
  tri::io::ExporterPLY<MyMesh>::Save(polyMesh,"poly.ply",tri::io::Mask::IOM_VERTCOLOR| tri::io::Mask::IOM_EDGEINDEX ,false);

//  tri::io::ImporterPLY<MyMesh>::Open(baseMesh,argv[1]);
//  tri::UpdateTopology<MyMesh>::VertexFace(baseMesh);
//  tri::PoissonSampling<MyMesh>(baseMesh,pointVec,sampleNum,radius,radiusVariance);
//  tri::VoronoiProcessing<MyMesh>::SeedToVertexConversion(baseMesh,pointVec,seedVec);
//  tri::IsotropicDistance<MyMesh> id(baseMesh,radiusVariance);
//  tri::VoronoiProcessing<MyMesh, tri::IsotropicDistance<MyMesh> >::VoronoiRelaxing(baseMesh, seedVec, iterNum,id,vpp);
//  tri::VoronoiProcessing<MyMesh, tri::IsotropicDistance<MyMesh> >::ConvertVoronoiDiagramToMesh(baseMesh,outMesh,polyMesh,seedVec, id, vpp);

//  tri::io::ExporterPLY<MyMesh>::Save(outMesh,"outW.ply",tri::io::Mask::IOM_VERTCOLOR );
//  tri::io::ExporterPLY<MyMesh>::Save(polyMesh,"polyW.ply",tri::io::Mask::IOM_VERTCOLOR | tri::io::Mask::IOM_EDGEINDEX,false);
//  tri::io::ExporterDXF<MyMesh>::Save(polyMesh,"outW.dxf");
  printf("Completed! %i iterations in %f sec for %lu seeds \n",iterNum,float(t1-t0)/CLOCKS_PER_SEC,seedVec.size());
  return 0;
}
