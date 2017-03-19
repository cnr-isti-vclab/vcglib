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
#include<vcg/complex/algorithms/voronoi_processing.h>
#include<vcg/complex/algorithms/voronoi_volume_sampling.h>


using namespace vcg;
using namespace std;

class MyEdge;
class MyFace;
class MyVertex;
struct MyUsedTypes : public UsedTypes<	Use<MyVertex>   ::AsVertexType,
                                        Use<MyEdge>     ::AsEdgeType,
                                        Use<MyFace>     ::AsFaceType>{};

class MyVertex  : public Vertex<MyUsedTypes,  vertex::Coord3f, vertex::Normal3f, vertex::VFAdj, vertex::Qualityf, vertex::Color4b, vertex::BitFlags  >{};
class MyFace    : public Face< MyUsedTypes,   face::VertexRef, face::Normal3f, face::BitFlags, face::Mark, face::VFAdj, face::FFAdj > {};
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
class EmMesh    : public tri::TriMesh< vector<EmVertex>, vector<EmEdge>, vector<EmFace>   > {};


int main( int argc, char **argv )
{
  MyMesh mSur; // the original surface
  MyMesh mOff; // the offsetted surface
  if(argc < 3 )
  {
    printf("Usage trimesh_voronoisampling origSurf offsetSurf radius \n"
           "radius is expressed as a perc bbox diag  (e.g. 0.1 -> 10 % of bbox diag) \n");
     return -1;
  }
  
  tri::io::ImporterPLY<MyMesh>::Open(mSur,argv[1]);
  tri::io::ImporterPLY<MyMesh>::Open(mOff,argv[2]);
  tri::UpdateNormal<MyMesh>::PerVertexNormalizedPerFaceNormalized(mOff);
  tri::UpdateNormal<MyMesh>::PerVertexNormalizedPerFaceNormalized(mSur);
  
  std::vector<Point3f> sampleSurVec, sampleOffVec;
//  MontecarloSampling(mSur, sampleSurVec,50000);
  MontecarloSampling(mOff, sampleOffVec,50000);
  sampleSurVec.insert( sampleSurVec.end(), sampleOffVec.begin(), sampleOffVec.end() );
  
  
  printf("Read %i vn %i fn \n",mOff.vn, mOff.fn);
  
  float poissonRadius = mOff.bbox.Diag()*atof(argv[3]); 
//  float poissonRadius = mOff.bbox.Diag()/10; 
  printf("Poisson Radius %f\n",poissonRadius);
  
  float sampleSurfRadius=mOff.bbox.Diag()/100.0f;
  int montecarloSampleNum = 100000;
 
  MyMesh  seedM; 
  VoronoiVolumeSampling<MyMesh> vvs(mOff, seedM);
  printf("Sampling Surface at a radius %f ",sampleSurfRadius);
  vvs.Init(sampleSurfRadius);
  tri::BuildMeshFromCoordVector(vvs.seedDomainMesh, sampleSurVec);
  printf("Sampled\n");
  vvs.BuildVolumeSampling(montecarloSampleNum,0,poissonRadius);

  
  tri::io::ExporterPLY<MyMesh>::Save(vvs.seedDomainMesh,"seedDomainMesh.ply");
  tri::io::ExporterPLY<MyMesh>::Save(vvs.poissonSurfaceMesh,"poissonSurfaceMesh.ply");
  tri::io::ExporterPLY<MyMesh>::Save(vvs.montecarloVolumeMesh,"montecarloVolumeMesh.ply");
  tri::io::ExporterPLY<MyMesh>::Save(seedM,"seedMesh0.ply");
//  vvs.restrictedRelaxationFlag=true;
//  vvs.BarycentricRelaxVoronoiSamples(10);
  vvs.QuadricRelaxVoronoiSamples(10);
  
  tri::UpdateColor<MyMesh>::PerVertexQualityRamp(seedM);
  tri::io::ExporterPLY<MyMesh>::Save(seedM,"seedMesh1.ply",tri::io::Mask::IOM_VERTCOLOR | tri::io::Mask::IOM_VERTQUALITY );
  vvs.QuadricRelaxVoronoiSamples(10);
  
  tri::UpdateColor<MyMesh>::PerVertexQualityRamp(seedM);
  tri::io::ExporterPLY<MyMesh>::Save(seedM,"seedMesh2.ply",tri::io::Mask::IOM_VERTCOLOR | tri::io::Mask::IOM_VERTQUALITY );
  printf("\n Saved %i points \n",seedM.vn);
  // Second Pipeline 
  return true;
  
  return 0;
}
