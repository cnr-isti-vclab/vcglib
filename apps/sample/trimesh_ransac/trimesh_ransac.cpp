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
#include<wrap/io_trimesh/import.h>
#include<wrap/io_trimesh/export.h>
#include<vcg/complex/algorithms/update/topology.h>
#include<vcg/complex/algorithms/update/normal.h>
#include<vcg/complex/algorithms/update/color.h>
#include<vcg/complex/algorithms/inertia.h>
#include<vcg/complex/algorithms/point_sampling.h>
#include<vcg/complex/algorithms/ransac_matching.h>

#include<vcg/space/index/kdtree/kdtree.h>
#include<vcg/space/point_matching.h>

using namespace vcg;

class MyVertex; class MyEdge; class MyFace;
struct MyUsedTypes : public UsedTypes<Use<MyVertex>   ::AsVertexType,
                                           Use<MyEdge>     ::AsEdgeType,
                                           Use<MyFace>     ::AsFaceType>{};

class MyVertex  : public Vertex< MyUsedTypes, vertex::Coord3f, vertex::Normal3f, vertex::Color4b, vertex::Qualityf, vertex::BitFlags  >{};
class MyFace    : public Face<   MyUsedTypes, face::FFAdj,  face::Color4b, face::Normal3f, face::VertexRef, face::BitFlags > {};
class MyEdge    : public Edge<   MyUsedTypes, edge::VertexRef, edge::VEAdj, edge::EEAdj, edge::BitFlags> {};

class MyMesh    : public tri::TriMesh< std::vector<MyVertex>, std::vector<MyFace> , std::vector<MyEdge>  > {};


int main( int argc, char **argv )
{
  setvbuf(stdout, NULL, _IONBF, 0);
  if(argc<3)
  {
    printf("Usage alignmeshmesh <meshfilename2.ply> <meshfilename2.ply>\n");
    return -1;
  }
  
  MyMesh fixM,movM;

  tri::io::Importer<MyMesh>::Open(fixM,argv[1]);
  printf( "Mesh0 has %i vert and %i faces\n", fixM.VN(), fixM.FN() );
  tri::io::Importer<MyMesh>::Open(movM,argv[2]);
  printf( "Mesh1 has %i vert and %i faces\n", movM.VN(), movM.FN() );

  math::MarsenneTwisterRNG rnd;
  
  tri::UpdateBounding<MyMesh>::Box(fixM);
  tri::UpdateBounding<MyMesh>::Box(movM);
  
  
  float featureRad = fixM.bbox.Diag() * 0.005;
  
//  RansacFramework<MyMesh,BaseFeatureSet<MeshType> >::EvalNormalVariation(movM,featureRad*5.0);
//  tri::io::ExporterPLY<MyMesh>::Save(movM,"mv0.ply",tri::io::Mask::IOM_VERTQUALITY + tri::io::Mask::IOM_VERTCOLOR);  
//  RansacFramework<MyMesh,BaseFeatureSet>::EvalNormalVariation(movM,featureRad*10.0);
//  tri::io::ExporterPLY<MyMesh>::Save(movM,"mv1.ply",tri::io::Mask::IOM_VERTQUALITY + tri::io::Mask::IOM_VERTCOLOR);  
  
//  RansacFramework<MyMesh,BaseFeatureSet>::EvalNormalVariation(fixM,featureRad*5.0);
//  tri::io::ExporterPLY<MyMesh>::Save(fixM,"fv0.ply",tri::io::Mask::IOM_VERTQUALITY + tri::io::Mask::IOM_VERTCOLOR);  
//  RansacFramework<MyMesh,BaseFeatureSet>::EvalNormalVariation(fixM,featureRad*10.0);
//  tri::io::ExporterPLY<MyMesh>::Save(fixM,"fv1.ply",tri::io::Mask::IOM_VERTQUALITY + tri::io::Mask::IOM_VERTCOLOR);  

  
  Point3f delta = math::GeneratePointInUnitBallUniform<float>(rnd) * fixM.bbox.Diag();
  tri::UpdatePosition<MyMesh>::Translate(movM,delta);
  tri::UpdateNormal<MyMesh>::PerVertexNormalizedPerFaceNormalized(fixM);
  tri::UpdateNormal<MyMesh>::PerVertexNormalizedPerFaceNormalized(movM);
  tri::io::ExporterPLY<MyMesh>::Save(movM,"out.ply");
  int randSeed = clock();
  
  RansacFramework<MyMesh,BaseFeatureSet<MyMesh> > Ran;
  RansacFramework<MyMesh,BaseFeatureSet<MyMesh> >::Param pp;
  pp.samplingRadiusPerc=0.005;
  pp.evalSize=50;
  pp.inlierRatioThr = 0.4;
  Ran.Init(fixM,movM,pp);
   
  std::vector<RansacFramework<MyMesh,BaseFeatureSet<MyMesh> >::Candidate> cVec;
  Ran.Process_SearchEvaluateTriple(cVec,pp);
  
  MyMesh out0; tri::Append<MyMesh,MyMesh>::MeshCopy(out0,movM);   
  tri::UpdatePosition<MyMesh>::Matrix(out0,cVec[0].Tr);
  tri::io::ExporterPLY<MyMesh>::Save(out0,"out0.ply");  
  
  
  MyMesh inlierMesh0;
  Ran.DumpInlier(inlierMesh0,cVec[0],pp);
  tri::io::ExporterPLY<MyMesh>::Save(inlierMesh0,"inlier0.ply");
  
  MyMesh out1; tri::Append<MyMesh,MyMesh>::MeshCopy(out1,movM);
  tri::UpdatePosition<MyMesh>::Matrix(out1,cVec[1].Tr);
  tri::io::ExporterPLY<MyMesh>::Save(out1,"out1.ply");  
  MyMesh inlierMesh1;
  Ran.DumpInlier(inlierMesh1,cVec[1],pp);
  tri::io::ExporterPLY<MyMesh>::Save(inlierMesh1,"inlier1.ply");

  return 0;
}

