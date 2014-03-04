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
#include<vcg/complex/algorithms/voronoi_processing.h>


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
class MyEdge    : public Edge< MyUsedTypes>{};
class MyMesh    : public tri::TriMesh< vector<MyVertex>, vector<MyFace> , vector<MyEdge>  > {};



int main( int argc, char **argv )
{
  MyMesh baseMesh, clusteredMesh;
  if(argc < 4 )
  {
    printf("Usage trimesh_voronoiclustering mesh region_num iterNum\n");
     return -1;
  }
  int seed = atoi(argv[2]);
  int iterNum   = atoi(argv[3]);
  printf("Reading %s and sampling %i \n",argv[1],seed);
  int ret= tri::io::ImporterPLY<MyMesh>::Open(baseMesh,argv[1]);
  if(ret!=0)
  {
    printf("Unable to open %s for '%s'\n",argv[1],tri::io::ImporterPLY<MyMesh>::ErrorMsg(ret));
    return -1;
  }

  int randSeed=time(0);
  tri::UpdateTopology<MyMesh>::VertexFace(baseMesh);
  std::vector<MyVertex *> seedVec;
  tri::ClusteringSampler<MyMesh> cs(seedVec);
  tri::SurfaceSampling<MyMesh, vcg::tri::ClusteringSampler<MyMesh> >::SamplingRandomGenerator().initialize(randSeed);
  tri::SurfaceSampling<MyMesh, vcg::tri::ClusteringSampler<MyMesh> >::VertexUniform(baseMesh,cs,seed);
  tri::VoronoiProcessingParameter vpp;
  tri::EuclideanDistance<MyMesh> df;
  tri::VoronoiProcessing<MyMesh>::VoronoiRelaxing(baseMesh, seedVec, iterNum, df, vpp);
  tri::VoronoiProcessing<MyMesh>::TopologicalVertexColoring(baseMesh, seedVec);
  tri::VoronoiProcessing<MyMesh>::ConvertDelaunayTriangulationToMesh(baseMesh,clusteredMesh,seedVec);

  tri::io::ExporterPLY<MyMesh>::Save(baseMesh,"base.ply",tri::io::Mask::IOM_VERTCOLOR );
  tri::io::ExporterPLY<MyMesh>::Save(clusteredMesh,"clustered.ply");
  return 0;
}
