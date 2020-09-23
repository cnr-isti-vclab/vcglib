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

#include<wrap/io_trimesh/import_off.h>
#include<wrap/io_trimesh/export_off.h>

#include<vcg/complex/algorithms/point_sampling.h>
#include<vcg/complex/algorithms/create/platonic.h>

using namespace vcg;
using namespace std;

class MyEdge;
class MyFace;
class MyVertex;
struct MyUsedTypes : public UsedTypes<	Use<MyVertex>   ::AsVertexType,
                                        Use<MyEdge>     ::AsEdgeType,
                                        Use<MyFace>     ::AsFaceType>{};

class MyVertex  : public Vertex<MyUsedTypes, vertex::Coord3f, vertex::Normal3f, vertex::BitFlags  >{};
class MyFace    : public Face< MyUsedTypes, face::FFAdj,  face::Normal3f, face::VertexRef, face::BitFlags > {};
class MyEdge    : public Edge<MyUsedTypes>{};
class MyMesh    : public tri::TriMesh< vector<MyVertex>, vector<MyFace> , vector<MyEdge>  > {};

int main( int argc, char **argv )
{
  if(argc<2)
  {
    printf("Usage trimesh_base <meshfilename.obj> radius\n");
    return -1;
  }

  MyMesh m;

  if(tri::io::ImporterOFF<MyMesh>::Open(m,argv[1])!=0)
  {
    printf("Error reading file  %s\n",argv[1]);
    exit(0);
  }
  tri::SurfaceSampling<MyMesh,tri::TrivialSampler<MyMesh> >::SamplingRandomGenerator().initialize(time(0));

  //----------------------------------------------------------------------
  // Advanced Sample
  // Make a feature dependent Poisson Disk sampling
  MyMesh MontecarloSurfaceMesh;

  std::vector<Point3f> sampleVec;
  tri::TrivialSampler<MyMesh> mps(sampleVec);
  int t0=clock();
  tri::SurfaceSampling<MyMesh,tri::TrivialSampler<MyMesh> >::Montecarlo(m,mps,10000000);
  int t1=clock();

  printf("Computed a montecarlo distribution of %lu vertices radius is %6.3f\n",sampleVec.size(),float(t1-t0)/float(CLOCKS_PER_SEC));

  return 0;
}
