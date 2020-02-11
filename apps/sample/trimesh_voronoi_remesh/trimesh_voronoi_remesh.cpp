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
#include<wrap/io_trimesh/import.h>
#include<wrap/io_trimesh/export_ply.h>
#include <vcg/complex/algorithms/voronoi_remesher.h>

using namespace vcg;
using namespace std;

class MyVertex;
class MyFace;
class MyEdge;

struct MyUsedTypes : public vcg::UsedTypes<
        vcg::Use<MyVertex>::AsVertexType,
        vcg::Use<MyFace>::AsFaceType,
        vcg::Use<MyEdge>::AsEdgeType> {};

class MyVertex: public vcg::Vertex<MyUsedTypes,
         vcg::vertex::Coord3f,  vcg::vertex::Normal3f,
         vcg::vertex::Color4b,  vcg::vertex::Qualityd,
         vcg::vertex::VFAdj,    vcg::vertex::VEAdj,
         vcg::vertex::BitFlags, vcg::vertex::Mark> {};
class MyFace   : public vcg::Face<MyUsedTypes,
         vcg::face::VertexRef, vcg::face::Normal3f,
         vcg::face::Color4b, vcg::face::BitFlags,
         vcg::face::VFAdj, vcg::face::FFAdj,
         vcg::face::Mark> {};
class MyEdge: public vcg::Edge<MyUsedTypes,
         vcg::edge::VertexRef, vcg::edge::BitFlags,
         vcg::edge::EEAdj,    vcg::edge::VEAdj> {};

class MyMesh   : public vcg::tri::TriMesh< std::vector<MyVertex>, std::vector<MyFace> > {};

int main( int argc, char **argv )
{
  MyMesh startMesh;
  if(argc < 2 )
  {
    printf("Usage trimesh_voro mesh samplingRadius\n"
           "samplingRadius is in the same unit of the mesh and is approximately the expected edge length");
    return -1;
  }
  printf("Reading %s  \n",argv[1]);
  int ret = tri::io::Importer<MyMesh>::Open(startMesh,argv[1]);
  if(ret!=0)
  {
    printf("Unable to open %s: '%s'\n",argv[1],tri::io::ImporterPLY<MyMesh>::ErrorMsg(ret));
    return -1;
  }
  tri::UpdateBounding<MyMesh>::Box(startMesh);

  float samplingRadius = startMesh.bbox.Diag() * 0.005f;
  if (argc == 3)
  {
	  try {
		samplingRadius = stof(string(argv[2]));
	  } catch (exception &) {}
  }
  std::cout << "Remeshing using sampling radius: " << samplingRadius << std::endl;
  auto remeshed = Remesher<MyMesh>::Remesh(startMesh, samplingRadius, 70.0);
  
  
  tri::io::ExporterPLY<MyMesh>::Save(*remeshed,"Full.ply",tri::io::Mask::IOM_VERTCOLOR|tri::io::Mask::IOM_WEDGTEXCOORD );
  return 0;
}
