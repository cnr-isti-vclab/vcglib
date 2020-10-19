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
#include<wrap/io_trimesh/export_stl.h>
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
  if(argc < 1 )
  {
    return -1;
  }
  printf("Reading %s  \n",argv[1]);
  int ret = tri::io::Importer<MyMesh>::Open(startMesh,argv[1]);
  if(ret!=0)
  {
    printf("Unable to open %s: '%s'\n",argv[1],tri::io::ImporterPLY<MyMesh>::ErrorMsg(ret));
    return -1;
  }
		
  float srcoeff = 0.005f,bordercrangle=70,internalcrangle=30;

  printf("Remeshing, Enter samplingRadius coefficient: \n");

  cin>>srcoeff;

  printf("Enter Border Crease Angle in degrees: \n");

  cin>>bordercrangle;

  printf("Enter Internal Crease Angle in degrees: \n");

  cin>>internalcrangle;

  int dv=tri::Clean<MyMesh>::RemoveDuplicateVertex(startMesh);
  printf("Removed in startMesh %i duplicated vertices\n",dv);

  int unref = tri::Clean<MyMesh>::RemoveUnreferencedVertex(startMesh);
  printf("Removed %i unreferenced vertices from mesh\n",unref);

  tri::UpdateBounding<MyMesh>::Box(startMesh);

  float samplingRadius = startMesh.bbox.Diag() * srcoeff;

  printf("startMesh.bbox.Diag() = %i\n",startMesh.bbox.Diag());

  std::cout << "Remeshing using sampling radius: " << samplingRadius << std::endl;

  auto remeshed = Remesher<MyMesh>::Remesh(startMesh,samplingRadius,bordercrangle,internalcrangle);

  dv=tri::Clean<MyMesh>::RemoveDuplicateVertex(*remeshed);
  printf("Removed in remeshed %i duplicated vertices\n",dv);

  unref = tri::Clean<MyMesh>::RemoveUnreferencedVertex(*remeshed);
  printf("Removed %i unreferenced vertices from remeshed \n",unref);
  
  tri::io::ExporterPLY<MyMesh>::Save(*remeshed,"Remeshed.ply",tri::io::Mask::IOM_VERTCOLOR|tri::io::Mask::IOM_WEDGTEXCOORD );

  tri::io::ExporterSTL<MyMesh>::Save(*remeshed,"Remeshed.stl",tri::io::Mask::IOM_VERTCOLOR|tri::io::Mask::IOM_WEDGTEXCOORD );

  printf("The end of the application, pls. enter any written symbol(rune) and press Enter to exit.");

  char c;
  cin>>c;
	
  return 0;
}
