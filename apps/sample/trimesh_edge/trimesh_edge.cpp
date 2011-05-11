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
#include <QtOpenGL/qgl.h>

#include<vcg/simplex/vertex/base.h>
#include<vcg/simplex/vertex/component.h>

#include <vcg/complex/used_types.h>

#include<vcg/simplex/face/base.h>
#include<vcg/simplex/face/component.h>

#include<vcg/simplex/face/topology.h>
#include<vcg/complex/complex.h>
#include<vcg/complex/append.h>

// input output
#include<wrap/io_trimesh/import.h>
#include<wrap/io_trimesh/export.h>

// topology computation
#include<vcg/complex/algorithms/update/topology.h>
#include<vcg/complex/algorithms/update/bounding.h>

// normals
#include<vcg/complex/algorithms/update/normal.h>
#include <vcg/complex/algorithms/intersection.h>
#include <wrap/gl/glu_tessellator_cap.h>

using namespace vcg;
using namespace std;

class MyEdge;
class MyFace;
class MyVertex;
struct MyUsedTypes : public UsedTypes<	Use<MyVertex>   ::AsVertexType,
                                        Use<MyEdge>     ::AsEdgeType,
                                        Use<MyFace>     ::AsFaceType>{};

class MyVertex  : public Vertex<MyUsedTypes, vertex::Coord3f, vertex::Normal3f, vertex::BitFlags  >{};
class MyFace    : public Face< MyUsedTypes, face::FFAdj,  face::VertexRef, face::BitFlags > {};
class MyEdge    : public Edge<MyUsedTypes,edge::EVAdj,edge::BitFlags,edge::EEAdj>{};
class MyMesh    : public tri::TriMesh< vector<MyVertex>, vector<MyFace> , vector<MyEdge>  > {};

int main( int argc, char **argv )
{
  if(argc<2)
  {
    printf("Usage trimesh_base <meshfilename.ply>\n");
    return -1;
  }

  MyMesh m,em,cm,full;

  if(tri::io::ImporterPLY<MyMesh>::Open(m,argv[1])!=0)
  {
    printf("Error reading file  %s\n",argv[1]);
    exit(0);
  }

  tri::UpdateFlags<MyMesh>::FaceBorderFromFF(m);
  tri::UpdateNormals<MyMesh>::PerVertexNormalized(m);
  tri::UpdateBounding<MyMesh>::Box(m);

  printf("Input mesh  vn:%i fn:%i\n",m.vn,m.fn);
  printf( "Mesh has %i vert and %i faces\n", m.vn, m.fn );
  srand(time(0));
  Plane3f slicingPlane;
  Point3f planeCenter = m.bbox.Center();

  for(int i=0;i<10;++i)
  {
    cm.Clear();
    em.Clear();
    Point3f planeDir = Point3f(-0.5f+float(rand())/RAND_MAX,-0.5f+float(rand())/RAND_MAX,-0.5f+float(rand())/RAND_MAX);
    planeDir.Normalize();
    printf("slicing dir %5.2f %5.2f %5.2f\n",planeDir[0],planeDir[1],planeDir[2]);

    slicingPlane.Init(planeCenter+planeDir*0.3f*m.bbox.Diag()*float(rand())/RAND_MAX,planeDir);

    vcg::IntersectionPlaneMesh<MyMesh, MyMesh, float>(m, slicingPlane, em );
    tri::Clean<MyMesh>::RemoveDuplicateVertex(em);
    vcg::tri::CapEdgeMesh(em,cm);

    printf("  edge mesh vn %5i en %5i fn %5i\n",em.vn,em.en,em.fn);
    printf("sliced mesh vn %5i en %5i fn %5i\n",cm.vn,cm.en,cm.fn);

    tri::Append<MyMesh,MyMesh>::Mesh(full,cm);
  }

  tri::io::ExporterPLY<MyMesh>::Save(full,"out.ply",false);

  return 0;
}
