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

#include<vcg/simplex/vertex/base.h>
#include<vcg/simplex/vertex/component.h>

#include <vcg/complex/used_types.h>

#include<vcg/simplex/face/base.h>
#include<vcg/simplex/face/component.h>

#include<vcg/simplex/face/topology.h>
#include<vcg/complex/complex.h>

// input output
#include<wrap/io_trimesh/import.h>
#include<wrap/io_trimesh/export.h>

// topology computation
#include<vcg/complex/algorithms/update/topology.h>

// normals
#include<vcg/complex/algorithms/update/normal.h>
#include <vcg/complex/algorithms/intersection.h>

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

  MyMesh m,em;

  if(tri::io::ImporterPLY<MyMesh>::Open(m,argv[1])!=0)
  {
    printf("Error reading file  %s\n",argv[1]);
    exit(0);
  }

  tri::UpdateFlags<MyMesh>::FaceBorderFromFF(m);
  tri::UpdateNormals<MyMesh>::PerVertexNormalized(m);

  printf("Input mesh  vn:%i fn:%i\n",m.vn,m.fn);
  printf( "Mesh has %i vert and %i faces\n", m.vn, m.fn );

  Plane3f slicingPlane;
  Point3f planeCenter = m.bbox.Center();
  slicingPlane.Init(planeCenter,Point3f(0,0,1));

  vcg::IntersectionPlaneMesh<MyMesh, MyMesh, float>(m, slicingPlane, em );
  tri::Clean<MyMesh>::RemoveDuplicateVertex(em);
  std::vector< std::vector<Point3f> > outlines;
  std::vector<Point3f> outline;
  vcg::tri::UpdateFlags<MyMesh>::EdgeClearV(em);
  int nv=0;

  tri::UpdateTopology<MyMesh>::EdgeEdge(em);

  for(size_t i=0;i<em.edge.size();i++)
  {
    if (!em.edge[i].IsV())
    {
      printf("Edge %i (%i %i) ee0=%i ee1=%i\n",i,tri::Index(em,em.edge[i].V(0)),tri::Index(em,em.edge[i].V(1)),
             tri::Index(em,em.edge[i].EEp(0)),tri::Index(em,em.edge[i].EEp(1)));
      MyEdge* startE=&(em.edge[i]);
      int startI = 0;
      int curI = startI;
      MyEdge* curE = startE;
      MyEdge* nextE; int nextI;
      do
      {
        curE->SetV();
        outline.push_back(curE->V(curI)->P());
        nextE=curE->EEp((curI+1)%2);
        nextI=curE->EEi((curI+1)%2);
        curE=nextE;
        curI=nextI;
        nv++;
      }
      while(curE != startE);

      outlines.push_back(outline);
      printf("Found one outline of %i vertices\n\n",outline.size());

      outline.clear();
    }
  }
  printf("Found %i outlines for a total of %i vertices",outlines.size(),nv);

  return 0;
}
