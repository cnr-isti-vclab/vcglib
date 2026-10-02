/****************************************************************************
* VCGLib                                                            o o     *
* Visual and Computer Graphics Library                            o     o   *
*                                                                _   O  _   *
* Copyright(C) 2004-2026                                           \/)\/    *
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

/*! \file edgemesh_stream_order.cpp
\ingroup code_sample

\brief Stream orders (Strahler, Hack) of a tree-shaped edge mesh

The vertices use the optional (OCF) VE adjacency: enable it, update it, then root the tree
and compute the orders. The tree is

             d   e
             |  /
       a     c
       |    /
       j---'
       |
       r

where the branch towards a goes straight on and the one towards c is longer, so the two
rules for Hack's main stem disagree at j.
*/

#include <cstdio>
#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/stream_order.h>

using namespace vcg;

class MyVertex;
class MyEdge;
class MyFace;
struct MyUsedTypes : public UsedTypes<Use<MyVertex>::AsVertexType,
                                      Use<MyEdge>::AsEdgeType,
                                      Use<MyFace>::AsFaceType> {};

class MyVertex : public Vertex<MyUsedTypes, vertex::InfoOcf, vertex::Coord3f,
                               vertex::BitFlags, vertex::VEAdjOcf> {};
class MyEdge   : public Edge<MyUsedTypes, edge::VertexRef, edge::VEAdj, edge::BitFlags> {};
class MyFace   : public Face<MyUsedTypes, face::VertexRef> {};
class MyMesh   : public tri::TriMesh<vertex::vector_ocf<MyVertex>, std::vector<MyEdge>,
                                     std::vector<MyFace>> {};

typedef tri::StreamOrder<MyMesh> Order;

int main()
{
  const char *name[] = { "r", "j", "a", "c", "d", "e" };
  const Point3f pts[] = { {0,0,0}, {0,1,0}, {0,1.5f,0}, {1.5f,2,0}, {1.5f,3.5f,0}, {3,2.5f,0} };
  const int edges[][2] = { {0,1}, {1,2}, {1,3}, {3,4}, {3,5} };

  MyMesh m;
  tri::Allocator<MyMesh>::AddVertices(m, 6);
  for (int i = 0; i < 6; ++i) m.vert[i].P() = pts[i];
  tri::Allocator<MyMesh>::AddEdges(m, 5);
  for (int i = 0; i < 5; ++i)
  {
    m.edge[i].V(0) = &m.vert[edges[i][0]];
    m.edge[i].V(1) = &m.vert[edges[i][1]];
  }
  m.vert.EnableVEAdjacency();
  tri::UpdateTopology<MyMesh>::VertexEdge(m);

  MyMesh::VertexPointer root = &m.vert[0];
  const std::vector<int> strahler = Order::Strahler(m, root);
  const std::vector<int> straight = Order::Hack(m, root, Order::MinDeviationAngle);
  const std::vector<int> longest  = Order::Hack(m, root, Order::LongestPath);

  std::printf("vertex  Strahler  Hack (straightest)  Hack (longest)\n");
  for (int i = 0; i < 6; ++i)
    std::printf("  %s        %d             %d                 %d\n", name[i], strahler[i], straight[i], longest[i]);
  return 0;
}
