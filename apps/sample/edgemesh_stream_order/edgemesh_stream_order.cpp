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

\brief Stream orders (Strahler, Hack) and junction collapse on a tree-shaped edge mesh

The vertices use the optional (OCF) VE adjacency, so a vertex type shared with large
triangle meshes pays for it only while an edge-graph algorithm runs. The program checks
its own results and exits non-zero if any is wrong.
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

static int failures = 0;
#define CHECK(cond) do { if (!(cond)) { ++failures; std::printf("FAILED line %d: %s\n", __LINE__, #cond); } } while (0)

// Build an edge mesh from points and index pairs, in the given order.
static void Build(MyMesh &m, const std::vector<Point3f> &pts, const std::vector<std::pair<int,int>> &edges)
{
  m.Clear();
  m.vert.EnableVEAdjacency();
  tri::Allocator<MyMesh>::AddVertices(m, pts.size());
  for (size_t i = 0; i < pts.size(); ++i) m.vert[i].P() = pts[i];
  tri::Allocator<MyMesh>::AddEdges(m, edges.size());
  for (size_t i = 0; i < edges.size(); ++i)
  {
    m.edge[i].V(0) = &m.vert[edges[i].first];
    m.edge[i].V(1) = &m.vert[edges[i].second];
  }
  tri::UpdateTopology<MyMesh>::VertexEdge(m);
}

//            a   b     d   e
//             \  |      \ /
//              \ |       c
//               \|      /
//                j-----'
//                |
//                r
// j has three children: two tips (Strahler 1) and c, whose two tips make it 2.
// So j is 2 -- one child at the maximum, not two -- and so is the root.
static void TestStrahler()
{
  const std::vector<Point3f> pts = {
    {0,0,0}, {0,1,0}, {-1,2,0}, {0,2,0}, {1,1.5f,0}, {1,2.5f,0}, {2,2,0} };
  enum { r, j, a, b, c, d, e };
  std::vector<std::pair<int,int>> edges = { {r,j}, {j,a}, {j,b}, {j,c}, {c,d}, {c,e} };

  // The answer must not depend on the order the children are stored in; the
  // incremental update this replaces got 3 for j when c was processed after both tips.
  for (int pass = 0; pass < 2; ++pass)
  {
    MyMesh m;
    Build(m, pts, edges);
    const std::vector<int> s = Order::Strahler(m, &m.vert[r]);
    CHECK(s[a] == 1 && s[b] == 1 && s[d] == 1 && s[e] == 1);
    CHECK(s[c] == 2);
    CHECK(s[j] == 2);
    CHECK(s[r] == 2);
    std::reverse(edges.begin(), edges.end());
  }
}

//   up (0,1.5)             side (3,2)
//     |                  /
//     j (0,1) ----------'
//     |
//     r (0,0)
// The straight continuation is short, the side branch is long: the two rules disagree.
static void TestHackRules()
{
  const std::vector<Point3f> pts = { {0,0,0}, {0,1,0}, {0,1.5f,0}, {3,2,0} };
  enum { r, j, up, side };
  MyMesh m;
  Build(m, pts, { {r,j}, {j,up}, {j,side} });

  const std::vector<int> angle = Order::Hack(m, &m.vert[r], Order::MinDeviationAngle);
  CHECK(angle[r] == 1 && angle[j] == 1);
  CHECK(angle[up] == 1 && angle[side] == 2);

  const std::vector<int> longest = Order::Hack(m, &m.vert[r], Order::LongestPath);
  CHECK(longest[up] == 2 && longest[side] == 1);
}

// Rooted at a junction there is no incoming direction; exactly one child may continue
// the stem, where the old code gave order 1 to all of them.
static void TestRootAtJunction()
{
  const std::vector<Point3f> pts = { {0,0,0}, {1,0,0}, {0,3,0}, {-1,0,0} };
  MyMesh m;
  Build(m, pts, { {0,1}, {0,2}, {0,3} });
  const std::vector<int> h = Order::Hack(m, &m.vert[0], Order::MinDeviationAngle);
  CHECK(h[0] == 1);
  CHECK(h[2] == 1);                 // the longest arm continues the stem
  CHECK(h[1] == 2 && h[3] == 2);
}

static void TestLoopIsRejected()
{
  const std::vector<Point3f> pts = { {0,0,0}, {1,0,0}, {1,1,0} };
  MyMesh m;
  Build(m, pts, { {0,1}, {1,2}, {2,0} });
  bool thrown = false;
  try { Order::Strahler(m, &m.vert[0]); }
  catch (const vcg::MissingPreconditionException &) { thrown = true; }
  CHECK(thrown);
}

// Collapse the edge j-c of the first tree onto j: c's two tips move to j, which
// becomes a junction of five edges, and the result is still a valid tree.
static void TestCollapseAtJunction()
{
  const std::vector<Point3f> pts = {
    {0,0,0}, {0,1,0}, {-1,2,0}, {0,2,0}, {1,1.5f,0}, {1,2.5f,0}, {2,2,0} };
  enum { r, j, a, b, c, d, e };
  MyMesh m;
  Build(m, pts, { {r,j}, {j,a}, {j,b}, {j,c}, {c,d}, {c,e} });

  MyEdge *jc = &m.edge[3];
  edge::VEEdgeCollapseToVertex(m, jc, jc->V(0) == &m.vert[j] ? 0 : 1);
  CHECK(m.VN() == 6 && m.EN() == 5);
  std::vector<MyVertex *> star;
  edge::VVStarVE(&m.vert[j], star);
  CHECK(star.size() == 5);

  // Compaction must carry the optional VE pointers along with the vertices.
  tri::Allocator<MyMesh>::CompactEveryVector(m);
  CHECK(m.vert.size() == 6);
  size_t junctions = 0;
  for (auto &v : m.vert)
  {
    edge::VVStarVE(&v, star);
    if (star.size() == 5) ++junctions;
  }
  CHECK(junctions == 1);
  const std::vector<int> s = Order::Strahler(m, &m.vert[r]);
  CHECK(s[r] == 2); // four tips on one junction
}

int main()
{
  TestStrahler();
  TestHackRules();
  TestRootAtJunction();
  TestLoopIsRejected();
  TestCollapseAtJunction();
  if (failures == 0) std::printf("edgemesh_stream_order: all checks passed\n");
  return failures == 0 ? 0 : 1;
}
