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
#ifndef VCG_REEB_GRAPH_H
#define VCG_REEB_GRAPH_H

#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/mesh_assert.h>
#include <vcg/complex/algorithms/update/topology.h>
#include <algorithm>

namespace vcg {
namespace tri {

/**
 * @brief Reeb graph of the vertex quality of a triangle mesh.
 *
 * The Reeb graph collapses every connected component of every level set of a field to
 * a point. Its nodes are the critical vertices (minima, maxima, saddles), its arcs the
 * families of contours swept between them. For a closed orientable surface of genus g it
 * has exactly g independent cycles.
 *
 * The graph is built with the on-line algorithm of Pascucci, Scorzelli, Bremer and
 * Mascarenhas ("Robust on-line computation of Reeb graphs: simplicity and speed", ACM TOG
 * 2007): every edge starts as an arc, every triangle glues the paths of its edges by a
 * merge sort along the field, and a vertex is dropped from the graph as soon as all its
 * triangles are in and it has one arc below and one above.
 *
 * The field is the vertex quality, e.g. a height from UpdateQuality::VertexFromPlane.
 * Vertices are ordered by quality with ties broken by index, so the field does not need to
 * be Morse; it must be finite.
 *
 * The graph is written to an edge mesh: a vertex per node, placed at its mesh vertex, and
 * an edge per arc, from its lower node V(0) to its upper node V(1); VE adjacency is filled
 * when the edge mesh has it. This class keeps the map from the mesh to the graph: the node
 * or the arc of every vertex and the arcs crossed by every edge, as indices of the edge
 * mesh's vertices and edges.
 *
 * Requires per-vertex quality and FF adjacency on an edge-manifold mesh.
 */
template <class MeshType>
class ReebGraph
{
public:
  typedef typename MeshType::FaceType     FaceType;

  std::vector<int> nodeVert;                     ///< mesh vertex of each node
  std::vector<int> vertNode;                     ///< node of each vertex, -1 for regular vertices
  std::vector<int> vertArc;                      ///< arc of each regular vertex, -1 for nodes
  std::vector<int> rank;                         ///< position of each vertex in the total order
  std::vector<int> faceEdge;                     ///< edge of side z of face fi, at 3*fi+z
  std::vector<std::pair<int,int> > edgeVert;     ///< lower and upper vertex of each edge
  std::vector<std::vector<int> > edgePath;       ///< arcs crossed by each edge, going up

  bool Below(int v, int w) const { return rank[v] < rank[w]; }

  template <class EdgeMeshType>
  void Compute(MeshType &m, EdgeMeshType &graph)
  {
    RequirePerVertexQuality(m);
    RequireFFAdjacency(m);
    MeshAssert<MeshType>::FFTwoManifoldEdge(m);
    MeshAssert<MeshType>::VertexQualityFinite(m);

    // Total order of the vertices: by quality, ties by index.
    std::vector<int> order;
    for (size_t i = 0; i < m.vert.size(); ++i)
      if (!m.vert[i].IsD()) order.push_back(int(i));
    std::sort(order.begin(), order.end(), [&](int a, int b) {
      const auto qa = m.vert[a].cQ(), qb = m.vert[b].cQ();
      return qa < qb || (qa == qb && a < b); });
    rank.assign(m.vert.size(), -1);
    for (size_t i = 0; i < order.size(); ++i) rank[order[i]] = int(i);

    // One index per edge, shared by the two faces across it.
    faceEdge.assign(3 * m.face.size(), -1);
    edgeVert.clear();
    std::vector<int> edgeFaces;
    std::vector<int> vertFaces(m.vert.size(), 0);
    for (FaceType &f : m.face) if (!f.IsD())
    {
      const size_t fi = tri::Index(m, f);
      for (int z = 0; z < 3; ++z)
      {
        ++vertFaces[tri::Index(m, f.V(z))];
        if (faceEdge[3 * fi + z] >= 0) continue;
        int v0 = int(tri::Index(m, f.V0(z))), v1 = int(tri::Index(m, f.V1(z)));
        if (Below(v1, v0)) std::swap(v0, v1);
        faceEdge[3 * fi + z] = int(edgeVert.size());
        edgeFaces.push_back(face::IsBorder(f, z) ? 1 : 2);
        if (!face::IsBorder(f, z))
          faceEdge[3 * tri::Index(m, f.FFp(z)) + f.FFi(z)] = int(edgeVert.size());
        edgeVert.push_back(std::make_pair(v0, v1));
      }
    }

    wArc.clear();
    up.assign(m.vert.size(), std::vector<int>());
    down.assign(m.vert.size(), std::vector<int>());
    path.assign(edgeVert.size(), std::vector<int>());
    pathHi.assign(edgeVert.size(), std::vector<int>());
    active.assign(edgeVert.size(), false);
    std::vector<int> label(m.vert.size(), -1);

    // Sweep the triangles by their highest vertex, so that contours are glued close to
    // the level they belong to and vertices leave the graph soon after the sweep passes.
    std::vector<int> faceOrder;
    std::vector<int> faceTop(m.face.size(), -1);
    for (FaceType &f : m.face) if (!f.IsD())
    {
      const int fi = int(tri::Index(m, f));
      for (int z = 0; z < 3; ++z) faceTop[fi] = std::max(faceTop[fi], rank[tri::Index(m, f.V(z))]);
      faceOrder.push_back(fi);
    }
    std::sort(faceOrder.begin(), faceOrder.end(), [&](int a, int b) { return faceTop[a] < faceTop[b]; });

    for (int fi : faceOrder)
    {
      FaceType &f = m.face[fi];
      int v[3], e[3];
      for (int z = 0; z < 3; ++z) { v[z] = int(tri::Index(m, f.V(z))); e[z] = faceEdge[3 * fi + z]; }
      for (int z = 0; z < 3; ++z)
        if (path[e[z]].empty()) NewEdge(e[z]);
      // Sides of the sorted triangle: side z joins V(z) and V(z+1).
      int lo = 0, mid = 1, hi = 2;
      if (Below(v[mid], v[lo])) std::swap(lo, mid);
      if (Below(v[hi], v[mid])) std::swap(mid, hi);
      if (Below(v[mid], v[lo])) std::swap(lo, mid);
      auto side = [&](int a, int b) { return (a + 1) % 3 == b ? e[a] : e[b]; };
      Glue(side(lo, hi), side(lo, mid), side(mid, hi));

      for (int z = 0; z < 3; ++z)
        if (--edgeFaces[e[z]] == 0) Finish(e[z]);
      for (int z = 0; z < 3; ++z)
        if (--vertFaces[v[z]] == 0) RemoveRegular(v[z], label);
    }

    // Arcs merged after a vertex was complete can leave it with one arc on each side.
    for (int vi : order) RemoveRegular(vi, label);

    // The surviving arcs and nodes become the edges and vertices of the graph.
    std::vector<int> arcId(wArc.size(), -1);
    int na = 0;
    for (size_t a = 0; a < wArc.size(); ++a) if (wArc[a].merged < 0) arcId[a] = na++;
    nodeVert.clear();
    vertNode.assign(m.vert.size(), -1);
    vertArc.assign(m.vert.size(), -1);
    for (int vi : order)
    {
      if (!up[vi].empty() || !down[vi].empty()) { vertNode[vi] = int(nodeVert.size()); nodeVert.push_back(vi); }
      else if (label[vi] >= 0) vertArc[vi] = arcId[Resolve(label[vi], 2 * rank[vi])];
    }
    graph.Clear();
    tri::Allocator<EdgeMeshType>::AddVertices(graph, nodeVert.size());
    for (size_t i = 0; i < nodeVert.size(); ++i)
      graph.vert[i].P() = EdgeMeshType::CoordType::Construct(m.vert[nodeVert[i]].cP());
    tri::Allocator<EdgeMeshType>::AddEdges(graph, na);
    for (size_t a = 0; a < wArc.size(); ++a)
      if (arcId[a] >= 0)
      {
        graph.edge[arcId[a]].V(0) = &graph.vert[vertNode[wArc[a].lo]];
        graph.edge[arcId[a]].V(1) = &graph.vert[vertNode[wArc[a].hi]];
      }
    if (tri::HasVEAdjacency(graph)) tri::UpdateTopology<EdgeMeshType>::VertexEdge(graph);

    // Map every edge to the final arcs it crosses.
    edgePath.assign(edgeVert.size(), std::vector<int>());
    auto node = [&](int a, int k) { return int(tri::Index(graph, graph.edge[a].cV(k))); };
    for (size_t ei = 0; ei < edgeVert.size(); ++ei)
    {
      // An edge is monotone: between two points of one arc it can only follow that arc.
      const int u = edgeVert[ei].first, w = edgeVert[ei].second;
      const int au = vertArc[u], aw = vertArc[w];
      const int a = au >= 0 && (au == aw || (aw < 0 && node(au, 1) == vertNode[w])) ? au
                  : aw >= 0 && au < 0 && node(aw, 0) == vertNode[u] ? aw : -1;
      if (a >= 0) { edgePath[ei].push_back(a); continue; }
      const int top = 2 * rank[edgeVert[ei].second];
      int pos = 2 * rank[edgeVert[ei].first] + 1;
      size_t k = 0;
      while (pos < top)
      {
        while (2 * rank[pathHi[ei][k]] < pos) ++k;
        const int y = Resolve(path[ei][k], pos);
        edgePath[ei].push_back(arcId[y]);
        pos = 2 * rank[wArc[y].hi] + 1;
      }
    }
    wArc.clear(); up.clear(); down.clear(); path.clear(); pathHi.clear(); active.clear();
  }

private:
  // Working arc: the active edges whose path contains it, the arc it was merged into,
  // and the arcs split off below it, as (lower node before the split, new lower arc).
  struct WArc
  {
    int lo, hi, merged;
    std::vector<int> edges;
    std::vector<std::pair<int,int> > below;
  };
  std::vector<WArc> wArc;
  std::vector<std::vector<int> > up, down;  // working arcs at each vertex
  std::vector<std::vector<int> > path;      // working path of each edge, frozen when finished
  std::vector<std::vector<int> > pathHi;    // upper node of each arc of a frozen path
  std::vector<bool> active;

  int NewArc(int lo, int hi)
  {
    wArc.push_back(WArc{lo, hi, -1, std::vector<int>(), std::vector<std::pair<int,int> >()});
    const int a = int(wArc.size()) - 1;
    up[lo].push_back(a);
    down[hi].push_back(a);
    return a;
  }

  void NewEdge(int e)
  {
    const int a = NewArc(edgeVert[e].first, edgeVert[e].second);
    path[e].push_back(a);
    wArc[a].edges.push_back(e);
    active[e] = true;
  }

  void Finish(int e)
  {
    active[e] = false;
    for (int a : path[e]) pathHi[e].push_back(wArc[a].hi);
  }

  static void Erase(std::vector<int> &v, int x) { v.erase(std::find(v.begin(), v.end(), x)); }

  // Glue the path of the long side lh of a triangle to the paths of its two short sides.
  void Glue(int lh, int lm, int mh)
  {
    const std::vector<int> p = path[lh];
    std::vector<int> q = path[lm];
    q.insert(q.end(), path[mh].begin(), path[mh].end());
    size_t i = 0, j = 0;
    while (i < p.size() && j < q.size())
    {
      const int a = p[i], b = q[j];
      assert(wArc[a].lo == wArc[b].lo);
      if (a == b) { ++i; ++j; }
      else if (wArc[a].hi == wArc[b].hi)
      {
        // Keep the arc on more edges, to rewrite fewer paths.
        if (wArc[a].edges.size() < wArc[b].edges.size()) Merge(b, a); else Merge(a, b);
        ++i; ++j;
      }
      else if (Below(wArc[a].hi, wArc[b].hi)) { Split(b, a); ++i; }
      else { Split(a, b); ++j; }
    }
  }

  // Arc b runs between the same nodes as a and is the same family of contours.
  void Merge(int a, int b)
  {
    for (int e : wArc[b].edges) if (active[e])
    {
      *std::find(path[e].begin(), path[e].end(), b) = a;
      wArc[a].edges.push_back(e);
    }
    Erase(up[wArc[b].lo], b);
    Erase(down[wArc[b].hi], b);
    wArc[b].merged = a;
    wArc[b].edges.clear();
  }

  // Arc a starts where b starts and ends below it: a becomes the lower part of b.
  void Split(int b, int a)
  {
    WArc &B = wArc[b];
    Erase(up[B.lo], b);
    B.below.push_back(std::make_pair(B.lo, a));
    B.lo = wArc[a].hi;
    up[B.lo].push_back(b);
    std::vector<int> still;
    for (int e : B.edges) if (active[e])
    {
      path[e].insert(std::find(path[e].begin(), path[e].end(), b), a);
      wArc[a].edges.push_back(e);
      still.push_back(e);
    }
    B.edges.swap(still);
  }

  // Once all its triangles are in, a vertex with one arc below and one above is regular:
  // its two arcs become one. Every edge still active crosses it, so its paths go d, u.
  void RemoveRegular(int v, std::vector<int> &label)
  {
    if (up[v].size() != 1 || down[v].size() != 1) return;
    const int d = down[v][0], u = up[v][0];
    for (int e : wArc[u].edges) if (active[e]) Erase(path[e], u);
    wArc[d].hi = wArc[u].hi;
    *std::find(down[wArc[u].hi].begin(), down[wArc[u].hi].end(), u) = d;
    wArc[u].merged = d;
    wArc[u].edges.clear();
    up[v].clear(); down[v].clear();
    label[v] = d;
  }

  // The live arc that now holds position pos (2*rank of a vertex, plus one just above it)
  // of an arc that held it when it was recorded.
  int Resolve(int x, int pos) const
  {
    for (;;)
    {
      const WArc &X = wArc[x];
      if (pos < 2 * rank[X.lo])
      {
        // Lower nodes grow along the history: take the last split below pos.
        x = std::partition_point(X.below.begin(), X.below.end(),
              [&](const std::pair<int,int> &s) { return 2 * rank[s.first] < pos; })[-1].second;
      }
      else if (X.merged >= 0) x = X.merged;
      else return x;
    }
  }
};

} // end namespace tri
} // end namespace vcg

#endif // VCG_REEB_GRAPH_H
