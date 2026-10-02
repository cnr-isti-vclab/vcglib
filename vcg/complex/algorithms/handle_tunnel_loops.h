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
#ifndef VCG_HANDLE_TUNNEL_LOOPS_H
#define VCG_HANDLE_TUNNEL_LOOPS_H

#include <vcg/complex/complex.h>
#include <vcg/complex/append.h>
#include <vcg/complex/algorithms/clean.h>
#include <vcg/complex/algorithms/stat.h>
#include <vcg/complex/algorithms/update/quality.h>
#include <vcg/complex/algorithms/update/topology.h>
#include <vcg/complex/algorithms/reeb_graph.h>
#include <vcg/math/disjoint_set.h>
#include <vcg/math/random_generator.h>
#include <vcg/space/planar_polygon_tessellation.h>
#include <bitset>
#include <map>
#include <numeric>
#include <queue>
#include <set>

/** \file handle_tunnel_loops.h
 * \brief Handle and tunnel loops of a surface: tri::HandleTunnelLoops.
 */

namespace vcg {
namespace tri {

/// \cond
namespace handle_tunnel_detail {
class Vertex; class Face;
struct Types : public UsedTypes<Use<Vertex>::AsVertexType, Use<Face>::AsFaceType> {};
class Vertex : public vcg::Vertex<Types, vertex::Coord3d, vertex::Qualityd, vertex::BitFlags> {};
class Face : public vcg::Face<Types, face::VertexRef, face::FFAdj, face::BitFlags> {};
class Mesh : public TriMesh<std::vector<Vertex>, std::vector<Face> > {};
class GVertex; class GEdge;
struct GTypes : public UsedTypes<Use<GVertex>::AsVertexType, Use<GEdge>::AsEdgeType> {};
class GVertex : public vcg::Vertex<GTypes, vertex::Coord3d, vertex::BitFlags> {};
class GEdge : public vcg::Edge<GTypes, edge::VertexRef, edge::BitFlags> {};
class Graph : public TriMesh<std::vector<GVertex>, std::vector<GEdge> > {};
}
/// \endcond

/** \ingroup trimesh
 * \brief Handle and tunnel loops of an orientable surface.
 *
 * A closed surface of genus g bounds an inside I and an outside O. A \b handle loop bounds
 * a surface in I but not in O: it circles a handle, like the meridian of a torus, which
 * goes around the tube. A \b tunnel loop bounds in O but not in I: it circles a tunnel,
 * like the longitude of a torus, which goes around the hole. The handle loops span a
 * g-dimensional half of the homology of the surface and the tunnel loops the other half.
 * Compute() finds g short handle loops and g short tunnel loops spanning them, as closed
 * vertex-edge paths on the mesh. They are what to cut to remove a handle or to fill a
 * tunnel, e.g. to repair the topology of a scanned or reconstructed surface.
 *
 * \par Usage
 * \code
 * tri::HandleTunnelLoops<MyMesh> ht;
 * ht.Compute(m);
 * for (const auto &loop : ht.handles)      // ht.genus loops
 *   for (int vi : loop) { ... m.vert[vi] ... }
 * MyEdgeMesh em;                           // to save or draw them
 * tri::HandleTunnelLoops<MyMesh>::LoopsToEdgeMesh(m, ht.handles, em);
 * \endcode
 * See also the sample trimesh_handle_tunnel.cpp.
 *
 * \par Requirements
 * - Components: vertex coordinates and face vertex references only. The algorithm works on
 *   an internal copy with its own adjacency, so no adjacency has to be up to date and the
 *   mesh is not modified.
 * - The surface must be connected, orientable and edge-manifold. Small holes are allowed:
 *   each is closed with a fan around its centroid, and the loops are routed along the
 *   hole boundary instead of across it. Larger or knotted holes make the handle/tunnel
 *   classification depend on how they would be closed.
 * - A violated requirement throws vcg::MissingPreconditionException.
 *
 * \par Output
 * #genus, and the loops in #handles and #tunnels, as indices of mesh vertices. A sphere
 * gives no loops.
 *
 * \par Algorithm
 * It follows T. K. Dey, F. Fan and Y. Wang, "An efficient computation of handle and tunnel
 * loops via Reeb graphs", ACM TOG 32(4), 2013, without any volume mesh:
 * 1. the Reeb graph (tri::ReebGraph) of a random height function gives g loops on the
 *    surface, each with a dual loop in the level set just above its lowest point;
 * 2. the local shape at that lowest point tells whether each loop is non-trivial inside or
 *    outside; pushing each loop off the surface on that side gives bases of the homology
 *    of I and of O;
 * 3. inverting the matrix of linking numbers between the loops on the surface and the
 *    pushed ones turns the 2g loops into a handle basis and a tunnel basis;
 * 4. both bases are then tightened: among the loops made of an edge and two shortest paths
 *    to a few base points, the shortest independent ones of each family are kept, and the
 *    base points are moved onto them, until the total length stops decreasing for
 *    Param::patience rounds.
 *
 * Tightening is a heuristic: the loops are short, not the shortest. It dominates the
 * running time, which grows with the size and the genus of the mesh: about 20 s for genus
 * 4 and 640K triangles, 10 s for genus 27 and 92K triangles, 5 minutes for genus 100 and
 * 287K triangles.
 *
 * Before tightening (Param::maxIter = 0) a basis element is a sum of several of the 2g
 * loops and often has more than one component: each becomes its own Loop, so there can be
 * more than g of them. Tightened elements are single loops.
 */
template <class MeshType>
class HandleTunnelLoops
{
public:
  /// A closed vertex-edge path: indices of mesh vertices, consecutive ones joined by an
  /// edge, the last joined to the first.
  typedef std::vector<int> Loop;

  /// Options of Compute(); the defaults suit most meshes.
  struct Param
  {
    unsigned int seed = 0;  ///< seeds the height direction and the projection directions
    int maxIter = 100;      ///< tightening rounds; 0 keeps the loops found on the Reeb graph
    int patience = 10;      ///< stop tightening after this many rounds without a shorter basis
    int attempts = 10;      ///< height directions tried before giving up on degenerate ones
  };

  std::vector<Loop> handles;  ///< the handle loops, around the handles
  std::vector<Loop> tunnels;  ///< the tunnel loops, around the tunnels
  int genus = 0;              ///< genus of the surface, with its holes closed

  /**
   * \brief Compute #genus, #handles and #tunnels of \a m.
   *
   * \param m    the surface; only read
   * \param par  options
   * \throws vcg::MissingPreconditionException if the mesh has non-manifold edges, a hole
   *         whose boundary touches itself, more than one connected component, or no
   *         consistent orientation, or if every height direction tried was degenerate.
   */
  void Compute(MeshType &m, const Param &par = Param())
  {
    handles.clear(); tunnels.clear(); genus = 0;
    w.Clear();
    tri::Append<WMesh, MeshType>::MeshCopy(w, m);
    std::vector<int> toM;
    for (size_t i = 0; i < m.vert.size(); ++i) if (!m.vert[i].IsD()) toM.push_back(int(i));
    tri::UpdateTopology<WMesh>::FaceFace(w);
    tri::MeshAssert<WMesh>::FFTwoManifoldEdge(w);
    CloseHoles();
    if (tri::Clean<WMesh>::CountConnectedComponents(w) != 1)
      throw vcg::MissingPreconditionException("HandleTunnelLoops: the mesh is not connected.");
    bool oriented, orientable;
    tri::Clean<WMesh>::OrientCoherentlyMesh(w, oriented, orientable);
    if (!orientable)
      throw vcg::MissingPreconditionException("HandleTunnelLoops: the mesh is not orientable.");
    if (tri::Stat<WMesh>::ComputeMeshVolume(w) < 0) tri::Clean<WMesh>::FlipMesh(w);
    tri::UpdateTopology<WMesh>::FaceFace(w);

    math::MarsenneTwisterRNG rnd(par.seed);
    for (int attempt = 0; !InitialBases(rnd); ++attempt)
      if (attempt + 1 == par.attempts)
        throw vcg::MissingPreconditionException("HandleTunnelLoops: no generic height direction found.");
    if (genus == 0) return;
    Annotate();

    std::vector<Loop> out[2];
    if (par.maxIter > 0) Tighten(par.maxIter, par.patience, rnd, out);
    else
      for (int t = 0; t < 2; ++t)
        for (const std::vector<int> &c : chain[t])
          for (const Loop &l : ChainLoops(c)) out[t].push_back(l);
    for (int t = 0; t < 2; ++t)
      for (const Loop &l : out[t])
      {
        Loop o = ToInput(l, toM);
        if (!o.empty()) (t == 0 ? handles : tunnels).push_back(o);
      }
  }

  /**
   * \brief Append \a loops to the edge mesh \a em, one closed polyline each.
   *
   * \param m     the mesh the loops are on
   * \param loops loops of \a m, e.g. #handles or #tunnels
   * \param em    an edge mesh with vertex coordinates and edge vertex references
   */
  template <class EdgeMeshType>
  static void LoopsToEdgeMesh(const MeshType &m, const std::vector<Loop> &loops, EdgeMeshType &em)
  {
    for (const Loop &l : loops)
    {
      const size_t base = em.vert.size();
      tri::Allocator<EdgeMeshType>::AddVertices(em, l.size());
      tri::Allocator<EdgeMeshType>::AddEdges(em, l.size());
      for (size_t k = 0; k < l.size(); ++k)
      {
        em.vert[base + k].P() = EdgeMeshType::CoordType::Construct(m.vert[l[k]].cP());
        typename EdgeMeshType::EdgeType &e = em.edge[em.edge.size() - l.size() + k];
        e.V(0) = &em.vert[base + k];
        e.V(1) = &em.vert[base + (k + 1) % l.size()];
      }
    }
  }

private:
  typedef handle_tunnel_detail::Mesh WMesh;
  typedef WMesh::FaceType WFace;
  typedef std::vector<uint64_t> Bits;
  typedef planar_polygon_detail::IndexedPoint2 IPoint;

  WMesh w;                                  // closed, outward oriented copy of the input
  ReebGraph<WMesh> rg;
  handle_tunnel_detail::Graph graph;        // the Reeb graph, an edge per arc going up
  int firstFan;                             // vertices added to close holes start here
  std::vector<Loop> holes;                  // boundary of each closed hole, in order
  std::vector<int> order;                   // vertices by rank
  std::vector<Point3d> faceNrm;             // outward unit normal of each face
  std::vector<double> lift;                 // push-off distance over each face
  std::vector<std::vector<int> > vertEdges, vertFaces;
  std::vector<double> edgeLen;              // edge weights
  std::vector<std::vector<int> > chain[2];  // initial handle and tunnel bases, as edge sets
  std::vector<int> base;                    // lowest point of each Reeb graph loop
  int words;                                // 64-bit words of a homology class

  // ---- Setup ---------------------------------------------------------------------------

  // Close each hole with a fan of triangles around its centroid.
  void CloseHoles()
  {
    firstFan = int(w.vert.size());
    holes.clear();
    std::map<int,int> next;                 // boundary sides, as seen walking the hole
    for (WFace &f : w.face)
      for (int z = 0; z < 3; ++z)
        if (face::IsBorder(f, z))
          if (!next.insert(std::make_pair(int(tri::Index(w, f.V1(z))), int(tri::Index(w, f.V0(z))))).second)
            throw vcg::MissingPreconditionException("HandleTunnelLoops: a hole boundary touches itself.");
    while (!next.empty())
    {
      Loop h;
      for (int v = next.begin()->first; next.count(v); )
      {
        h.push_back(v);
        const int u = next[v];
        next.erase(v);
        v = u;
      }
      holes.push_back(h);
    }
    if (holes.empty()) return;
    tri::Allocator<WMesh>::AddVertices(w, holes.size());
    size_t nf = 0;
    for (size_t k = 0; k < holes.size(); ++k)
    {
      Point3d c(0, 0, 0);
      for (int v : holes[k]) c += w.vert[v].P();
      w.vert[firstFan + k].P() = c / double(holes[k].size());
      nf += holes[k].size();
    }
    tri::Allocator<WMesh>::AddFaces(w, nf);
    size_t fi = w.face.size() - nf;
    for (size_t k = 0; k < holes.size(); ++k)
      for (size_t i = 0; i < holes[k].size(); ++i, ++fi)
      {
        w.face[fi].V(0) = &w.vert[holes[k][i]];
        w.face[fi].V(1) = &w.vert[holes[k][(i + 1) % holes[k].size()]];
        w.face[fi].V(2) = &w.vert[firstFan + k];
      }
    tri::UpdateTopology<WMesh>::FaceFace(w);
  }

  int Lo(int a) const { return int(tri::Index(graph, graph.edge[a].cV(0))); }
  int Hi(int a) const { return int(tri::Index(graph, graph.edge[a].cV(1))); }

  int Other(int e, int v) const { return rg.edgeVert[e].first ^ rg.edgeVert[e].second ^ v; }

  int EdgeOf(int v, int u) const
  {
    for (int e : vertEdges[v]) if (Other(e, v) == u) return e;
    assert(false);
    return -1;
  }

  // Adjacency, weights and normals, once the Reeb graph has numbered the edges.
  void Setup()
  {
    const int nv = int(w.vert.size());
    vertFaces.assign(nv, std::vector<int>());
    faceNrm.clear(); lift.clear();
    for (WFace &f : w.face)
    {
      for (int z = 0; z < 3; ++z) vertFaces[tri::Index(w, f.V(z))].push_back(int(tri::Index(w, f)));
      // A small fraction of the distance from the centroid to the farthest side's line.
      const Point3d n = (f.P(1) - f.P(0)) ^ (f.P(2) - f.P(0));
      double longest = 0;
      for (int z = 0; z < 3; ++z) longest = std::max(longest, Distance(f.P0(z), f.P1(z)));
      faceNrm.push_back(n / std::max(n.Norm(), std::numeric_limits<double>::min()));
      lift.push_back(longest > 0 ? 1e-3 * n.Norm() / (3 * longest) : 0);
    }
    vertEdges.assign(nv, std::vector<int>());
    edgeLen.resize(rg.edgeVert.size());
    double total = 0;
    for (size_t e = 0; e < rg.edgeVert.size(); ++e)
    {
      const int a = rg.edgeVert[e].first, b = rg.edgeVert[e].second;
      vertEdges[a].push_back(int(e));
      vertEdges[b].push_back(int(e));
      edgeLen[e] = Distance(w.vert[a].P(), w.vert[b].P());
      total += edgeLen[e];
    }
    // Loops should not cross the fans that close the holes.
    for (size_t e = 0; e < rg.edgeVert.size(); ++e)
      if (std::max(rg.edgeVert[e].first, rg.edgeVert[e].second) >= firstFan) edgeLen[e] = total;
  }

  // ---- Loops ---------------------------------------------------------------------------

  // Remove the spikes (a, b, a) of a closed walk, also across its ends.
  static void CancelBacktracks(Loop &l)
  {
    Loop s;
    for (int v : l)
    {
      if (!s.empty() && s.back() == v) continue;
      if (s.size() >= 2 && s[s.size() - 2] == v) { s.pop_back(); continue; }
      s.push_back(v);
    }
    size_t b = 0;
    while (s.size() - b >= 3 && (s[b] == s.back() || s[b + 1] == s.back()))
    {
      if (s[b] == s.back()) s.pop_back();
      else { s.pop_back(); ++b; }
    }
    l.assign(s.begin() + b, s.end());
    if (l.size() < 3) l.clear();
  }

  static void Join(Loop &l, const Loop &seg)
  {
    l.insert(l.end(), seg.begin() + (!l.empty() && l.back() == seg.front() ? 1 : 0), seg.end());
  }

  std::vector<int> Edges(const Loop &l) const
  {
    std::vector<int> e;
    for (size_t i = 0; i < l.size(); ++i) e.push_back(EdgeOf(l[i], l[(i + 1) % l.size()]));
    return e;
  }

  // Edges used an odd number of times: the sum of the loops over Z2.
  std::vector<int> Sum(const std::vector<const Loop*> &loops) const
  {
    std::vector<int> all;
    for (const Loop *l : loops) { std::vector<int> e = Edges(*l); all.insert(all.end(), e.begin(), e.end()); }
    std::sort(all.begin(), all.end());
    std::vector<int> odd;
    for (size_t i = 0; i < all.size(); )
    {
      size_t j = i;
      while (j < all.size() && all[j] == all[i]) ++j;
      if ((j - i) % 2) odd.push_back(all[i]);
      i = j;
    }
    return odd;
  }

  // An edge set where every vertex has even degree, as one closed walk per connected part
  // (Hierholzer), leaving out the parts that bound.
  std::vector<Loop> ChainLoops(const std::vector<int> &c) const
  {
    std::map<int, std::vector<int> > adj;
    for (int e : c) { adj[rg.edgeVert[e].first].push_back(e); adj[rg.edgeVert[e].second].push_back(e); }
    std::set<int> left(c.begin(), c.end());
    std::vector<Loop> loops;
    while (!left.empty())
    {
      Loop stack(1, rg.edgeVert[*left.begin()].first), l;
      while (!stack.empty())
      {
        std::vector<int> &a = adj[stack.back()];
        while (!a.empty() && !left.count(a.back())) a.pop_back();
        if (a.empty()) { l.push_back(stack.back()); stack.pop_back(); continue; }
        left.erase(a.back());
        stack.push_back(Other(a.back(), stack.back()));
      }
      l.pop_back();
      CancelBacktracks(l);
      if (!l.empty() && !Zero(Class(Edges(l)))) loops.push_back(l);
    }
    return loops;
  }

  // Map a loop back to the input mesh, going around the holes instead of across them.
  Loop ToInput(const Loop &l, const std::vector<int> &toM) const
  {
    Loop o;
    for (size_t i = 0; i < l.size(); ++i)
    {
      if (l[i] < firstFan) { o.push_back(l[i]); continue; }
      const Loop &h = holes[l[i] - firstFan];
      const int n = int(h.size());
      const int x = int(std::find(h.begin(), h.end(), l[(i + l.size() - 1) % l.size()]) - h.begin());
      const int y = int(std::find(h.begin(), h.end(), l[(i + 1) % l.size()]) - h.begin());
      const int step = (y - x + n) % n <= n / 2 ? 1 : n - 1;
      for (int k = (x + step) % n; k != y; k = (k + step) % n) o.push_back(h[k]);
    }
    CancelBacktracks(o);
    for (int &v : o) v = toM[v];
    return o;
  }

  // ---- Initial bases from the Reeb graph -----------------------------------------------

  // Preimage of a Reeb arc: a path from its lower to its upper node through the triangles
  // the arc crosses, snapped to the vertices on its left.
  Loop ArcPath(int a) const
  {
    const int va = rg.nodeVert[Lo(a)], vb = rg.nodeVert[Hi(a)];
    auto crosses = [&](int e) { const std::vector<int> &p = rg.edgePath[e]; return std::find(p.begin(), p.end(), a) != p.end(); };
    auto reaches = [&](int f) {
      for (int z = 0; z < 3; ++z)
        if (int(tri::Index(w, w.face[f].cV(z))) == vb)
          for (int s : {z, (z + 2) % 3})
          {
            const int e = rg.faceEdge[3 * f + s];
            if (rg.edgeVert[e].second == vb && rg.edgePath[e].back() == a) return true;
          }
      return false;
    };
    std::map<int,int> from;                 // face -> 3 * previous face + side crossed
    std::queue<int> q;
    for (int f : vertFaces[va])
      for (int z = 0; z < 3; ++z)
        if (int(tri::Index(w, w.face[f].cV(z))) == va)
          for (int s : {z, (z + 2) % 3})
          {
            const int e = rg.faceEdge[3 * f + s];
            if (rg.edgeVert[e].first == va && rg.edgePath[e].front() == a && from.insert(std::make_pair(f, -1)).second) q.push(f);
          }
    while (!q.empty() && !reaches(q.front()))
    {
      const int f = q.front(); q.pop();
      for (int z = 0; z < 3; ++z)
      {
        const WFace &F = w.face[f];
        if (face::IsBorder(F, z) || !crosses(rg.faceEdge[3 * f + z])) continue;
        if (from.insert(std::make_pair(int(tri::Index(w, F.cFFp(z))), 3 * f + z)).second) q.push(int(tri::Index(w, F.cFFp(z))));
      }
    }
    if (q.empty()) return Loop();
    Loop l(1, vb);
    for (int c = from[q.front()]; c >= 0; c = from[c / 3])
      if (l.back() != int(tri::Index(w, w.face[c / 3].cV1(c % 3)))) l.push_back(int(tri::Index(w, w.face[c / 3].cV1(c % 3))));
    if (l.back() != va) l.push_back(va);
    std::reverse(l.begin(), l.end());
    return l;
  }

  // The contour just above the lower vertex of edge e0 that crosses e0, traced with the
  // upper side on its left: its upper vertices and its points projected along the height.
  struct Contour { Loop above; std::vector<IPoint> pts; std::set<int> edges; };
  Contour Trace(int e0, const Point3d &du, const Point3d &dv) const
  {
    const int p = rg.edgeVert[e0].first, r = rg.rank[p];
    const double h = 0.5 * (w.vert[p].Q() + w.vert[order[r + 1]].Q());
    auto below = [&](const WFace &f, int z) { return rg.rank[tri::Index(w, f.cV(z))] <= r; };
    int f = -1, z = 0;
    for (int g : vertFaces[p])
      for (int s = 0; s < 3; ++s)
        if (rg.faceEdge[3 * g + s] == e0 && int(tri::Index(w, w.face[g].cV0(s))) == p) { f = g; z = s; }
    Contour c;
    do
    {
      const WFace &F = w.face[f];
      const int e = rg.faceEdge[3 * f + z];
      const Point3d &a = F.cP0(z), &b = F.cP1(z);
      const double ha = F.cV0(z)->cQ(), hb = F.cV1(z)->cQ();
      const Point3d x = a + (b - a) * std::min(1.0, std::max(0.0, (h - ha) / (hb - ha)));
      c.above.push_back(int(tri::Index(w, F.cV1(z))));
      c.pts.push_back(IPoint{Point2d(x * du, x * dv), int(c.pts.size())});
      c.edges.insert(e);
      const WFace &G = *F.cFFp(z);
      const int zg = F.cFFi(z);
      f = int(tri::Index(w, G));
      for (z = 0; z < 3; ++z)
        if (z != zg && below(G, z) && !below(G, (z + 1) % 3)) break;
      assert(z < 3);
    } while (rg.faceEdge[3 * f + z] != e0);
    return c;
  }

  int UpperEdge(int v, int a) const
  {
    for (int e : vertEdges[v]) if (rg.edgeVert[e].first == v && rg.edgePath[e].front() == a) return e;
    return -1;
  }

  static void Frame(const Point3d &d, Point3d &u, Point3d &v)
  {
    u = d ^ (std::abs(d[0]) < 0.6 ? Point3d(1, 0, 0) : Point3d(0, 1, 0));
    u.Normalize();
    v = d ^ u;
  }

  // A loop pushed off the surface, inside (side -1) or outside (+1). Pushing a vertex
  // along any one direction can cross a folded fan of triangles, so the copy never passes
  // over a vertex: it goes around each one on the left, through the centroids of the
  // triangles lifted along their normal and the midpoints of the edges lifted along the
  // bisector of the two triangles. Every segment then stays over a single triangle.
  std::vector<Point3d> PushOff(const Loop &l, double side) const
  {
    std::vector<Point3d> out;
    auto overEdge = [&](int f, int z) {
      const WFace &F = w.face[f];
      const int g = int(tri::Index(w, F.cFFp(z)));
      Point3d b = faceNrm[f] + faceNrm[g];
      b /= std::max(b.Norm(), std::numeric_limits<double>::min());
      return (F.cP0(z) + F.cP1(z)) * 0.5 + b * (side * std::min(lift[f], lift[g]));
    };
    for (size_t i = 0; i < l.size(); ++i)
    {
      const int u = l[(i + l.size() - 1) % l.size()], v = l[i], x = l[(i + 1) % l.size()];
      int f = -1, z = 0;                    // the triangle on the left of u->v
      for (int g : vertFaces[v])
        for (int s = 0; s < 3; ++s)
          if (int(tri::Index(w, w.face[g].cV0(s))) == u && int(tri::Index(w, w.face[g].cV1(s))) == v) { f = g; z = s; }
      out.push_back(overEdge(f, z));
      for (;;)
      {
        const WFace &F = w.face[f];
        out.push_back((F.cP(0) + F.cP(1) + F.cP(2)) / 3.0 + faceNrm[f] * (side * lift[f]));
        const int zn = (z + 1) % 3;         // the side leaving v
        if (int(tri::Index(w, F.cV1(zn))) == x) break;
        out.push_back(overEdge(f, zn));
        z = F.cFFi(zn);
        f = int(tri::Index(w, F.cFFp(zn)));
      }
    }
    return out;
  }

  // A closed polyline projected on a plane, with the depth of its points.
  struct Poly { std::vector<Point2d> p; std::vector<double> z; Box2d box; };
  static Poly Project(const std::vector<Point3d> &pts, const Point3d &du, const Point3d &dv, const Point3d &dw)
  {
    Poly q;
    for (const Point3d &x : pts)
    {
      q.p.push_back(Point2d(x * du, x * dv));
      q.z.push_back(x * dw);
      q.box.Add(q.p.back());
    }
    return q;
  }

  // Linking number over Z2: the parity of the crossings where a passes over b.
  // -1 on a degenerate projection.
  static int Link(const Poly &a, const Poly &b)
  {
    // Inclusive, unlike Box2::Collide: a segment along a projection axis on the border of
    // the overlap can still cross.
    auto overlap = [](const Box2d &x, const Box2d &y) {
      return x.min.X() <= y.max.X() && y.min.X() <= x.max.X() && x.min.Y() <= y.max.Y() && y.min.Y() <= x.max.Y(); };
    if (!overlap(a.box, b.box)) return 0;
    Box2d ov = a.box; ov.Intersect(b.box);
    auto near = [&](const Poly &q) {
      std::vector<size_t> s;
      for (size_t i = 0; i < q.p.size(); ++i)
      {
        Box2d sb; sb.Add(q.p[i]); sb.Add(q.p[(i + 1) % q.p.size()]);
        if (overlap(sb, ov)) s.push_back(i);
      }
      return s;
    };
    const std::vector<size_t> sa = near(a), sb = near(b);
    int count = 0;
    for (size_t i : sa)
      for (size_t j : sb)
      {
        const size_t i1 = (i + 1) % a.p.size(), j1 = (j + 1) % b.p.size();
        using planar_polygon_detail::Orient2D;
        const long double o1 = Orient2D(a.p[i], a.p[i1], b.p[j]), o2 = Orient2D(a.p[i], a.p[i1], b.p[j1]);
        if ((o1 > 0 && o2 > 0) || (o1 < 0 && o2 < 0)) continue;
        const long double o3 = Orient2D(b.p[j], b.p[j1], a.p[i]), o4 = Orient2D(b.p[j], b.p[j1], a.p[i1]);
        if ((o3 > 0 && o4 > 0) || (o3 < 0 && o4 < 0)) continue;
        if (o1 == 0 || o2 == 0 || o3 == 0 || o4 == 0) return -1;
        const double s = double(o3 / (o3 - o4)), t = double(o1 / (o1 - o2));
        const double za = a.z[i] + s * (a.z[i1] - a.z[i]), zb = b.z[j] + t * (b.z[j1] - b.z[j]);
        if (za == zb) return -1;
        if (za > zb) ++count;
      }
    return count % 2;
  }

  // Steps 1-5 of the algorithm for one height direction; false if it is degenerate.
  bool InitialBases(math::MarsenneTwisterRNG &rnd)
  {
    const Point3d d = math::GeneratePointOnUnitSphereUniform<double>(rnd);
    tri::UpdateQuality<WMesh>::VertexFromPlane(w, Plane3d(0, d));
    rg.Compute(w, graph);
    Setup();
    order.assign(w.vert.size(), -1);
    int nv = 0;
    for (size_t i = 0; i < w.vert.size(); ++i)
    {
      order[rg.rank[i]] = int(i);
      if (!vertFaces[i].empty()) ++nv;
    }
    genus = (2 - nv + int(rg.edgeVert.size()) - w.fn) / 2;
    if (genus == 0) return true;

    // Maximum spanning tree of the Reeb graph, each arc weighted by its lower node: every
    // other arc closes a loop whose lowest point is its lower node, a splitting saddle.
    const int nn = int(rg.nodeVert.size()), na = int(graph.edge.size());
    std::vector<int> byLow(na);
    std::iota(byLow.begin(), byLow.end(), 0);
    std::sort(byLow.begin(), byLow.end(), [&](int a, int b) {
      return rg.Below(rg.nodeVert[Lo(b)], rg.nodeVert[Lo(a)]); });
    std::vector<int> node(nn);
    DisjointSet<int> ds;
    for (int &x : node) ds.MakeSet(&x);
    std::vector<std::vector<std::pair<int,int> > > tree(nn);
    std::vector<int> cut;
    for (int a : byLow)
    {
      const int lo = Lo(a), hi = Hi(a);
      if (ds.FindSet(&node[lo]) == ds.FindSet(&node[hi])) { cut.push_back(a); continue; }
      ds.Union(&node[lo], &node[hi]);
      tree[lo].push_back(std::make_pair(hi, a));
      tree[hi].push_back(std::make_pair(lo, a));
    }
    if (int(cut.size()) != genus) return false;

    Point3d du, dv;
    Frame(d, du, dv);
    std::vector<Loop> arcPath(na);
    auto arcLoop = [&](int a) -> const Loop & { if (arcPath[a].empty()) arcPath[a] = ArcPath(a); return arcPath[a]; };
    std::vector<Loop> loops(2 * genus);     // the Reeb graph loops, then their duals
    std::vector<double> side(2 * genus);    // -1 to push a loop inside, +1 outside
    base.clear();
    for (int k = 0; k < genus; ++k)
    {
      // The loop: up the cut arc, then down the tree back to its lower node.
      const int p = Lo(cut[k]);
      std::vector<std::pair<int,int> > from(nn, std::make_pair(-1, -1));
      std::queue<int> q;
      q.push(p);
      from[p].first = p;
      while (!q.empty())
      {
        const int x = q.front(); q.pop();
        for (const std::pair<int,int> &y : tree[x])
          if (from[y.first].first < 0) { from[y.first] = std::make_pair(x, y.second); q.push(y.first); }
      }
      Loop &l = loops[k];
      if (arcLoop(cut[k]).empty()) return false;
      Join(l, arcLoop(cut[k]));
      int last = -1;
      for (int x = Hi(cut[k]); x != p; x = from[x].first)
      {
        last = from[x].second;
        Loop seg = arcLoop(last);
        if (seg.empty()) return false;
        if (Hi(last) == x) std::reverse(seg.begin(), seg.end());
        Join(l, seg);
      }
      l.pop_back();
      CancelBacktracks(l);

      // Its dual: the contour just above the lowest point on the cut arc. Seen along the
      // height, a contour not enclosing the other one at the saddle winds counterclockwise
      // when the region it bounds is inside near it; then the loop is non-trivial inside.
      const int pv = rg.nodeVert[p], e1 = UpperEdge(pv, cut[k]), e2 = UpperEdge(pv, last);
      if (e1 < 0 || e2 < 0 || l.empty()) return false;
      const Contour b1 = Trace(e1, du, dv), b2 = Trace(e2, du, dv);
      if (b1.edges.count(e2)) return false;
      const Point2d o(w.vert[pv].P() * du, w.vert[pv].P() * dv);
      auto far = [&](const Contour &c) {
        size_t best = 0;
        for (size_t i = 0; i < c.pts.size(); ++i)
          if (Distance(c.pts[i].point, o) > Distance(c.pts[best].point, o)) best = i;
        return c.pts[best].point;
      };
      const bool b2in1 = planar_polygon_detail::PointInContour(far(b2), b1.pts, 0);
      const bool inside = planar_polygon_detail::SignedDoubleArea(b2in1 ? b2.pts : b1.pts) > 0;
      loops[genus + k] = b1.above;
      CancelBacktracks(loops[genus + k]);
      if (loops[genus + k].empty()) return false;
      side[k] = inside ? -1 : 1;
      side[genus + k] = -side[k];
      base.push_back(pv);
    }

    // Linking numbers between the loops and their pushed-off copies, and their inverse:
    // row k of the inverse combines the loops into one linked only with pushed copy k.
    const Point3d dw = math::GeneratePointOnUnitSphereUniform<double>(rnd);
    Point3d pu, pv;
    Frame(dw, pu, pv);
    const int n2 = 2 * genus;
    std::vector<Poly> on, off;
    for (int k = 0; k < n2; ++k)
    {
      std::vector<Point3d> pts;
      for (int v : loops[k]) pts.push_back(w.vert[v].P());
      on.push_back(Project(pts, pu, pv, dw));
      off.push_back(Project(PushOff(loops[k], side[k]), pu, pv, dw));
    }
    std::vector<Bits> L(n2, Bits((n2 + 63) / 64, 0));
    for (int i = 0; i < n2; ++i)
      for (int j = 0; j < n2; ++j)
      {
        const int lk = Link(on[i], off[j]);
        if (lk < 0) return false;
        if (lk) Flip(L[i], j);
      }
    if (!Invert(L)) return false;
    // A combination linked with no copy pushed outside bounds inside (Lemma 3.3 of the
    // paper, whose text swaps the two families): linked with one pushed inside, a handle.
    chain[0].clear(); chain[1].clear();
    for (int k = 0; k < n2; ++k)
    {
      std::vector<const Loop*> sel;
      for (int j = 0; j < n2; ++j) if (Bit(L[k], j)) sel.push_back(&loops[j]);
      chain[side[k] < 0 ? 0 : 1].push_back(Sum(sel));
    }
    return int(chain[0].size()) == genus;
  }

  // ---- Z2 linear algebra ---------------------------------------------------------------

  static bool Bit(const Bits &b, int i) { return (b[i >> 6] >> (i & 63)) & 1; }
  static void Flip(Bits &b, int i) { b[i >> 6] ^= uint64_t(1) << (i & 63); }
  static void Xor(Bits &a, const Bits &b) { for (size_t k = 0; k < a.size(); ++k) a[k] ^= b[k]; }
  static bool Zero(const Bits &b) { for (uint64_t x : b) if (x) return false; return true; }
  static bool Odd(const Bits &a, const Bits &b)
  {
    uint64_t x = 0;
    for (size_t k = 0; k < a.size(); ++k) x ^= a[k] & b[k];
    return std::bitset<64>(x).count() % 2 == 1;
  }

  // Invert a square matrix given by rows; false if it is singular.
  static bool Invert(std::vector<Bits> &a)
  {
    const int n = int(a.size());
    std::vector<Bits> inv(n, Bits(a[0].size(), 0));
    for (int i = 0; i < n; ++i) Flip(inv[i], i);
    for (int c = 0; c < n; ++c)
    {
      int r = c;
      while (r < n && !Bit(a[r], c)) ++r;
      if (r == n) return false;
      std::swap(a[r], a[c]); std::swap(inv[r], inv[c]);
      for (int k = 0; k < n; ++k)
        if (k != c && Bit(a[k], c)) { Xor(a[k], a[c]); Xor(inv[k], inv[c]); }
    }
    a.swap(inv);
    return true;
  }

  // Classes inserted so far, reduced: a class is new if it does not reduce to zero.
  struct Span
  {
    std::vector<Bits> rows;
    std::vector<int> pivot;
    bool Insert(Bits x)
    {
      for (size_t i = 0; i < rows.size(); ++i) if (Bit(x, pivot[i])) Xor(x, rows[i]);
      size_t k = 0;
      while (k < x.size() && x[k] == 0) ++k;
      if (k == x.size()) return false;
      int b = 0;
      while (!((x[k] >> b) & 1)) ++b;
      rows.push_back(x);
      pivot.push_back(int(64 * k) + b);
      return true;
    }
  };

  // ---- Tightening ----------------------------------------------------------------------

  std::vector<uint64_t> ann;                // homology class of each edge, words per edge

  Bits Class(const std::vector<int> &edges) const
  {
    Bits x(words, 0);
    for (int e : edges) for (int k = 0; k < words; ++k) x[k] ^= ann[e * words + k];
    return x;
  }

  // Classes of the edges against a cohomology basis from a tree-cotree decomposition: each
  // of the 2g edges outside both trees, plus the dual path joining its faces in the dual
  // tree, crosses exactly one loop of the corresponding homology basis.
  void Annotate()
  {
    const int ne = int(rg.edgeVert.size()), nf = int(w.face.size());
    words = (2 * genus + 63) / 64;
    ann.assign(size_t(ne) * words, 0);
    std::vector<int> parent(nf, -2), parentEdge(nf, -1);
    std::vector<char> dual(ne, 0), primal(ne, 0);
    std::queue<int> q;
    parent[0] = -1; q.push(0);
    while (!q.empty())
    {
      const int f = q.front(); q.pop();
      for (int z = 0; z < 3; ++z)
      {
        const int g = int(tri::Index(w, w.face[f].cFFp(z)));
        if (parent[g] != -2) continue;
        parent[g] = f; parentEdge[g] = rg.faceEdge[3 * f + z]; dual[parentEdge[g]] = 1;
        q.push(g);
      }
    }
    std::vector<char> seen(w.vert.size(), 0);
    std::queue<int> qv;
    seen[rg.edgeVert[0].first] = 1; qv.push(rg.edgeVert[0].first);
    while (!qv.empty())
    {
      const int v = qv.front(); qv.pop();
      for (int e : vertEdges[v])
        if (!dual[e] && !seen[Other(e, v)]) { seen[Other(e, v)] = 1; primal[e] = 1; qv.push(Other(e, v)); }
    }
    std::vector<std::vector<int> > edgeFace(ne);
    for (int f = 0; f < nf; ++f) for (int z = 0; z < 3; ++z) edgeFace[rg.faceEdge[3 * f + z]].push_back(f);
    int k = 0;
    for (int e = 0; e < ne; ++e)
    {
      if (dual[e] || primal[e]) continue;
      ann[size_t(e) * words + k / 64] ^= uint64_t(1) << (k % 64);
      for (int f : edgeFace[e])
        for (; parent[f] >= 0; f = parent[f])
          ann[size_t(parentEdge[f]) * words + k / 64] ^= uint64_t(1) << (k % 64);
      ++k;
    }
    assert(k == 2 * genus);
  }

  // Shortest path tree from r: distances and the edge to the parent.
  void Dijkstra(int r, std::vector<double> &dist, std::vector<int> &parentEdge, std::vector<int> &settled) const
  {
    dist.assign(w.vert.size(), std::numeric_limits<double>::max());
    parentEdge.assign(w.vert.size(), -1);
    settled.clear();
    typedef std::pair<double,int> QE;
    std::priority_queue<QE, std::vector<QE>, std::greater<QE> > q;
    dist[r] = 0; q.push(QE(0, r));
    while (!q.empty())
    {
      const QE t = q.top(); q.pop();
      if (t.first > dist[t.second]) continue;
      settled.push_back(t.second);
      for (int e : vertEdges[t.second])
      {
        const int u = Other(e, t.second);
        if (t.first + edgeLen[e] < dist[u]) { dist[u] = t.first + edgeLen[e]; parentEdge[u] = e; q.push(QE(dist[u], u)); }
      }
    }
  }

  // A basis element: a canonical loop (edge plus the tree paths of its ends to a root)
  // or an initial loop, with its length and class.
  struct Elem { double len; Bits a; int root, edge; std::vector<Loop> loops; };

  bool InFamily(const Bits &x, const std::vector<Bits> &masks) const
  {
    for (const Bits &m : masks) if (Odd(x, m)) return false;
    return true;
  }

  // Shortest canonical loop of every class of each family, over the tree from r.
  void Canonical(int r, const std::vector<Bits> masks[2], std::map<Bits, Elem> best[2]) const
  {
    std::vector<double> dist;
    std::vector<int> par, settled;
    Dijkstra(r, dist, par, settled);
    std::vector<uint64_t> A(w.vert.size() * words, 0);
    for (int v : settled)
      if (v != r)
        for (int k = 0; k < words; ++k)
          A[size_t(v) * words + k] = A[size_t(Other(par[v], v)) * words + k] ^ ann[size_t(par[v]) * words + k];
    Bits x(words);
    for (size_t e = 0; e < rg.edgeVert.size(); ++e)
    {
      const int u = rg.edgeVert[e].first, v = rg.edgeVert[e].second;
      if (par[u] == int(e) || par[v] == int(e) || dist[u] == std::numeric_limits<double>::max()) continue;
      for (int k = 0; k < words; ++k) x[k] = A[size_t(u) * words + k] ^ A[size_t(v) * words + k] ^ ann[e * words + k];
      if (Zero(x)) continue;
      for (int t = 0; t < 2; ++t)
        if (InFamily(x, masks[t]))
        {
          const Elem c{dist[u] + dist[v] + edgeLen[e], x, r, int(e), std::vector<Loop>()};
          typename std::map<Bits, Elem>::iterator it = best[t].find(x);
          if (it == best[t].end()) best[t].insert(std::make_pair(x, c));
          else if (c.len < it->second.len) it->second = c;
        }
    }
  }

  // Both families grow from the same base points: tunnel loops run through the handles
  // and handle loops around the tunnels, so each family seeds the other.
  void Tighten(int maxIter, int patience, math::MarsenneTwisterRNG &rnd, std::vector<Loop> out[2])
  {
    std::vector<Bits> M;
    for (int t = 0; t < 2; ++t) for (const std::vector<int> &c : chain[t]) M.push_back(Class(c));
    const std::vector<Bits> classes = M;
    if (!Invert(M)) throw std::logic_error("HandleTunnelLoops: the initial loops are not a homology basis.");
    std::vector<Bits> masks[2];
    std::vector<Elem> cur[2];
    std::map<Bits, Elem> best[2];
    double total[2];
    for (int t = 0; t < 2; ++t)
    {
      // A class is in the family when its coordinates on the other family vanish.
      for (int j = (1 - t) * genus; j < (2 - t) * genus; ++j)
      {
        Bits col(words, 0);
        for (int i = 0; i < 2 * genus; ++i) if (Bit(M[i], j)) Flip(col, i);
        masks[t].push_back(col);
      }
      for (int i = 0; i < genus; ++i)
      {
        double len = 0;
        for (int e : chain[t][i]) len += edgeLen[e];
        cur[t].push_back(Elem{len, classes[t * genus + i], -1, -1, ChainLoops(chain[t][i])});
      }
      total[t] = std::numeric_limits<double>::max();
    }
    std::set<int> done;
    std::vector<int> roots = base;
    for (int it = 0, stall = 0; it < maxIter && stall < patience; ++it)
    {
      for (int r : roots) if (done.insert(r).second) Canonical(r, masks, best);
      bool shorter = false;
      for (int t = 0; t < 2; ++t)
      {
        std::vector<const Elem*> cand;
        for (const Elem &e : cur[t]) cand.push_back(&e);
        for (typename std::map<Bits, Elem>::const_iterator i = best[t].begin(); i != best[t].end(); ++i) cand.push_back(&i->second);
        std::stable_sort(cand.begin(), cand.end(), [](const Elem *a, const Elem *b) { return a->len < b->len; });
        Span span;
        std::vector<Elem> next;
        for (const Elem *c : cand)
          if (int(next.size()) < genus && span.Insert(c->a)) next.push_back(*c);
        Materialize(next);
        double sum = 0;
        for (const Elem &e : next) sum += e.len;
        cur[t].swap(next);
        if (sum < total[t] * (1 - 1e-9)) shorter = true;
        total[t] = std::min(total[t], sum);
      }
      stall = shorter ? 0 : stall + 1;
      roots.clear();
      for (int t = 0; t < 2; ++t)
        for (const Elem &e : cur[t])
          for (const Loop &l : e.loops)
            for (int k = 0; k < 2; ++k) roots.push_back(l[rnd.generate(unsigned(l.size()))]);
    }
    for (int t = 0; t < 2; ++t)
      for (const Elem &e : cur[t]) out[t].insert(out[t].end(), e.loops.begin(), e.loops.end());
  }

  // Build the paths of the canonical loops, one shortest path tree per root.
  void Materialize(std::vector<Elem> &el) const
  {
    std::map<int, std::vector<Elem*> > byRoot;
    for (Elem &e : el) if (e.loops.empty()) byRoot[e.root].push_back(&e);
    std::vector<double> dist;
    std::vector<int> par, settled;
    for (typename std::map<int, std::vector<Elem*> >::iterator i = byRoot.begin(); i != byRoot.end(); ++i)
    {
      Dijkstra(i->first, dist, par, settled);
      for (Elem *e : i->second)
      {
        Loop up, down;
        for (int x = rg.edgeVert[e->edge].first; ; x = Other(par[x], x)) { up.push_back(x); if (x == i->first) break; }
        for (int x = rg.edgeVert[e->edge].second; x != i->first; x = Other(par[x], x)) down.push_back(x);
        Loop l(up.rbegin(), up.rend());
        l.insert(l.end(), down.begin(), down.end());
        CancelBacktracks(l);
        e->loops.push_back(l);
      }
    }
  }
};

} // end namespace tri
} // end namespace vcg

#endif // VCG_HANDLE_TUNNEL_LOOPS_H
