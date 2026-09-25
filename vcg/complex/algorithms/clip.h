/****************************************************************************
* VCGLib                                                            o o     *
* Visual and Computer Graphics Library                            o     o   *
*                                                                _   O  _   *
* Copyright(C) 2004-2016                                           \/)\/    *
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

#ifndef __VCGLIB_TRI_CLIP
#define __VCGLIB_TRI_CLIP
#include <algorithm>
#include <unordered_map>
#include <vector>
#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/clean.h>
#include <vcg/complex/algorithms/hole.h>
#include <vcg/complex/algorithms/refine.h>
#include <vcg/complex/algorithms/update/bounding.h>
#include <vcg/complex/algorithms/update/normal.h>
#include <vcg/complex/algorithms/update/topology.h>
#include <vcg/space/planar_polygon_tessellation.h>

namespace vcg
{
namespace tri
{


template <class TriMeshType>
int CapPlanarBoundary(TriMeshType &m, const Plane3<typename TriMeshType::ScalarType> &plane);

/** \brief How far from a plane a vertex of \a m may lie and still count as lying on it.

  Rounding, not geometry: 1e-5 of the bounding-box diagonal, or, for a piece small against
  its distance from the origin, 16 units in the last place of its largest coordinate --
  float coordinates are coarser there than any fraction of the piece. It is what
  ClipMeshWithPlane puts on the plane rather than splitting a hair away from, and what
  CapPlanarBoundary takes to lie on it; compare against it to tell whether an edge of a cut
  is still open, with the distance computed by PlaneDistance. Needs a current bounding box.
 */
template <class TriMeshType>
typename TriMeshType::ScalarType PlaneTolerance(const TriMeshType &m)
{
  typedef typename TriMeshType::ScalarType ScalarType;
  ScalarType extent = 0;
  for (int i = 0; i < 3; ++i)
    extent = std::max(extent, std::max(math::Abs(m.bbox.min[i]), math::Abs(m.bbox.max[i])));
  return std::max(ScalarType(1e-5) * m.bbox.Diag(),
                  ScalarType(16) * std::numeric_limits<ScalarType>::epsilon() * extent);
}

/// Signed distance of \a p from \a plane, whose normal must be unit length. Evaluated in
/// double, so that the test against PlaneTolerance adds no rounding of its own.
template <class ScalarType>
ScalarType PlaneDistance(const Plane3<ScalarType> &plane, const Point3<ScalarType> &p)
{
  return ScalarType(Point3d::Construct(plane.Direction()) * Point3d::Construct(p)
                    - double(plane.Offset()));
}

/** \brief Split a mesh along a plane and discard one halfspace.

  Faces crossing the plane are split along it, so what is left ends on a clean planar
  boundary instead of on the nearest pre-existing face borders. The halfspace the plane's
  normal points into is the one kept; a face lying in the plane itself goes with the other
  side, since it bounds nothing on the kept one and a cap would lie on top of it.

  The split is the refine framework driven by the distance to the plane, and it is exact:
  every crossing edge is split where it meets the plane. The one exception is what makes
  it safe when the plane meets vertices the mesh already has: a vertex within
  PlaneTolerance of the plane is put exactly on it and the edges meeting there are left
  whole, so no vertex is duplicated a hair away from one already there. Split there
  instead, those edges leave a cluster of near-coincident vertices whose outline zig-zags
  at rounding level, and cannot be capped.

  The distance rides in a temporary per-vertex attribute rather than in the vertex quality,
  so a mesh carrying a scalar of its own keeps it; the refine framework interpolates
  quality, colour and texture coordinates onto the new vertices in the usual way.

  With \a capCut the boundary the cut opened is filled, and only that boundary: see
  CapPlanarBoundary.

  \param m        the mesh, modified in place.
  \param plane    the cutting plane; its normal need not be unit length.
  \param capCut   fill the boundary the cut opened.
  \param tolerance how close to an endpoint, as a fraction of the edge, a crossing may fall
                  before that endpoint is moved onto the plane instead of the edge being
                  split. Zero splits every crossing edge, which is exact; a vertex within
                  numerical noise of the plane is put on it whatever this says. A larger
                  value leaves fewer thin triangles along the cut, but moves vertices by up
                  to that fraction of an edge, and on long edges -- the cap of an earlier
                  cut, a CAD model -- that is a distance that shows.
  \return true when the mesh was clipped; false when the plane leaves it untouched, either
          because nothing was on the discarded side or because the plane misses it entirely.
          A \a capCut that could not be performed does not make this false.
 */
template <class TriMeshType>
bool ClipMeshWithPlane(
    TriMeshType &m,
    const Plane3<typename TriMeshType::ScalarType> &plane,
    bool capCut = false,
    typename TriMeshType::ScalarType tolerance = 0)
{
  typedef typename TriMeshType::ScalarType ScalarType;
  typedef typename TriMeshType::CoordType CoordType;
  typedef typename TriMeshType::VertexType VertexType;
  typedef typename TriMeshType::VertexPointer VertexPointer;

  if (m.FN() == 0) return false;

  Plane3<ScalarType> pl = plane;
  pl.Normalize();
  const CoordType n = pl.Direction();   // a point p survives when PlaneDistance(pl, p) >= 0

  typename TriMeshType::template PerVertexAttributeHandle<ScalarType> dist =
      Allocator<TriMeshType>::template GetPerVertexAttribute<ScalarType>(m, "vcg::ClipMeshWithPlane::dist");

  // A vertex within numerical noise of the plane -- the tolerance CapPlanarBoundary uses to
  // decide what lies on it -- counts as lying on it, and is put there exactly once it is
  // clear there is something to cut. Split a hair away from it instead, its crossing edges
  // leave a cluster of split points that zig-zag at rounding level, and an outline that does
  // that cannot be tessellated. The move is below anything visible, which is the point: it
  // is the only snap the cut needs.
  UpdateBounding<TriMeshType>::Box(m);
  const ScalarType noise = PlaneTolerance(m);
  bool anyBelow = false, anyAbove = false;
  for (typename TriMeshType::VertexIterator vi = m.vert.begin(); vi != m.vert.end(); ++vi) {
    if ((*vi).IsD()) continue;
    ScalarType s = PlaneDistance(pl, (*vi).cP());
    if (math::Abs(s) <= noise) s = 0;
    dist[&*vi] = s;
    if (s < 0) anyBelow = true; else anyAbove = true;
  }
  if (!anyBelow) {                     // nothing to remove
    Allocator<TriMeshType>::template DeletePerVertexAttribute<ScalarType>(m, dist);
    return false;
  }
  for (typename TriMeshType::VertexIterator vi = m.vert.begin(); vi != m.vert.end(); ++vi)
    if (!(*vi).IsD() && dist[&*vi] == 0) (*vi).P() -= n * PlaneDistance(pl, (*vi).cP());

  if (anyAbove) {
    RequireFFAdjacency(m);
    UpdateTopology<TriMeshType>::FaceFace(m);
    AttributeMidPointFunctor<TriMeshType, ScalarType> midPoint(&m, dist);
    AttributeEdgePredicate<TriMeshType, ScalarType> crossesPlane(dist, tolerance);
    RefineE<TriMeshType, AttributeMidPointFunctor<TriMeshType, ScalarType>,
            AttributeEdgePredicate<TriMeshType, ScalarType> >(m, midPoint, crossesPlane, false);
  }

  // Every edge still crossing the plane is one the predicate declined to split: either one
  // end lies within `tolerance` of the crossing, or the edge is non-manifold and RefineE
  // cannot split it at all. Either way the nearer end stands in for the split vertex, so it
  // goes onto the plane -- for a non-manifold edge however far that is, since the only
  // alternative is an outline that leaves the plane there and cannot be capped. Deciding it
  // here, per edge, matters: judged afterwards, the declined edge has usually been deleted
  // with the discarded side, and the edges left at the vertex can be too short to show that
  // it was ever close enough -- one such vertex breaks the whole outline of the cut, since
  // a single gap is all it takes to stop a loop being a loop.
  if (anyAbove) {
    std::vector<bool> snap(m.vert.size(), false);
    for (typename TriMeshType::FaceIterator fi = m.face.begin(); fi != m.face.end(); ++fi) {
      if ((*fi).IsD()) continue;
      for (int e = 0; e < 3; ++e) {
        VertexPointer a = (*fi).V0(e);
        VertexPointer b = (*fi).V1(e);
        const ScalarType da = dist[a], db = dist[b];
        if (da * db >= 0) continue;
        snap[size_t(Index(m, (da / (da - db) < ScalarType(0.5)) ? a : b))] = true;
      }
    }
    for (size_t vi = 0; vi < m.vert.size(); ++vi) {
      if (!snap[vi] || m.vert[vi].IsD()) continue;
      VertexType &v = m.vert[vi];
      v.P() -= n * dist[&v];
      dist[&v] = 0;
    }
  }

  // After refining and snapping, no face straddles the plane, so each one goes or stays
  // whole. A face lying in the plane itself goes too. The snap flattens a sliver onto the
  // plane now and then -- split points and snapped vertices both sit at exactly zero -- and
  // such a face bounds nothing on the kept side; kept, it would lie under the cap, leaving
  // edges that four faces share.
  for (typename TriMeshType::FaceIterator fi = m.face.begin(); fi != m.face.end(); ++fi) {
    if ((*fi).IsD()) continue;
    const ScalarType d0 = dist[(*fi).V(0)], d1 = dist[(*fi).V(1)], d2 = dist[(*fi).V(2)];
    if (d0 + d1 + d2 < 0 || (d0 == 0 && d1 == 0 && d2 == 0))
      Allocator<TriMeshType>::DeleteFace(m, *fi);
  }
  Clean<TriMeshType>::RemoveUnreferencedVertex(m);
  Allocator<TriMeshType>::CompactEveryVector(m);
  Allocator<TriMeshType>::template DeletePerVertexAttribute<ScalarType>(m, dist);

  if (capCut && m.FN() > 0)
    CapPlanarBoundary(m, pl);
  return true;
}

/** \brief Fill the boundary loops of \a m that lie on \a plane, and only those.

  Used by ClipMeshWithPlane; separate because filling the planar boundary of a mesh that was
  cut some other way is the same job. Holes the mesh already had are left alone.

  All the loops on the plane are tessellated together under the even-odd rule, so
  concentric outlines -- a torus sliced through the plane of its central circle -- cap as a
  ring, and disjoint ones cap separately. When that fails, each outer outline is tried alone
  with the holes inside it, outer and hole told apart by how the mesh winds them: one bad
  loop then does not cost the others their caps, and two shells of the mesh that pass
  through each other, whose outlines cross, get a cap each. What the planar tessellator
  still cannot take is ear-cut on the mesh boundary instead: a single loop it rejects,
  typically one that crosses itself because the surface it came from does, and an outline
  pinched where it touches itself. Ear cutting always closes such a hole, but where the
  outline crosses itself it cannot avoid folding a few triangles over others. It needs
  per-vertex and per-face normals; without them those holes are left open.

  Also left open: an outline that is not closed on the plane (it runs into a hole the mesh
  already had, which filling it would fill as well), an outer outline with holes that the
  tessellator rejects (ear cutting would pave the holes over), and a hole needing ear
  cutting that has a vertex on a non-manifold edge, which the walk around the hole cannot
  cross.

  Returns the number of holes filled.
 */
template <class TriMeshType>
int CapPlanarBoundary(
    TriMeshType &m,
    const Plane3<typename TriMeshType::ScalarType> &plane)
{
  typedef typename TriMeshType::ScalarType ScalarType;
  typedef typename TriMeshType::CoordType CoordType;
  typedef typename TriMeshType::FaceType FaceType;
  typedef typename TriMeshType::VertexPointer VertexPointer;

  if (m.FN() == 0) return 0;

  Plane3<ScalarType> pl = plane;
  pl.Normalize();
  const CoordType n = pl.Direction();

  RequireFFAdjacency(m);
  UpdateTopology<TriMeshType>::FaceFace(m);

  UpdateBounding<TriMeshType>::Box(m);
  const ScalarType tol = PlaneTolerance(m);
  const auto onPlane = [&](const VertexPointer v) {
    return math::Abs(PlaneDistance(pl, v->cP())) <= tol;
  };

  // The loops are found as a plain graph of border edges lying on the plane, not by walking
  // them with Pos::NextB. NextB turns around each vertex until it meets the next border
  // edge, and around a non-manifold fan it can come back to the edge it started from; the
  // usual guard is to refuse any mesh with a non-manifold edge anywhere, which on a scan --
  // where a handful of them is normal -- means never capping at all. Here a defect costs
  // only the part of the outline it touches.
  struct BorderEdge { VertexPointer a, b; FaceType *f; int z; };
  std::vector<BorderEdge> edges;
  for (size_t fi = 0; fi < m.face.size(); ++fi) {
    FaceType &f = m.face[fi];
    if (f.IsD()) continue;
    for (int e = 0; e < 3; ++e) {
      if (!face::IsBorder(f, e)) continue;
      VertexPointer a = f.V0(e), b = f.V1(e);
      if (onPlane(a) && onPlane(b)) edges.push_back(BorderEdge{a, b, &f, e});
    }
  }
  if (edges.size() < 3) return 0;

  std::unordered_map<VertexPointer, std::vector<size_t> > incident;
  for (size_t i = 0; i < edges.size(); ++i) {
    incident[edges[i].a].push_back(i);
    incident[edges[i].b].push_back(i);
  }

  // Connected components of that graph. One whose every vertex has two edges is a simple
  // loop, the planar tessellator's input. One whose every vertex has an even number of them
  // is still a closed outline, pinched where it touches itself; it is ear-cut. Any other is
  // not closed on the plane at all -- it runs into a hole the mesh already had, or breaks
  // where the plane met a non-manifold edge it could not split -- and is left open, since
  // filling it would fill whatever it runs into as well.
  std::vector<int> component(edges.size(), -1);
  std::vector< std::vector<VertexPointer> > loops;     // simple loops, in order
  std::vector<size_t> loopSeed;                        // an edge of each loop
  std::vector<bool> earCut(edges.size(), false);       // by component: ear-cut it
  for (size_t seed = 0; seed < edges.size(); ++seed) {
    if (component[seed] >= 0) continue;
    std::vector<size_t> members(1, seed);
    component[seed] = int(seed);
    bool simple = true, closed = true;
    for (size_t k = 0; k < members.size(); ++k) {
      const VertexPointer ends[2] = { edges[members[k]].a, edges[members[k]].b };
      for (VertexPointer p : ends) {
        const std::vector<size_t> &here = incident[p];
        if (here.size() != 2) simple = false;
        if (here.size() % 2 != 0) closed = false;
        for (size_t other : here)
          if (component[other] < 0) { component[other] = int(seed); members.push_back(other); }
      }
    }
    if (members.size() < 3 || !closed) continue;
    if (!simple) { earCut[seed] = true; continue; }
    // Every vertex has two edges, so following the one we did not arrive by closes the loop.
    std::vector<VertexPointer> loop;
    size_t along = seed;
    VertexPointer at = edges[seed].a;
    do {
      loop.push_back(at);
      at = (edges[along].a == at) ? edges[along].b : edges[along].a;
      const std::vector<size_t> &here = incident[at];
      along = (here[0] == along) ? here[1] : here[0];
    } while (at != edges[seed].a && loop.size() <= members.size());
    loops.push_back(loop);
    loopSeed.push_back(seed);
  }

  // Projected into the plane's own frame: the 3D entry point re-derives the plane and then
  // requires the points to sit on it far more tightly than float coordinates can express.
  CoordType u = (math::Abs(n[0]) < ScalarType(0.9)) ? CoordType(1, 0, 0) : CoordType(0, 1, 0);
  u = (u - n * (n * u)).Normalize();
  const CoordType v = (n ^ u).Normalize();
  std::vector< std::vector<Point2<ScalarType> > > contours(loops.size());
  for (size_t i = 0; i < loops.size(); ++i)
    for (VertexPointer p : loops[i])
      contours[i].push_back(Point2<ScalarType>(u * p->cP(), v * p->cP()));

  std::vector<int> planar;              // triangle corners, indices into planarPoints
  std::vector<VertexPointer> planarPoints;
  int planarLoops = 0;
  const auto tessellate = [&](const std::vector<size_t> &group) {
    std::vector< std::vector<Point2<ScalarType> > > part;
    for (size_t i : group) part.push_back(contours[i]);
    std::vector<int> tri;
    if (!TessellatePlanarContours2(part, tri) || tri.size() < 3) return false;
    const int base = int(planarPoints.size());
    for (size_t i : group) planarPoints.insert(planarPoints.end(), loops[i].begin(), loops[i].end());
    for (int t : tri) planar.push_back(base + t);
    planarLoops += int(group.size());
    return true;
  };

  std::vector<size_t> all(loops.size());
  for (size_t i = 0; i < loops.size(); ++i) all[i] = i;
  // The usual case, and the only one that gets nesting exactly right: every loop at once,
  // under the even-odd rule.
  if (!loops.empty() && !tessellate(all)) {
    // One bad loop -- typically a knot, where a surface that folds over itself makes the
    // true section cross itself -- would otherwise lose the whole cap, and so would two
    // shells of one mesh that pass through each other, whose outlines cross. So each outer
    // outline is tried alone with the holes inside it. Outer and hole are told apart by how
    // the mesh winds them, not by nesting: two crossing outlines each lie partly inside the
    // other, so nesting cannot say which is which, but on a consistently oriented mesh an
    // outer outline always winds against a hole's. The largest outline is outer, and says
    // which way that is.
    const CoordType c = m.bbox.Center();
    std::vector<ScalarType> winding(edges.size(), ScalarType(0));    // by component
    for (size_t e = 0; e < edges.size(); ++e)
      winding[size_t(component[e])] += ((edges[e].b->cP() - c) ^ (edges[e].a->cP() - c)) * n;
    const auto area = [&](size_t i) { return winding[loopSeed[i]]; };
    size_t largest = 0;
    for (size_t i = 1; i < loops.size(); ++i)
      if (math::Abs(area(i)) > math::Abs(area(largest))) largest = i;
    const auto inside = [&](const Point2<ScalarType> &p, const std::vector<Point2<ScalarType> > &cc) {
      bool in = false;
      for (size_t k = 0, l = cc.size() - 1; k < cc.size(); l = k++)
        if (((cc[k][1] > p[1]) != (cc[l][1] > p[1]))
            && (p[0] < (cc[l][0] - cc[k][0]) * (p[1] - cc[k][1]) / (cc[l][1] - cc[k][1]) + cc[k][0]))
          in = !in;
      return in;
    };
    std::vector<size_t> outer, hole;
    for (size_t i = 0; i < loops.size(); ++i)
      ((area(i) > 0) == (area(largest) > 0) ? outer : hole).push_back(i);
    std::vector< std::vector<size_t> > holesOf(loops.size());
    for (size_t h : hole) {
      // A hole belongs to the smallest outer outline around it; one with none around it is
      // an outline the mesh winds the other way, and is capped on its own.
      size_t best = loops.size();
      for (size_t o : outer)
        if (inside(contours[h][0], contours[o])
            && (best == loops.size() || math::Abs(area(o)) < math::Abs(area(best))))
          best = o;
      if (best == loops.size()) outer.push_back(h);
      else holesOf[best].push_back(h);
    }
    // A single outline the tessellator still rejects is ear-cut, which works on the mesh's
    // own boundary and does not care that its projection crosses itself. One with holes is
    // left open: ear cutting would pave the holes over.
    for (size_t o : outer) {
      std::vector<size_t> group(1, o);
      group.insert(group.end(), holesOf[o].begin(), holesOf[o].end());
      if (!tessellate(group) && group.size() == 1) earCut[loopSeed[o]] = true;
    }
  }

  // The ear filler walks the boundary with Pos, which cannot cross a non-manifold edge: it
  // asserts, and without asserts it wanders off into the wrong fan. An outline with a vertex
  // on such an edge is left open instead.
  if (std::find(earCut.begin(), earCut.end(), true) != earCut.end()) {
    std::vector<bool> nonManifold(m.vert.size(), false);
    for (size_t fi = 0; fi < m.face.size(); ++fi) {
      if (m.face[fi].IsD()) continue;
      for (int e = 0; e < 3; ++e)
        if (!face::IsManifold(m.face[fi], e)) {
          nonManifold[size_t(Index(m, m.face[fi].V0(e)))] = true;
          nonManifold[size_t(Index(m, m.face[fi].V1(e)))] = true;
        }
    }
    for (size_t e = 0; e < edges.size(); ++e)
      if (nonManifold[size_t(Index(m, edges[e].a))] || nonManifold[size_t(Index(m, edges[e].b))])
        earCut[size_t(component[e])] = false;
  }
  std::vector<size_t> toEar;                           // the edges of those components
  for (size_t e = 0; e < edges.size(); ++e)
    if (earCut[size_t(component[e])]) toEar.push_back(e);

  // Ear cutting first: it walks the mesh's current boundary, so it runs before any planar
  // cap face is added, and it is handed each hole directly, so it fills those and no other.
  // The ears need face normals, and vertex normals to tell a convex corner from a reflex
  // one. The mesh's own are wrong for that -- along a cut they lie nearly in the plane, and
  // judged by them close to half of the ears come out folded -- so the outline's vertices
  // are given the cap's normal for the duration, and have theirs back afterwards.
  int earHoles = 0;
  if (!toEar.empty() && HasPerVertexNormal(m) && HasPerFaceNormal(m)) {
    // The ear faces take the mesh's own orientation, so the way they will face is read off
    // the outline as the mesh winds it, not assumed from the plane.
    const CoordType c = m.bbox.Center();
    CoordType winding(0, 0, 0);
    for (size_t e : toEar) winding += (edges[e].b->cP() - c) ^ (edges[e].a->cP() - c);
    const CoordType facing = (winding * n > 0) ? n : -n;
    std::unordered_map<VertexPointer, CoordType> saved;
    for (size_t e : toEar) {
      const VertexPointer ends[2] = { edges[e].a, edges[e].b };
      for (VertexPointer p : ends)
        if (saved.emplace(p, p->N()).second) p->N() = facing;
    }
    UpdateNormal<TriMeshType>::PerFaceNormalized(m);

    // Every edge is offered, not one per component: a pinched outline is two holes or one
    // depending on how the faces around the pinch are connected, and an edge whose hole has
    // been filled already is no longer a border, so each hole is filled once.
    typedef typename face::Pos<FaceType> PosType;
    std::vector<PosType> holes;
    for (size_t e : toEar) holes.push_back(PosType(edges[e].f, edges[e].z, edges[e].f->V(edges[e].z)));
    std::vector<FaceType **> kept;                     // adding faces can move the others
    for (PosType &h : holes) kept.push_back(&h.f);
    for (PosType &h : holes) {
      if (!h.IsBorder()) continue;
      Hole<TriMeshType>::template FillHoleEar<MinimumWeightEar<TriMeshType> >(m, h, kept);
      ++earHoles;
    }
    for (const auto &kv : saved) kv.first->N() = kv.second;
    // Around a pinch the filler needs fewer faces than it allocated, and deletes the rest.
    Allocator<TriMeshType>::CompactFaceVector(m);
  }

  if (planar.size() >= 3) {
    const size_t firstNew = m.face.size();
    for (size_t t = 0; t + 2 < planar.size(); t += 3) {
      const int a = planar[t], b = planar[t + 1], c = planar[t + 2];
      if (a == b || b == c || a == c) continue;
      Allocator<TriMeshType>::AddFace(m, planarPoints[size_t(a)], planarPoints[size_t(b)],
                                      planarPoints[size_t(c)]);
    }
    // The tessellator winds by the basis, which need not face the way the cut does: the cap
    // closes the side that was discarded, so it points against the plane normal.
    CoordType facing(0, 0, 0);
    for (size_t fi = firstNew; fi < m.face.size(); ++fi) facing += TriangleNormal(m.face[fi]);
    if (facing * n > 0)
      for (size_t fi = firstNew; fi < m.face.size(); ++fi) std::swap(m.face[fi].V(1), m.face[fi].V(2));
  }
  return planarLoops + earHoles;
}

  template <typename MESH_TYPE>
  class GenericVertexInterpolator
  {
  public:
    typedef typename MESH_TYPE::VertexType VertexType;
    typedef GenericVertexInterpolator<MESH_TYPE> ClassType;
    typedef typename VertexType::CoordType CoordType;
    typedef typename CoordType::ScalarType ScalarType;
    GenericVertexInterpolator(MESH_TYPE &_m) : m(_m) {}
  private:
    MESH_TYPE &m;
  public:
    inline void operator () (const VertexType & v0, const VertexType & v1, const VertexType & v2, const ScalarType & a, const ScalarType & b, VertexType & r) const
    {
      // position
      r.P() = v0.cP() + (v1.cP() - v0.cP()) * a + (v2.cP() - v0.cP()) * b;

      // normal
      if (tri::HasPerVertexNormal(m))
      {
        r.N() = v0.cN() + (v1.cN() - v0.cN()) * a + (v2.cN() - v0.cN()) * b;
      }

      // color
      if (tri::HasPerVertexColor(m))
      {
        vcg::Point4<ScalarType> vc[3];
        vc[0].Import(v0.cC());
        vc[1].Import(v1.cC());
        vc[2].Import(v2.cC());
        const vcg::Point4<ScalarType> rc = (vc[0] + (vc[1] - vc[0]) * a + (vc[2] - vc[0]) * b);
        r.C()[0] = (typename vcg::Color4b::ScalarType)(rc[0]);
        r.C()[1] = (typename vcg::Color4b::ScalarType)(rc[1]);
        r.C()[2] = (typename vcg::Color4b::ScalarType)(rc[2]);
        r.C()[3] = (typename vcg::Color4b::ScalarType)(rc[3]);
      }

      // texcoord
      if (tri::HasPerVertexTexCoord(m))
      {
        const short nt = 1; //typename VertexType::TextureType::N();
        for (short i=0; i<nt; ++i)
        {
          r.T().t(i) = v0.cT().t(i) + (v1.cT().t(i) - v0.cT().t(i)) * a + (v2.cT().t(i) - v0.cT().t(i)) * b;
        }
      }
    }
  };
template <typename TRIMESHTYPE>
class TriMeshClipper
{
public:

  typedef TriMeshClipper<TRIMESHTYPE> ClassType;
  typedef TRIMESHTYPE TriMeshType;
  typedef typename TriMeshType::FaceType FaceType;
  typedef typename FaceType::VertexType VertexType;
  typedef typename VertexType::CoordType CoordType;
  typedef typename CoordType::ScalarType ScalarType;

  /*
    static inline void Box(const Box3<ScalarType> & b, VERTEXINTEPOLATOR & vInterp, TriMeshType & m);

    Clip mesh "m" against an axis aligned box (in place version);

    Notes:
      1) faces marked as deleted are skipped;
      2) faces completely outside box are marked as deleted;
      3) faces completely inside box are left unchanged;
      4) faces intersecting box's sides are marked as deleted:
         they are replaced with proper tesselation; new vertices and faces
         are created, so reallocation could occour; previously saved pointers
         could not to be valid anymore, thus they should be updated;
      5) vInterp functor must implement a n operator with signature
             void operator () (const VERTEX & v0, const VERTEX & v1, const VERTEX & v2, const Scalar & a, const Scalar & b, VERTEX & r);
         its semantic is to intepolate vertex attribute across triangle; a typical implementation is;
           r.P() = v0.P() + a * (v1.P() - v0.P()) + b * (v2.P() - v0.P());  // interpolate position
           r.N() = v0.N() + a * (v1.N() - v0.N()) + b * (v2.N() - v0.N());  // interpolate normal
           ...    // interpolate other vertex attributes
  */
  template <class ScalarType>
      class VertexClipInfo
  {
  public:
//			typedef VertexClipInfo ClassType;

    ScalarType fU;
    ScalarType fV;
    unsigned int idx;
    unsigned int tref;
  };
  typedef typename std::vector< VertexClipInfo<ScalarType> > VertexClipInfoVec;
  class TriangleInfo
  {
  public:
    typedef TriangleInfo ClassType;

    unsigned int v[3];
    unsigned int idx;
  };

  typedef std::vector<TriangleInfo> TriangleInfoVec;

  class EdgeIsect
  {
  public:
    CoordType p;
    unsigned int idx;
  };

  template <typename VERTEXINTEPOLATOR>
  static inline void Box(const Box3<ScalarType> & b, VERTEXINTEPOLATOR & vInterp, TriMeshType & m)
  {
    std::vector<unsigned int> facesToDelete;
    ClassType::Box(b, vInterp, m, facesToDelete);
    for (size_t i=0; i<facesToDelete.size(); ++i)
    {
      m.face[facesToDelete[i]].SetD();
    }
  }

  class EdgeIntersections
  {
  public:
    unsigned int n;
    EdgeIsect isects[6];

    EdgeIntersections(void)
    {
      this->n = 0;
    }
  };

  typedef std::unordered_map<unsigned int, EdgeIntersections> UIntHMap;
  typedef typename UIntHMap::iterator UIntHMap_i;
  typedef typename UIntHMap::value_type UIntHMap_v;

  typedef std::unordered_map<unsigned int, UIntHMap> EdgeMap;
  typedef typename EdgeMap::iterator EdgeMap_i;
  typedef typename EdgeMap::value_type EdgeMap_v;

  typedef typename TriMeshType::FaceIterator FaceIterator;

  template <typename VERTEXINTEPOLATOR, typename FACEINDEXCONTAINER>
  static inline void Box(const Box3<ScalarType> & b, VERTEXINTEPOLATOR & vInterp, TriMeshType & m, FACEINDEXCONTAINER & facesToDelete)
  {
    if (m.fn <= 0)
    {
      return;
    }

    EdgeMap edges;
    VertexClipInfoVec vInfos;
    TriangleInfoVec tInfos;

    CoordType vTriangle[4];
    CoordType vClipped[64];

    CoordType pvP0[64];
    CoordType pvP1[64];

    unsigned int numDeletedTris = 0;
    unsigned int numTriangles = 0;
    unsigned int numVertices = 0;

    unsigned int vIdx = (unsigned int)(m.vn);
    unsigned int tIdx = (unsigned int)(m.fn);

    ScalarType boxOffsets[6];

    boxOffsets[0] =  b.min[0];
    boxOffsets[1] = -b.max[0];
    boxOffsets[2] =  b.min[1];
    boxOffsets[3] = -b.max[1];
    boxOffsets[4] =  b.min[2];
    boxOffsets[5] = -b.max[2];

    UIntHMap emptyMap;
    EdgeIntersections emptyIsects;

    const ScalarType eps = (ScalarType)(1e-6);

    for (FaceIterator it=m.face.begin(); it!=m.face.end(); ++it)
    {
      if ((*it).IsD())
      {
        continue;
      }

      unsigned int cc[3];

      cc[0] = ClassType::BoxClipCode(boxOffsets, (*it).V(0)->P());
      cc[1] = ClassType::BoxClipCode(boxOffsets, (*it).V(1)->P());
      cc[2] = ClassType::BoxClipCode(boxOffsets, (*it).V(2)->P());

      if ((cc[0] | cc[1] | cc[2]) == 0)
      {
        continue;
      }

      const unsigned int refT = (unsigned int)(std::distance(m.face.begin(), it));

      if ((cc[0] & cc[1] & cc[2]) != 0)
      {
        facesToDelete.push_back(refT);
        (*it).SetD();
        numDeletedTris++;
        continue;
      }

      facesToDelete.push_back(refT);

      vTriangle[0] = (*it).V(0)->P();
      vTriangle[1] = (*it).V(1)->P();
      vTriangle[2] = (*it).V(2)->P();
      vTriangle[3] = (*it).V(0)->P();

      unsigned int n, n0, n1;

      ClipPolygonLine(0, b.min[0], vTriangle,  4, pvP1,     n1);
      ClipPolygonLine(1, b.max[0], pvP1,      n1, pvP0,     n0);

      ClipPolygonLine(2, b.min[1], pvP0,      n0, pvP1,     n1);
      ClipPolygonLine(3, b.max[1], pvP1,      n1, pvP0,     n0);

      ClipPolygonLine(4, b.min[2], pvP0,      n0, pvP1,     n1);
      ClipPolygonLine(5, b.max[2], pvP1,      n1, vClipped,  n);

      assert(n < 64);

      unsigned int firstV, lastV;

      if (n > 2)
      {
        if (vClipped[0] == vClipped[n - 1])
        {
          n--;
        }

        const CoordType vU = vTriangle[1] - vTriangle[0];
        const CoordType vV = vTriangle[2] - vTriangle[0];

        const ScalarType tArea = (vU ^ vV).SquaredNorm();
        if (tArea < eps)
        {
          continue;
        }

        unsigned int tvidx[3];
        tvidx[0] = (*it).V(0) - &(*(m.vert.begin()));
        tvidx[1] = (*it).V(1) - &(*(m.vert.begin()));
        tvidx[2] = (*it).V(2) - &(*(m.vert.begin()));

        numTriangles += n - 2;

//				size_t vBegin = vInfos.size();

        VertexClipInfo<ScalarType> vnfo;
        TriangleInfo tnfo;

        unsigned int vmin[3];
        unsigned int vmax[3];

        if (tvidx[0] < tvidx[1])
        {
          vmin[0] = tvidx[0];
          vmax[0] = tvidx[1];
        }
        else
        {
          vmin[0] = tvidx[1];
          vmax[0] = tvidx[0];
        }

        if (tvidx[0] < tvidx[2])
        {
          vmin[1] = tvidx[0];
          vmax[1] = tvidx[2];
        }
        else
        {
          vmin[1] = tvidx[2];
          vmax[1] = tvidx[0];
        }

        if (tvidx[1] < tvidx[2])
        {
          vmin[2] = tvidx[1];
          vmax[2] = tvidx[2];
        }
        else
        {
          vmin[2] = tvidx[2];
          vmax[2] = tvidx[1];
        }

        for (unsigned int i=0; i<n; ++i)
        {
          vnfo.tref = refT;

          const CoordType vP = vClipped[i] - vTriangle[0];

          ScalarType tAreaU = (vU ^ vP).SquaredNorm();
          ScalarType tAreaV = (vP ^ vV).SquaredNorm();

          vnfo.fU = (ScalarType)(sqrt(tAreaU / tArea));
          vnfo.fV = (ScalarType)(sqrt(tAreaV / tArea));

          if (vClipped[i] == vTriangle[0])
          {
            vnfo.idx = tvidx[0];
          }
          else if (vClipped[i] == vTriangle[1])
          {
            vnfo.idx = tvidx[1];
          }
          else if (vClipped[i] == vTriangle[2])
          {
            vnfo.idx = tvidx[2];
          }
          else if (vnfo.fV < eps)
          {
            std::pair<EdgeMap_i, bool> mi = edges.insert(std::make_pair(vmin[1], emptyMap));
            std::pair<UIntHMap_i, bool> hi = (*(mi.first)).second.insert(std::make_pair(vmax[1], emptyIsects));
            bool found = false;
            for (unsigned int s=0; s<(*(hi.first)).second.n; ++s)
            {
              if (vClipped[i] == (*(hi.first)).second.isects[s].p)
              {
                found = true;
                vnfo.idx = (*(hi.first)).second.isects[s].idx;
                break;
              }
            }
            if (!found)
            {
              vnfo.idx = vIdx++;
              numVertices++;
              vInfos.push_back(vnfo);

              (*(hi.first)).second.isects[(*(hi.first)).second.n].p = vClipped[i];
              (*(hi.first)).second.isects[(*(hi.first)).second.n].idx = vnfo.idx;
              (*(hi.first)).second.n++;
            }
          }
          else if (vnfo.fU < eps)
          {
            std::pair<EdgeMap_i, bool> mi = edges.insert(std::make_pair(vmin[0], emptyMap));
            std::pair<UIntHMap_i, bool> hi = (*(mi.first)).second.insert(std::make_pair(vmax[0], emptyIsects));
            bool found = false;
            for (unsigned int s=0; s<(*(hi.first)).second.n; ++s)
            {
              if (vClipped[i] == (*(hi.first)).second.isects[s].p)
              {
                found = true;
                vnfo.idx = (*(hi.first)).second.isects[s].idx;
                break;
              }
            }
            if (!found)
            {
              vnfo.idx = vIdx++;
              numVertices++;
              vInfos.push_back(vnfo);

              (*(hi.first)).second.isects[(*(hi.first)).second.n].p = vClipped[i];
              (*(hi.first)).second.isects[(*(hi.first)).second.n].idx = vnfo.idx;
              (*(hi.first)).second.n++;
            }
          }
          else if ((vnfo.fU + vnfo.fV) >= ((ScalarType)(1.0 - 1e-5)))
          {
            std::pair<EdgeMap_i, bool> mi = edges.insert(std::make_pair(vmin[2], emptyMap));
            std::pair<UIntHMap_i, bool> hi = (*(mi.first)).second.insert(std::make_pair(vmax[2], emptyIsects));
            bool found = false;
            for (unsigned int s=0; s<(*(hi.first)).second.n; ++s)
            {
              if (vClipped[i] == (*(hi.first)).second.isects[s].p)
              {
                found = true;
                vnfo.idx = (*(hi.first)).second.isects[s].idx;
                break;
              }
            }
            if (!found)
            {
              vnfo.idx = vIdx++;
              numVertices++;
              vInfos.push_back(vnfo);

              (*(hi.first)).second.isects[(*(hi.first)).second.n].p = vClipped[i];
              (*(hi.first)).second.isects[(*(hi.first)).second.n].idx = vnfo.idx;
              (*(hi.first)).second.n++;
            }
          }
          else
          {
            vnfo.idx = vIdx++;
            numVertices++;
            vInfos.push_back(vnfo);
          }

          if (i == 0)
          {
            firstV = vnfo.idx;
          }

          if (i > 1)
          {
            tnfo.idx = tIdx++;
            tnfo.v[0] = firstV;
            tnfo.v[1] = lastV;
            tnfo.v[2] = vnfo.idx;

            tInfos.push_back(tnfo);
          }

          lastV = vnfo.idx;
        }
      }
    }

    if (numTriangles == 0)
    {
      return;
    }

    const unsigned int vSize = (unsigned int)(m.vn);
    const unsigned int tSize = (unsigned int)(m.fn);

    typedef Allocator<TriMeshType> TriMeshAllocatorType;

    TriMeshAllocatorType::AddVertices(m, numVertices);
    TriMeshAllocatorType::AddFaces(m, numTriangles);

    unsigned int j = vSize;
    for (size_t i=0; i<vInfos.size(); ++i)
    {
      if (vInfos[i].idx >= vSize)
      {
        const unsigned int tref = vInfos[i].tref;
        vInterp(*(m.face[tref].V(0)), *(m.face[tref].V(1)), *(m.face[tref].V(2)), vInfos[i].fV, vInfos[i].fU, m.vert[j]);
        j++;
      }
    }

    j = tSize;
    for (size_t i=0; i<tInfos.size(); ++i)
    {
      m.face[j].V(0) = &(m.vert[tInfos[i].v[0]]);
      m.face[j].V(1) = &(m.vert[tInfos[i].v[1]]);
      m.face[j].V(2) = &(m.vert[tInfos[i].v[2]]);
      j++;
    }
  }


  /*
    static inline void Box(const Box3<ScalarType> & b, VERTEXINTEPOLATOR & vInterp, const TriMeshType & m, TriMeshType & r);

    Clip mesh "m" against an axis aligned box and put resulting data in mesh "r" (out of place version);

    Notes:
      1) input mesh is not modified;
      2) faces marked as deleted are skipped;
      3) vInterp functor must implement a n operator with signature
             void operator () (const VERTEX & v0, const VERTEX & v1, const VERTEX & v2, const Scalar & a, const Scalar & b, VERTEX & r);
         its semantic is to intepolate vertex attribute across triangle; a typical implementation is;
           r.P() = v0.P() + a * (v1.P() - v0.P()) + b * (v2.P() - v0.P());  // interpolate position
           r.N() = v0.N() + a * (v1.N() - v0.N()) + b * (v2.N() - v0.N());  // interpolate normal
           ...    // interpolate other vertex attributes
  */

  template <typename VERTEXINTEPOLATOR>
  static inline void Box(const Box3<ScalarType> & b, VERTEXINTEPOLATOR & vInterp, const TriMeshType & m, TriMeshType & r)
  {
    r.Clear();

    if (m.fn <= 0)
    {
      return;
    }

    class VertexClipInfo
    {
    public:
      typedef VertexClipInfo ClassType;

      ScalarType fU;
      ScalarType fV;
      unsigned int idx;
      unsigned int tref;
    };

    typedef std::vector<VertexClipInfo> VertexClipInfoVec;

    class TriangleInfo
    {
    public:
      typedef TriangleInfo ClassType;

      unsigned int v[3];
      unsigned int idx;
    };

    typedef std::vector<TriangleInfo> TriangleInfoVec;

    class EdgeIsect
    {
    public:
      CoordType p;
      unsigned int idx;
    };

    class EdgeIntersections
    {
    public:
      unsigned int n;
      EdgeIsect isects[6];

      EdgeIntersections(void)
      {
        this->n = 0;
      }
    };

    typedef std::unordered_map<unsigned int, EdgeIntersections> UIntHMap;
    typedef typename UIntHMap::iterator UIntHMap_i;
    typedef typename UIntHMap::value_type UIntHMap_v;

    typedef std::unordered_map<unsigned int, UIntHMap> EdgeMap;
    typedef typename EdgeMap::iterator EdgeMap_i;
    typedef typename EdgeMap::value_type EdgeMap_v;

    typedef std::unordered_map<unsigned int, unsigned int> UIHMap;
    typedef typename UIHMap::iterator UIHMap_i;

    typedef typename TriMeshType::ConstFaceIterator ConstFaceIterator;

    UIHMap origVertsMap;
    EdgeMap edges;
    VertexClipInfoVec vInfos;
    TriangleInfoVec tInfos;

    CoordType vTriangle[4];
    CoordType vClipped[64];

    CoordType pvP0[64];
    CoordType pvP1[64];

    unsigned int numDeletedTris = 0;
    unsigned int numTriangles = 0;
    unsigned int numVertices = 0;

    unsigned int vIdx = 0;
    unsigned int tIdx = 0;

    ScalarType boxOffsets[6];

    boxOffsets[0] =  b.min[0];
    boxOffsets[1] = -b.max[0];
    boxOffsets[2] =  b.min[1];
    boxOffsets[3] = -b.max[1];
    boxOffsets[4] =  b.min[2];
    boxOffsets[5] = -b.max[2];

    UIntHMap emptyMap;
    EdgeIntersections emptyIsects;

    const ScalarType eps = (ScalarType)(1e-6);

    for (ConstFaceIterator it=m.face.begin(); it!=m.face.end(); ++it)
    {
      if ((*it).IsD())
      {
        continue;
      }

      unsigned int cc[3];

      cc[0] = ClassType::BoxClipCode(boxOffsets, (*it).V(0)->P());
      cc[1] = ClassType::BoxClipCode(boxOffsets, (*it).V(1)->P());
      cc[2] = ClassType::BoxClipCode(boxOffsets, (*it).V(2)->P());

      if ((cc[0] | cc[1] | cc[2]) == 0)
      {
        TriangleInfo tnfo;
        VertexClipInfo vnfo;

        tnfo.idx = tIdx++;

        for (int i=0; i<3; ++i)
        {
          const unsigned int v = (*it).V(i) - &(*(m.vert.begin()));
          std::pair<UIHMap_i, bool> hi = origVertsMap.insert(std::make_pair(v, vIdx));

          if (hi.second)
          {
            vnfo.idx = v;
            vInfos.push_back(vnfo);
            tnfo.v[i] = vIdx++;
          }
          else
          {
            tnfo.v[i] = (*(hi.first)).second;
          }
        }

        tInfos.push_back(tnfo);

        continue;
      }

      if ((cc[0] & cc[1] & cc[2]) != 0)
      {
        numDeletedTris++;
        continue;
      }

      vTriangle[0] = (*it).V(0)->P();
      vTriangle[1] = (*it).V(1)->P();
      vTriangle[2] = (*it).V(2)->P();
      vTriangle[3] = (*it).V(0)->P();

      unsigned int n, n0, n1;

      ClipPolygonLine(0, b.min[0], vTriangle,  4, pvP1,     n1);
      ClipPolygonLine(1, b.max[0], pvP1,      n1, pvP0,     n0);

      ClipPolygonLine(2, b.min[1], pvP0,      n0, pvP1,     n1);
      ClipPolygonLine(3, b.max[1], pvP1,      n1, pvP0,     n0);

      ClipPolygonLine(4, b.min[2], pvP0,      n0, pvP1,     n1);
      ClipPolygonLine(5, b.max[2], pvP1,      n1, vClipped,  n);

      assert(n < 64);

      unsigned int firstV, lastV;


      if (n > 2)
      {
        if (vClipped[0] == vClipped[n - 1])
        {
          n--;
        }

        const CoordType vU = vTriangle[1] - vTriangle[0];
        const CoordType vV = vTriangle[2] - vTriangle[0];

        const ScalarType tArea = (vU ^ vV).SquaredNorm();
        if (tArea < eps)
        {
          continue;
        }

        unsigned int tvidx[3];
        tvidx[0] = (*it).V(0) - &(*(m.vert.begin()));
        tvidx[1] = (*it).V(1) - &(*(m.vert.begin()));
        tvidx[2] = (*it).V(2) - &(*(m.vert.begin()));

        unsigned int refT = (unsigned int)(std::distance(m.face.begin(), it));

        numTriangles += n - 2;

        VertexClipInfo vnfo;
        TriangleInfo tnfo;

        unsigned int vmin[3];
        unsigned int vmax[3];

        if (tvidx[0] < tvidx[1])
        {
          vmin[0] = tvidx[0];
          vmax[0] = tvidx[1];
        }
        else
        {
          vmin[0] = tvidx[1];
          vmax[0] = tvidx[0];
        }

        if (tvidx[0] < tvidx[2])
        {
          vmin[1] = tvidx[0];
          vmax[1] = tvidx[2];
        }
        else
        {
          vmin[1] = tvidx[2];
          vmax[1] = tvidx[0];
        }

        if (tvidx[1] < tvidx[2])
        {
          vmin[2] = tvidx[1];
          vmax[2] = tvidx[2];
        }
        else
        {
          vmin[2] = tvidx[2];
          vmax[2] = tvidx[1];
        }

        for (unsigned int i=0; i<n; ++i)
        {
          vnfo.tref = refT;

          const CoordType vP = vClipped[i] - vTriangle[0];

          ScalarType tAreaU = (vU ^ vP).SquaredNorm();
          ScalarType tAreaV = (vP ^ vV).SquaredNorm();

          vnfo.fU = (ScalarType)(sqrt(tAreaU / tArea));
          vnfo.fV = (ScalarType)(sqrt(tAreaV / tArea));

          unsigned int currVIdx;

          if (vClipped[i] == vTriangle[0])
          {
            std::pair<UIHMap_i, bool> hi = origVertsMap.insert(std::make_pair(tvidx[0], vIdx));
            if (hi.second)
            {
              vnfo.idx = tvidx[0];
              vInfos.push_back(vnfo);
              currVIdx = vIdx++;
            }
            else
            {
              currVIdx = (*(hi.first)).second;
            }
          }
          else if (vClipped[i] == vTriangle[1])
          {
            std::pair<UIHMap_i, bool> hi = origVertsMap.insert(std::make_pair(tvidx[1], vIdx));
            if (hi.second)
            {
              vnfo.idx = tvidx[1];
              vInfos.push_back(vnfo);
              currVIdx = vIdx++;
            }
            else
            {
              currVIdx = (*(hi.first)).second;
            }
          }
          else if (vClipped[i] == vTriangle[2])
          {
            std::pair<UIHMap_i, bool> hi = origVertsMap.insert(std::make_pair(tvidx[2], vIdx));
            if (hi.second)
            {
              vnfo.idx = tvidx[2];
              vInfos.push_back(vnfo);
              currVIdx = vIdx++;
            }
            else
            {
              currVIdx = (*(hi.first)).second;
            }
          }
          else if (vnfo.fV < eps)
          {
            std::pair<EdgeMap_i, bool> mi = edges.insert(std::make_pair(vmin[1], emptyMap));
            std::pair<UIntHMap_i, bool> hi = (*(mi.first)).second.insert(std::make_pair(vmax[1], emptyIsects));
            bool found = false;
            for (unsigned int s=0; s<(*(hi.first)).second.n; ++s)
            {
              if (vClipped[i] == (*(hi.first)).second.isects[s].p)
              {
                found = true;
                vnfo.idx = (unsigned int)(-1);
                currVIdx = (*(hi.first)).second.isects[s].idx;
                break;
              }
            }
            if (!found)
            {
              (*(hi.first)).second.isects[(*(hi.first)).second.n].p = vClipped[i];
              (*(hi.first)).second.isects[(*(hi.first)).second.n].idx = vIdx;
              (*(hi.first)).second.n++;

              vnfo.idx = (unsigned int)(-1);
              numVertices++;
              vInfos.push_back(vnfo);
              currVIdx = vIdx++;
            }
          }
          else if (vnfo.fU < eps)
          {
            std::pair<EdgeMap_i, bool> mi = edges.insert(std::make_pair(vmin[0], emptyMap));
            std::pair<UIntHMap_i, bool> hi = (*(mi.first)).second.insert(std::make_pair(vmax[0], emptyIsects));
            bool found = false;
            for (unsigned int s=0; s<(*(hi.first)).second.n; ++s)
            {
              if (vClipped[i] == (*(hi.first)).second.isects[s].p)
              {
                found = true;
                vnfo.idx = (unsigned int)(-1);
                currVIdx = (*(hi.first)).second.isects[s].idx;
                break;
              }
            }
            if (!found)
            {
              (*(hi.first)).second.isects[(*(hi.first)).second.n].p = vClipped[i];
              (*(hi.first)).second.isects[(*(hi.first)).second.n].idx = vIdx;
              (*(hi.first)).second.n++;

              vnfo.idx = (unsigned int)(-1);
              numVertices++;
              vInfos.push_back(vnfo);
              currVIdx = vIdx++;
            }
          }
          else if ((vnfo.fU + vnfo.fV) >= ((ScalarType)(1.0 - 1e-5)))
          {
            std::pair<EdgeMap_i, bool> mi = edges.insert(std::make_pair(vmin[2], emptyMap));
            std::pair<UIntHMap_i, bool> hi = (*(mi.first)).second.insert(std::make_pair(vmax[2], emptyIsects));
            bool found = false;
            for (unsigned int s=0; s<(*(hi.first)).second.n; ++s)
            {
              if (vClipped[i] == (*(hi.first)).second.isects[s].p)
              {
                found = true;
                vnfo.idx = (unsigned int)(-1);
                currVIdx = (*(hi.first)).second.isects[s].idx;
                break;
              }
            }
            if (!found)
            {
              (*(hi.first)).second.isects[(*(hi.first)).second.n].p = vClipped[i];
              (*(hi.first)).second.isects[(*(hi.first)).second.n].idx = vIdx;
              (*(hi.first)).second.n++;

              vnfo.idx = (unsigned int)(-1);
              numVertices++;
              vInfos.push_back(vnfo);
              currVIdx = vIdx++;
            }
          }
          else
          {
            vnfo.idx = (unsigned int)(-1);
            numVertices++;
            vInfos.push_back(vnfo);
            currVIdx = vIdx++;
          }

          if (i == 0)
          {
            firstV = currVIdx;
          }

          if (i > 1)
          {
            tnfo.idx = tIdx++;
            tnfo.v[0] = firstV;
            tnfo.v[1] = lastV;
            tnfo.v[2] = currVIdx;

            tInfos.push_back(tnfo);
          }

          lastV = currVIdx;
        }
      }
    }

    if (tInfos.empty())
    {
      return;
    }

    typedef Allocator<TriMeshType> TriMeshAllocatorType;

    TriMeshAllocatorType::AddVertices(r, (int)(vInfos.size()));
    TriMeshAllocatorType::AddFaces(r, (int)(tInfos.size()));

    for (size_t i=0; i<vInfos.size(); ++i)
    {
      if (vInfos[i].idx != ((unsigned int)(-1)))
      {
        r.vert[i] = m.vert[vInfos[i].idx];
      }
      else
      {
        const unsigned int tref = vInfos[i].tref;
        vInterp(*(m.face[tref].V(0)), *(m.face[tref].V(1)), *(m.face[tref].V(2)), vInfos[i].fV, vInfos[i].fU, r.vert[i]);
      }
    }

    for (size_t i=0; i<tInfos.size(); ++i)
    {
      r.face[i].V(0) = &(r.vert[tInfos[i].v[0]]);
      r.face[i].V(1) = &(r.vert[tInfos[i].v[1]]);
      r.face[i].V(2) = &(r.vert[tInfos[i].v[2]]);
    }
  }

protected:

  static inline unsigned int BoxClipCode(const ScalarType * offsets, const CoordType & p)
  {
    //const ScalarType eps = (ScalarType)(-1e-5);
    const ScalarType eps = (ScalarType)(0);
    unsigned int code = 0;

    code |= ((( p[0] - offsets[0]) < eps) ? (1 << 0) : (0));
    code |= (((-p[0] - offsets[1]) < eps) ? (1 << 1) : (0));
    code |= ((( p[1] - offsets[2]) < eps) ? (1 << 2) : (0));
    code |= (((-p[1] - offsets[3]) < eps) ? (1 << 3) : (0));
    code |= ((( p[2] - offsets[4]) < eps) ? (1 << 4) : (0));
    code |= (((-p[2] - offsets[5]) < eps) ? (1 << 5) : (0));

    return (code);
  }

  static inline unsigned int InRegion(int mode, const ScalarType & value, const CoordType & p_in)
  {
    //const ScalarType eps = (ScalarType)(-1e-5);
    const ScalarType eps = (ScalarType)(0);
    unsigned int flag = 0;

    switch(mode)
    {
      case 0:
        flag = p_in[0] + eps < value;
        break;
      case 1:
        flag = p_in[0] > value + eps;
        break;
      case 2:
        flag = p_in[1] + eps < value;
        break;
      case 3:
        flag = p_in[1] > value + eps;
        break;
      case 4:
        flag = p_in[2] + eps < value;
        break;
      case 5:
        flag = p_in[2] > value + eps;
        break;
      default:
        break;
    }

    return (flag);
  }

  static inline void CrossPoint(int mode, const ScalarType & value, const CoordType & SP, const CoordType & PP, CoordType & p_out)
  {
    switch(mode)
    {
      case 0:
      case 1:
        p_out[0] = value;
        if ((PP[0] - SP[0]) == ((ScalarType)(0)))
        {
          p_out[1] = PP[1];
          p_out[2] = PP[2];
        }
        else
        {
          p_out[1] = SP[1] + (value - SP[0]) * (PP[1] - SP[1]) / (PP[0] - SP[0]);
          p_out[2] = SP[2] + (value - SP[0]) * (PP[2] - SP[2]) / (PP[0] - SP[0]);
        }
        break;
      case 2:
      case 3:
        p_out[1] = value;
        if ((PP[1] - SP[1]) == ((ScalarType)(0)))
        {
          p_out[0] = PP[0];
          p_out[2] = PP[2];
        }
        else
        {
          p_out[0] = SP[0] + (value - SP[1]) * (PP[0] - SP[0]) / (PP[1] - SP[1]);
          p_out[2] = SP[2] + (value - SP[1]) * (PP[2] - SP[2]) / (PP[1] - SP[1]);
        }
        break;
      case 4:
      case 5:
        p_out[2] = value;
        if ((PP[2] - SP[2]) == ((ScalarType)(0)))
        {
          p_out[0] = PP[0];
          p_out[1] = PP[1];
        }
        else
        {
          p_out[0] = SP[0] + (value - SP[2]) * (PP[0] - SP[0]) / (PP[2] - SP[2]);
          p_out[1] = SP[1] + (value - SP[2]) * (PP[1] - SP[1]) / (PP[2] - SP[2]);
        }
        break;
      default:
        break;
    }
  }

  static inline void ClipPolygonLine(int mode, const ScalarType & value, CoordType * P_in, unsigned int n_in, CoordType * P_out, unsigned int & n_out)
  {
    unsigned int ps;
    CoordType * SP;
    CoordType * PP;

    n_out = 0;
    SP = &P_in[n_in-1];

    if (ClassType::InRegion(mode, value, *SP))
    {
      ps = 0;
    }
    else
    {
      ps = 2;
    }

    for(unsigned int i=0; i<n_in; ++i)
    {
      PP = &(P_in[i]);
      ps = (ps >> 1) | ((ClassType::InRegion(mode, value, *PP)) ? (0) : (2));

      switch(ps)
      {
        case 0:
          break;
        case 1:
          ClassType::CrossPoint(mode, value, *SP, *PP, P_out[n_out]);
          n_out++;
          break;
        case 2:
          ClassType::CrossPoint(mode, value, *SP, *PP, P_out[n_out]);
          n_out++;
          P_out[n_out] = *PP;
          n_out++;
          break;
        case 3:
          P_out[n_out] = *PP;
          n_out++;
          break;
        default:
          break;
      }

      SP = PP;
    }
  }

};

} // end namespace tri
} // end namespace vcg

#endif // __VCGLIB_TRI_CLIP
