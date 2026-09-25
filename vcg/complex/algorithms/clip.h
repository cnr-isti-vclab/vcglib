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
#include <unordered_map>
#include <vector>
#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/clean.h>
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

/** \brief Split a mesh along a plane and discard one halfspace.

  Faces crossing the plane are split along it, so what is left ends on a clean planar
  boundary instead of on the nearest pre-existing face borders. The halfspace the plane's
  normal points into is the one kept.

  The split is the refine framework driven by the distance to the plane, which is what
  makes it safe when the plane meets vertices the mesh already has: the edge predicate
  declines to split an edge whose crossing falls on, or within `tolerance` of, an endpoint,
  so no vertex is ever duplicated a hair away from one already there. Splitting every
  crossing edge unconditionally instead leaves pairs of near-coincident vertices along the
  cut, and a boundary that runs through both of them is not a simple loop -- it cannot be
  capped, and it confuses anything that walks it.

  The distance rides in a temporary per-vertex attribute rather than in the vertex quality,
  so a mesh carrying a scalar of its own keeps it; the refine framework interpolates
  quality, colour and texture coordinates onto the new vertices in the usual way.

  With \a capCut the boundary the cut opened is filled. All the loops it produced are
  tessellated together under the even-odd rule, so a cut that leaves concentric outlines --
  a torus sliced through the plane of its central circle -- is capped as a ring rather than
  as overlapping discs, and one that leaves disjoint outlines gets one cap each. Holes the
  mesh already had are left alone. Capping needs an edge-manifold mesh and a cut whose
  outline is a set of simple loops; when it is not, the geometry is still clipped and the
  cut is left open.

  \param m        the mesh, modified in place.
  \param plane    the cutting plane; its normal need not be unit length.
  \param capCut   fill the boundary the cut opened.
  \param tolerance how close to an endpoint a crossing may fall before the edge is left
                  unsplit, as a fraction of the edge. The refine framework's own default.
  \return true when the mesh was clipped; false when the plane leaves it untouched, either
          because nothing was on the discarded side or because the plane misses it entirely.
          A \a capCut that could not be performed does not make this false.
 */
template <class TriMeshType>
bool ClipMeshWithPlane(
    TriMeshType &m,
    const Plane3<typename TriMeshType::ScalarType> &plane,
    bool capCut = false,
    typename TriMeshType::ScalarType tolerance = 0.02)
{
  typedef typename TriMeshType::ScalarType ScalarType;
  typedef typename TriMeshType::CoordType CoordType;
  typedef typename TriMeshType::VertexType VertexType;
  typedef typename TriMeshType::FaceType FaceType;
  typedef typename TriMeshType::VertexPointer VertexPointer;

  if (m.FN() == 0) return false;

  Plane3<ScalarType> pl = plane;
  pl.Normalize();
  const CoordType n = pl.Direction();
  const ScalarType d = -pl.Offset();   // a point p survives when n*p + d >= 0

  typename TriMeshType::template PerVertexAttributeHandle<ScalarType> dist =
      Allocator<TriMeshType>::template GetPerVertexAttribute<ScalarType>(m, "vcg::ClipMeshWithPlane::dist");

  bool anyBelow = false, anyAbove = false;
  for (typename TriMeshType::VertexIterator vi = m.vert.begin(); vi != m.vert.end(); ++vi) {
    if ((*vi).IsD()) continue;
    const ScalarType s = n * (*vi).cP() + d;
    dist[&*vi] = s;
    if (s < 0) anyBelow = true; else anyAbove = true;
  }
  if (!anyBelow) {                     // nothing to remove
    Allocator<TriMeshType>::template DeletePerVertexAttribute<ScalarType>(m, dist);
    return false;
  }

  if (anyAbove) {
    RequireFFAdjacency(m);
    UpdateTopology<TriMeshType>::FaceFace(m);
    AttributeMidPointFunctor<TriMeshType, ScalarType> midPoint(&m, dist);
    AttributeEdgePredicate<TriMeshType, ScalarType> crossesPlane(dist, tolerance);
    RefineE<TriMeshType, AttributeMidPointFunctor<TriMeshType, ScalarType>,
            AttributeEdgePredicate<TriMeshType, ScalarType> >(m, midPoint, crossesPlane, false);
  }

  // After refining, no face straddles the plane, so each one goes or stays whole.
  for (typename TriMeshType::FaceIterator fi = m.face.begin(); fi != m.face.end(); ++fi) {
    if ((*fi).IsD()) continue;
    const ScalarType mean =
        (dist[(*fi).V(0)] + dist[(*fi).V(1)] + dist[(*fi).V(2)]) / ScalarType(3);
    if (mean < 0) Allocator<TriMeshType>::DeleteFace(m, *fi);
  }
  Clean<TriMeshType>::RemoveUnreferencedVertex(m);
  Allocator<TriMeshType>::CompactEveryVector(m);
  Allocator<TriMeshType>::template DeletePerVertexAttribute<ScalarType>(m, dist);

  // The predicate above declines to split an edge whose crossing lands within `tolerance`
  // of an endpoint, on the grounds that the endpoint is close enough to serve as the cut
  // vertex. Completing that decision means putting it on the plane. Without this the cut is
  // planar only to within a fraction of an edge, which is invisible on an axis-aligned cut
  // through a regular mesh -- the crossings land on vertices exactly -- and routine on an
  // oblique one, where anything needing a planar outline, the cap below included, then has
  // nothing it can use.
  //
  // Only border vertices move, and only by less than the deviation the predicate already
  // accepted. A vertex on a boundary the mesh already had is left alone unless it happens
  // to lie that close to the plane, where the move is below the tolerance either way.
  if (anyAbove) {
    RequireFFAdjacency(m);
    UpdateTopology<TriMeshType>::FaceFace(m);
    // Only vertices on a boundary are candidates: an interior one was never a crossing.
    std::vector<bool> onBorder(m.vert.size(), false);
    // The tolerance is a fraction of an edge, so an edge is the scale to judge it against --
    // and it has to be any incident edge, not just a border one. The edge whose split was
    // declined ran to the discarded side, so by now it is gone; judging the deviation
    // against the short edges that happen to remain at that vertex underestimates the
    // scale, leaves the vertex unsnapped, and is why oblique cuts still would not close.
    std::vector<ScalarType> reach(m.vert.size(), ScalarType(0));
    for (size_t fi = 0; fi < m.face.size(); ++fi) {
      if (m.face[fi].IsD()) continue;
      for (int e = 0; e < 3; ++e) {
        VertexPointer a = m.face[fi].V0(e);
        VertexPointer b = m.face[fi].V1(e);
        const ScalarType len = (a->cP() - b->cP()).Norm();
        const size_t ia = size_t(Index(m, a)), ib = size_t(Index(m, b));
        if (len > reach[ia]) reach[ia] = len;
        if (len > reach[ib]) reach[ib] = len;
        if (face::IsBorder(m.face[fi], e)) {
          onBorder[ia] = true;
          onBorder[ib] = true;
        }
      }
    }
    for (size_t vi = 0; vi < m.vert.size(); ++vi) {
      if (m.vert[vi].IsD() || !onBorder[vi]) continue;
      const ScalarType s = n * m.vert[vi].cP() + d;
      if (math::Abs(s) <= tolerance * reach[vi])
        m.vert[vi].P() -= n * s;
    }
  }

  if (capCut && m.FN() > 0)
    CapPlanarBoundary(m, pl);
  return true;
}

/** \brief Fill the boundary loops of \a m that lie on \a plane, and only those.

  Used by ClipMeshWithPlane; separate because filling the planar boundary of a mesh that was
  cut some other way is the same job. All the qualifying loops are tessellated together
  under the even-odd rule, so concentric outlines cap as a ring and disjoint ones cap
  separately. Returns the number of loops filled; zero when there were none, when the mesh
  is not edge-manifold, or when the outline is not a set of simple loops.
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
  const ScalarType d = -pl.Offset();

  RequireFFAdjacency(m);
  UpdateTopology<TriMeshType>::FaceFace(m);
  if (Clean<TriMeshType>::CountNonManifoldEdgeFF(m) > 0) return 0;

  UpdateBounding<TriMeshType>::Box(m);
  const ScalarType tol = ScalarType(1e-5) * std::max(ScalarType(1e-6), m.bbox.Diag());
  const auto onPlane = [&](const VertexPointer v) {
    return math::Abs(n * v->cP() + d) <= tol;
  };

  // Walk each border loop once. Termination is on the visited marks and not on returning to
  // the starting Pos: a Pos carries a vertex as well as a face and an edge, so one circuit
  // comes back to the same edge with the other endpoint and the loop would be walked twice.
  std::vector<bool> walked(m.face.size() * 3, false);
  std::vector< std::vector<VertexPointer> > loops;
  for (size_t fi = 0; fi < m.face.size(); ++fi) {
    if (m.face[fi].IsD()) continue;
    for (int e = 0; e < 3; ++e) {
      if (!face::IsBorder(m.face[fi], e) || walked[fi * 3 + size_t(e)]) continue;
      std::vector<VertexPointer> loop;
      bool allOnPlane = true;
      face::Pos<FaceType> pos(&m.face[fi], e, m.face[fi].V(e));
      for (;;) {
        const size_t at = size_t(Index(m, pos.F())) * 3 + size_t(pos.E());
        if (walked[at]) break;
        walked[at] = true;
        loop.push_back(pos.V());
        allOnPlane = allOnPlane && onPlane(pos.V());
        pos.NextB();
      }
      if (!allOnPlane || loop.size() < 3) continue;
      // A loop that visits a vertex twice is pinched; splitting it into simple cycles is
      // guesswork, so say nothing was capped rather than cap it wrongly.
      std::vector<VertexPointer> sorted = loop;
      std::sort(sorted.begin(), sorted.end());
      if (std::unique(sorted.begin(), sorted.end()) != sorted.end()) return 0;
      loops.push_back(loop);
    }
  }
  if (loops.empty()) return 0;

  // Projected into the plane's own frame: the 3D entry point re-derives the plane and then
  // requires the points to sit on it far more tightly than float coordinates can express.
  CoordType u = (math::Abs(n[0]) < ScalarType(0.9)) ? CoordType(1, 0, 0) : CoordType(0, 1, 0);
  u = (u - n * (n * u)).Normalize();
  const CoordType v = (n ^ u).Normalize();

  std::vector< std::vector<Point2<ScalarType> > > contours;
  std::vector<VertexPointer> flat;
  for (size_t i = 0; i < loops.size(); ++i) {
    std::vector<Point2<ScalarType> > contour;
    for (size_t k = 0; k < loops[i].size(); ++k) {
      const CoordType p = loops[i][k]->cP();
      contour.push_back(Point2<ScalarType>(u * p, v * p));
      flat.push_back(loops[i][k]);
    }
    contours.push_back(contour);
  }

  std::vector<int> tri;
  if (!TessellatePlanarContours2(contours, tri) || tri.size() < 3) return 0;

  const size_t firstNew = m.face.size();
  for (size_t t = 0; t + 2 < tri.size(); t += 3) {
    const int a = tri[t], b = tri[t + 1], c = tri[t + 2];
    if (a == b || b == c || a == c) continue;
    if (a < 0 || b < 0 || c < 0) continue;
    if (size_t(a) >= flat.size() || size_t(b) >= flat.size() || size_t(c) >= flat.size()) continue;
    Allocator<TriMeshType>::AddFace(m, flat[size_t(a)], flat[size_t(b)], flat[size_t(c)]);
  }
  if (m.face.size() == firstNew) return 0;

  // The tessellator winds by the basis, which need not face the way the cut does: the cap
  // closes the side that was discarded, so it points against the plane normal.
  UpdateNormal<TriMeshType>::PerFaceNormalized(m);
  if (m.face[firstNew].cN() * n > 0) {
    for (size_t fi = firstNew; fi < m.face.size(); ++fi)
      if (!m.face[fi].IsD()) std::swap(m.face[fi].V(1), m.face[fi].V(2));
  }
  return int(loops.size());
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
