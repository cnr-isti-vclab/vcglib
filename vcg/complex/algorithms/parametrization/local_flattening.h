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

#ifndef VCG_LOCAL_FLATTENING
#define VCG_LOCAL_FLATTENING

#include <vcg/complex/algorithms/clip.h>
#include <vcg/complex/algorithms/create/platonic.h>
#include <vcg/complex/algorithms/vertex_interpolation.h>
#include <vcg/complex/algorithms/closest.h>
#include <vcg/complex/algorithms/geodesic.h>
#include <vcg/complex/algorithms/update/texture.h>
#include <vcg/complex/algorithms/parametrization/poisson_solver.h>
#include <vcg/space/index/grid_static_ptr.h>

namespace vcg {
namespace tri {

/** Flattens a strip of a surface around a polyline drawn on it into a rectangle, as the local
  flattening operator of Pietroni, Corsini, Cignoni and Scopigno, "An Interactive Local
  Flattening Operator to Support Digital Investigations on Artwork Surfaces", IEEE TVCG 17(12),
  2011, its first four steps:

  1. the faces within a geodesic distance of the polyline a little over the half-width are
     taken, the connected piece the polyline lies on, which must be a disk;
  2. a conformal map (PoissonSolver) with the polyline's points softly pulled to U = 0 and V
     equal to their arc length: its gradients make a cross field that follows the polyline.
     The paper's eq. 1 pulls U only; constraining V as well, so that the whole polyline
     becomes the straight V axis, removes the nearly conformal deformations that U alone
     leaves free and that rotate the field away from the polyline (on a developable cylinder,
     a helical polyline: 1.7% average length error with U only and the ends pinned, 0.02%
     with both and no pins);
  3. a map aligned to that field, so that a unit of UV is a unit on the surface;
  4. the result cropped to the rectangle |U| <= halfWidth, 0 <= V <= length, the polyline
     running along V from its first end, exactly (ClipMeshWithUVLine).

  The surface keeps its scale: lengths read in UV are lengths on the surface, up to the
  distortion that Stats reports.

  Two options change what the patch is made of, both through a geometry image: the cropped
  surface sampled at the nodes of a regular grid laid exactly on the rectangle.
  - gridStep > 0 remeshes the patch as that grid, in parametric space, instead of keeping the
    surface's own connectivity; a node with no surface under it, in a hole, is left out.
  - smoothSigma > 0 smooths its shape by a Gaussian over the grid, as the paper's section 4.1 does
    to separate the overall shape from the detail. Without gridStep, each vertex of the cropped
    surface moves to the smoothed shape at its own UV, on a grid of about its edge length.
  A smooth patch is what casting along its normals onto the original needs to read the detail.
*/
template <class MeshType>
class LocalFlattening
{
  typedef typename MeshType::ScalarType ScalarType;
  typedef typename MeshType::CoordType CoordType;
  typedef typename MeshType::VertexType VertexType;
  typedef typename MeshType::FaceType FaceType;
  typedef Point2<ScalarType> UV;
  typedef GridStaticPtr<FaceType, ScalarType> GridType;

public:
  struct Param
  {
    ScalarType halfWidth = 0;   ///< half the width of the patch, in the surface's units
    ScalarType stiffness = 10;  ///< weight of each polyline sample's pull to the V axis (k_st)
    ScalarType gridStep = 0;    ///< > 0: remesh the patch as a regular grid of this step in UV
    ScalarType smoothSigma = 0; ///< > 0: smooth its shape by a Gaussian of this deviation
    CallBackPos *cb = nullptr;
  };

  /// The distortion of the patch, as the paper's Table 1: ratios of 3D to UV edge lengths
  struct Stats
  {
    ScalarType length = 0;      ///< the patch's extent along V
    ScalarType minRatio = 0, maxRatio = 0;
    ScalarType avgError = 0;    ///< mean of |ratio - 1|
  };

  /** \param m the surface, with face-face and vertex-face adjacency, per-vertex quality, per-face
           marks and up to date face normals. Its selection and visited flags are overwritten.
      \param poly a single open polyline (its edges), in the frame of \a m and lying on it.
      \param patch receives the rectangle; it needs per-vertex and per-wedge texture
           coordinates, per-face curvature directions, per-face marks and face-face adjacency.
      Throws on a polyline that is not one open curve or is too far from the surface, and on
      a region that is not a disk. */
  static Stats Flatten(MeshType &m, const MeshType &poly, MeshType &patch, const Param &par)
  {
    RequireFFAdjacency(m);
    RequireVFAdjacency(m);
    RequirePerVertexQuality(m);
    if (!(par.halfWidth > 0)) throw std::invalid_argument("LocalFlattening: the half-width must be positive");

    // the polyline: samples about an edge of m apart, from its first end to its last, with
    // their arc length
    std::vector<CoordType> samples;
    std::vector<ScalarType> arc;
    Resample(poly, AverageEdge(m), samples, arc);
    const CoordType ends[2] = { samples.front(), samples.back() };
    if (par.cb) par.cb(5, "Selecting the region");

    // 1. the faces within reach of the polyline, the connected piece the polyline lies on
    {
      GridType grid;
      grid.Set(m.face.begin(), m.face.end());
      std::vector<VertexType *> seeds;
      for (const CoordType &p : samples) {
        CoordType bary;
        FaceType *f = Closest(m, grid, p, par.halfWidth, bary);
        seeds.push_back(f->V(int(std::max_element(&bary[0], &bary[0] + 3) - &bary[0])));
      }
      std::sort(seeds.begin(), seeds.end());
      seeds.erase(std::unique(seeds.begin(), seeds.end()), seeds.end());
      std::vector<VertexType *> within;
      EuclideanDistance<MeshType> dd;
      std::vector<ScalarType> quality;  // the geodesic writes the quality: keep the mesh's own
      ForEachVertex(m, [&](VertexType &v) { quality.push_back(v.Q()); });
      // a margin over the half-width, so that the crop, and not the selection, makes the border
      Geodesic<MeshType>::Compute(m, seeds, dd, par.halfWidth * ScalarType(1.25) + 2 * AverageEdge(m), &within);
      UpdateFlags<MeshType>::VertexClearV(m);
      for (VertexType *v : within) v->SetV();
      size_t qi = 0;
      ForEachVertex(m, [&](VertexType &v) { v.Q() = quality[qi++]; });
      CoordType bary;
      SelectReachable(m, Closest(m, grid, ends[0], par.halfWidth, bary));
    }
    UpdateSelection<MeshType>::VertexFromFaceLoose(m);
    Append<MeshType, MeshType>::MeshCopy(patch, m, true);
    UpdateSelection<MeshType>::Clear(patch);
    UpdateNormal<MeshType>::PerFaceNormalized(patch);

    PoissonSolver<MeshType> PS(patch);
    if (!PS.IsFeasible())
      throw std::runtime_error("LocalFlattening: the region around the polyline is not a disk; make the half-width smaller");

    // the samples on the patch, by face and barycentric coordinates
    GridType grid;
    grid.Set(patch.face.begin(), patch.face.end());
    std::vector<std::pair<FaceType *, CoordType>> onPatch(samples.size());
    for (size_t i = 0; i < samples.size(); ++i)
      onPatch[i].first = Closest(patch, grid, samples[i], par.halfWidth, onPatch[i].second);
    std::pair<FaceType *, CoordType> endOnPatch[2];
    for (int k = 0; k < 2; ++k)
      endOnPatch[k].first = Closest(patch, grid, ends[k], par.halfWidth, endOnPatch[k].second);

    // 2. the conformal map with the polyline pulled to U = 0, and the cross field it gives
    if (par.cb) par.cb(15, "Computing the cross field");
    // the constraints along the whole polyline fix the map; stiff pins at its ends on top of
    // them, as eq. 1 has, add nothing but a singularity in the field around each pin
    PS.Init();
    for (size_t i = 0; i < onPatch.size(); ++i) {
      PS.AddPointConstraint(onPatch[i].first, onPatch[i].second, 0, 0, par.stiffness);
      PS.AddPointConstraint(onPatch[i].first, onPatch[i].second, 1, arc[i], par.stiffness);
    }
    if (!PS.SolvePoisson()) throw std::runtime_error("LocalFlattening: the cross field system is singular");
    ForEachFace(patch, [](FaceType &f) {
      f.PD1().Import(Gradient(f, 0));
      f.PD2().Import(Gradient(f, 1));
    });

    // 3. the map aligned to that field, at the surface's scale
    if (par.cb) par.cb(50, "Parametrizing along the field");
    PS.Init(true);
    PS.FixPoint(endOnPatch[0].first, endOnPatch[0].second, UV(0, 0));
    if (!PS.SolvePoisson()) throw std::runtime_error("LocalFlattening: the field-aligned system is singular");

    // 4. the polyline's first end at the origin and its last on +V, then the crop
    if (par.cb) par.cb(85, "Cropping");
    const UV s = At(endOnPatch[0]), d = At(endOnPatch[1]) - s;
    Stats st;
    st.length = d.Norm();
    if (!(st.length > 0)) throw std::runtime_error("LocalFlattening: the polyline's ends map to the same point");
    const UV dir = d / st.length;
    ForEachVertex(patch, [&](VertexType &v) {
      const UV q = v.T().P() - s;
      v.T().P() = UV(dir[1] * q[0] - dir[0] * q[1], dir[0] * q[0] + dir[1] * q[1]);
    });
    ClipMeshWithUVLine(patch, UV(1, 0), -par.halfWidth);
    ClipMeshWithUVLine(patch, UV(-1, 0), -par.halfWidth);
    ClipMeshWithUVLine(patch, UV(0, 1), ScalarType(0));
    ClipMeshWithUVLine(patch, UV(0, -1), -st.length);
    if (patch.FN() == 0) throw std::runtime_error("LocalFlattening: nothing is left after cropping");
    if (par.gridStep > 0 || par.smoothSigma > 0)
      SampleGeometryImage(patch, par.halfWidth, st.length, par.gridStep > 0 ? par.gridStep : AverageEdge(patch),
               par.gridStep > 0, par.smoothSigma);
    UpdateTexture<MeshType>::WedgeTexFromVertexTex(patch);
    UpdateNormal<MeshType>::PerFaceNormalized(patch);
    UpdateBounding<MeshType>::Box(patch);

    // the distortion, over the edges of every face
    st.minRatio = std::numeric_limits<ScalarType>::max();
    int n = 0;
    ForEachFace(patch, [&](FaceType &f) {
      for (int k = 0; k < 3; ++k) {
        const ScalarType uv = Distance(f.WT(k).P(), f.WT((k + 1) % 3).P());
        if (uv <= 0) continue;
        const ScalarType r = Distance(f.P(k), f.P((k + 1) % 3)) / uv;
        st.minRatio = std::min(st.minRatio, r);
        st.maxRatio = std::max(st.maxRatio, r);
        st.avgError += std::abs(r - 1);
        ++n;
      }
    });
    if (n > 0) st.avgError /= n;
    if (par.cb) par.cb(100, "Done");
    return st;
  }

private:
  /// sample the patch on a grid of about step over |U| <= w, 0 <= V <= L, smooth the samples by
  /// sigma if positive, then make the patch that grid (remesh) or move its vertices onto it
  static void SampleGeometryImage(MeshType &patch, ScalarType w, ScalarType L, ScalarType step, bool remesh, ScalarType sigma)
  {
    const int nu = std::max(2, int(std::lround(2 * w / step)) + 1), nv = std::max(2, int(std::lround(L / step)) + 1);
    const double du = 2.0 * w / (nu - 1), dv = double(L) / (nv - 1);
    const size_t n = size_t(nu) * size_t(nv);
    std::vector<int> nodeFace(n, -1);
    std::vector<Point3d> nodeBary(n), pos(n, Point3d(0, 0, 0));
    std::vector<double> mask(n, 0);
    std::vector<Point3d> raw;
    // Each node takes the nearest point of the patch in UV: the triangle that holds it, or, for a
    // node just outside every triangle, the nearest one within reach, a twentieth of a step. That
    // covers a node on an edge shared by two faces, which rounding can leave outside both, and the
    // border, which the crop leaves within rounding of the rectangle, on either side; a hole in
    // the surface wider than that stays a hole.
    const double reach = 0.05;
    std::vector<double> nodeDist(n, reach);
    for (size_t fi = 0; fi < patch.face.size(); ++fi) {
      const FaceType &f = patch.face[fi];
      if (f.IsD()) continue;
      Point2d t[3];
      for (int k = 0; k < 3; ++k)
        t[k] = Point2d((double(f.cV(k)->cT().U()) + w) / du, double(f.cV(k)->cT().V()) / dv);
      const double area = (t[1] - t[0]) ^ (t[2] - t[0]);
      if (area == 0) continue;
      const int i0 = std::max(0, int(std::floor(std::min({t[0][0], t[1][0], t[2][0]}) - reach)));
      const int i1 = std::min(nu - 1, int(std::ceil(std::max({t[0][0], t[1][0], t[2][0]}) + reach)));
      const int j0 = std::max(0, int(std::floor(std::min({t[0][1], t[1][1], t[2][1]}) - reach)));
      const int j1 = std::min(nv - 1, int(std::ceil(std::max({t[0][1], t[1][1], t[2][1]}) + reach)));
      for (int j = j0; j <= j1; ++j)
        for (int i = i0; i <= i1; ++i) {
          const size_t q = size_t(j) * nu + i;
          if (nodeDist[q] == 0) continue;
          const Point2d c(i, j);
          Point3d b;
          for (int k = 0; k < 3; ++k) b[k] = ((t[(k + 1) % 3] - c) ^ (t[(k + 2) % 3] - c)) / area;
          double d = 0;
          if (b[0] < 0 || b[1] < 0 || b[2] < 0) {  // outside: the nearest point of the triangle
            Point2d best;
            d = std::numeric_limits<double>::max();
            for (int k = 0; k < 3; ++k) {
              const Point2d a0 = t[k], e = t[(k + 1) % 3] - t[k];
              const double s = std::clamp(((c - a0) * e) / (e * e), 0.0, 1.0);
              const Point2d p = a0 + e * s;
              if (Distance(p, c) < d) { d = Distance(p, c); best = p; }
            }
            if (!(d < nodeDist[q])) continue;
            for (int k = 0; k < 3; ++k) b[k] = ((t[(k + 1) % 3] - best) ^ (t[(k + 2) % 3] - best)) / area;
            for (int k = 0; k < 3; ++k) b[k] = std::max(0.0, b[k]);
            b /= b[0] + b[1] + b[2];
          }
          nodeDist[q] = d;
          nodeFace[q] = int(fi);
          nodeBary[q] = b;
          pos[q] = Point3d(0, 0, 0);
          for (int k = 0; k < 3; ++k) pos[q] += Point3d::Construct(f.cP(k)) * b[k];
          mask[q] = 1;
        }
    }
    raw = pos;
    if (sigma > 0) {
      std::vector<Point3d> S = pos;
      std::vector<double> w8 = mask;
      Blur(nu, nv, sigma / du, sigma / dv, S, w8);
      for (size_t q = 0; q < n; ++q)
        if (w8[q] > 0) S[q] /= w8[q];
      auto at = [&](int i, int j) { return std::clamp(j, 0, nv - 1) * size_t(nu) + std::clamp(i, 0, nu - 1); };
      for (int j = 0; j < nv; ++j)
        for (int i = 0; i < nu; ++i) {
          const size_t q = size_t(j) * nu + i;
          if (!(mask[q] > 0)) continue;
          const Point3d nrm = ((S[at(i + 1, j)] - S[at(i - 1, j)]) ^ (S[at(i, j + 1)] - S[at(i, j - 1)])).normalized();
          pos[q] -= nrm * ((pos[q] - S[q]) * nrm);
        }
    }

    if (!remesh) {  // each vertex moved as the grid samples around it were, bilinearly
      ForEachVertex(patch, [&](VertexType &v) {
        const double x = std::clamp((double(v.T().U()) + w) / du, 0.0, nu - 1.0);
        const double y = std::clamp(double(v.T().V()) / dv, 0.0, nv - 1.0);
        const int i = std::min(int(x), nu - 2), j = std::min(int(y), nv - 2);
        Point3d sp(0, 0, 0);
        double sw = 0;
        for (int dj = 0; dj < 2; ++dj)
          for (int di = 0; di < 2; ++di) {
            const size_t q = size_t(j + dj) * nu + i + di;
            const double k = (di ? x - i : 1 - (x - i)) * (dj ? y - j : 1 - (y - j));
            sp += (pos[q] - raw[q]) * (k * mask[q]);
            sw += mask[q] * k;
          }
        if (sw > 0) v.P() += CoordType::Construct(sp / sw);
      });
      return;
    }

    // the grid's vertices, with the attributes of the surface where they lie, after the old ones
    const size_t oldVerts = patch.vert.size();
    std::vector<int> grid(n, -1);
    int count = 0;
    for (size_t q = 0; q < n; ++q)
      if (nodeFace[q] >= 0) grid[q] = count++;
    Allocator<MeshType>::AddVertices(patch, size_t(count));
    for (size_t q = 0; q < n; ++q) {
      if (grid[q] < 0) continue;
      VertexType &gv = patch.vert[oldVerts + size_t(grid[q])];
      const FaceType &f = patch.face[size_t(nodeFace[q])];
      VertexInterpolator<MeshType>::Barycentric(patch, gv, *f.cV(0), *f.cV(1), *f.cV(2), CoordType::Construct(nodeBary[q]));
      gv.P() = CoordType::Construct(pos[q]);
      gv.T().P() = Point2<ScalarType>(ScalarType(-w + (q % nu) * du), ScalarType((q / nu) * dv));
    }
    ForEachFace(patch, [&](FaceType &f) { Allocator<MeshType>::DeleteFace(patch, f); });
    for (size_t k = 0; k < oldVerts; ++k)
      if (!patch.vert[k].IsD()) Allocator<MeshType>::DeleteVertex(patch, patch.vert[k]);
    Allocator<MeshType>::CompactEveryVector(patch);
    SparseFaceGrid(patch, grid, nu, nv);
    if (patch.FN() == 0) throw std::runtime_error("LocalFlattening: the grid is too coarse for the patch");
  }

  /// separable Gaussian of positions weighted by a mask, in place (normalized convolution)
  static void Blur(int w, int h, double sx, double sy, std::vector<Point3d> &a, std::vector<double> &c)
  {
    for (int pass = 0; pass < 2; ++pass) {
      const double sg = pass == 0 ? sx : sy;
      const int r = std::max(1, int(std::ceil(3 * sg)));
      std::vector<double> k(2 * r + 1);
      for (int i = -r; i <= r; ++i) k[i + r] = std::exp(-0.5 * i * i / (sg * sg));
      const int len = pass == 0 ? w : h, lines = pass == 0 ? h : w;
      std::vector<Point3d> la(len);
      std::vector<double> lc(len);
      for (int l = 0; l < lines; ++l) {
        auto at = [&](int x) { return pass == 0 ? size_t(l) * w + x : size_t(x) * w + l; };
        for (int x = 0; x < len; ++x) { la[x] = a[at(x)]; lc[x] = c[at(x)]; }
        for (int x = 0; x < len; ++x) {
          Point3d sa(0, 0, 0);
          double sc = 0;
          for (int i = std::max(0, x - r); i <= std::min(len - 1, x + r); ++i) {
            sa += la[i] * k[i - x + r];
            sc += lc[i] * k[i - x + r];
          }
          a[at(x)] = sa;
          c[at(x)] = sc;
        }
      }
    }
  }

  static ScalarType AverageEdge(const MeshType &m)
  {
    double sum = 0; int n = 0;
    ForEachFace(m, [&](const FaceType &f) {
      for (int k = 0; k < 3; ++k) { sum += Distance(f.cP(k), f.cP((k + 1) % 3)); ++n; }
    });
    return n ? ScalarType(sum / n) : ScalarType(0);
  }

  /// samples of the polyline about step apart, walking it from its first end (the end with the
  /// lower index) to the other, and their arc length; throws if it is not one open curve
  static void Resample(const MeshType &poly, ScalarType step, std::vector<CoordType> &samples, std::vector<ScalarType> &arc)
  {
    std::vector<std::vector<int>> next(poly.vert.size());
    int edges = 0;
    ForEachEdge(poly, [&](const typename MeshType::EdgeType &e) {
      const int a = int(Index(poly, e.cV(0))), b = int(Index(poly, e.cV(1)));
      next[a].push_back(b);
      next[b].push_back(a);
      ++edges;
    });
    int first = -1, ends = 0;
    bool branches = false;
    for (size_t i = 0; i < next.size(); ++i) {
      if (next[i].size() > 2) branches = true;
      if (next[i].size() == 1) { ++ends; if (first < 0) first = int(i); }
    }
    // walk the chain; it is one open curve when the walk takes every edge
    int walked = 0;
    for (int prev = -1, cur = first; !branches && cur >= 0;) {
      int to = -1;
      for (int n : next[cur]) if (n != prev) to = n;
      if (to < 0) break;
      const CoordType p0 = poly.vert[cur].cP(), p1 = poly.vert[to].cP();
      const ScalarType len = Distance(p0, p1);
      const int n = std::max(1, int(std::ceil(len / step)));
      const ScalarType arc0 = arc.empty() ? 0 : arc.back();
      for (int k = (samples.empty() ? 0 : 1); k <= n; ++k) {
        samples.push_back(p0 + (p1 - p0) * (ScalarType(k) / n));
        arc.push_back(arc0 + len * ScalarType(k) / n);
      }
      ++walked;
      prev = cur;
      cur = to;
    }
    if (branches || ends != 2 || walked != edges || edges == 0)
      throw std::invalid_argument("LocalFlattening: the polyline must be a single open curve, without branches");
    if (!(arc.back() > 0)) throw std::invalid_argument("LocalFlattening: the polyline has no length");
  }

  /// the face of m closest to p, and the barycentric coordinates of the closest point in it
  static FaceType *Closest(MeshType &m, GridType &grid, const CoordType &p, ScalarType maxDist, CoordType &bary)
  {
    ScalarType dist;
    CoordType q;
    FaceType *f = GetClosestFaceBase(m, grid, p, maxDist, dist, q);
    if (!f) throw std::invalid_argument("LocalFlattening: the polyline is farther from the surface than the half-width");
    InterpolationParameters(*f, f->cN(), q, bary);
    return f;
  }

  /// select the faces with three visited vertices that can be reached from start through such faces
  static void SelectReachable(MeshType &m, FaceType *start)
  {
    UpdateSelection<MeshType>::FaceClear(m);
    auto inside = [](const FaceType *f) { return f->cV(0)->IsV() && f->cV(1)->IsV() && f->cV(2)->IsV(); };
    if (!inside(start)) return;
    std::vector<FaceType *> stack(1, start);
    start->SetS();
    while (!stack.empty()) {
      FaceType *f = stack.back();
      stack.pop_back();
      for (int k = 0; k < 3; ++k) {
        FaceType *g = f->FFp(k);
        if (g != f && !g->IsS() && inside(g)) { g->SetS(); stack.push_back(g); }
      }
    }
  }

  /// the 3D gradient on a face of its per-vertex U (coord 0) or V (coord 1)
  static CoordType Gradient(const FaceType &f, int coord)
  {
    const ScalarType a2 = DoubleArea(f);
    if (!(a2 > 0)) return CoordType(0, 0, 0);
    CoordType n = f.cN();
    n.Normalize();
    CoordType g(0, 0, 0);
    for (int i = 0; i < 3; ++i)
      g += (n ^ (f.cP((i + 2) % 3) - f.cP((i + 1) % 3))) * f.cV(i)->cT().P()[coord];
    return g / a2;
  }

  static UV At(const std::pair<FaceType *, CoordType> &s)
  {
    return s.first->V(0)->T().P() * s.second[0] + s.first->V(1)->T().P() * s.second[1] +
           s.first->V(2)->T().P() * s.second[2];
  }
};

} // namespace tri
} // namespace vcg
#endif
