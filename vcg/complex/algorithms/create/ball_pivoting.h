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
#ifndef BALL_PIVOTING_H
#define BALL_PIVOTING_H

#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/update/bounding.h>
#include <vcg/complex/algorithms/update/normal.h>

#include <algorithm>
#include <cmath>
#include <deque>
#include <unordered_map>
#include <vector>

/* Ball Pivoting Algorithm, after
   Bernardini, Mittleman, Rushmeier, Silva, Taubin,
   "The ball-pivoting algorithm for surface reconstruction", IEEE TVCG 1999.

   A port of BPA.jl (github.com/ctsilva/BPA.jl) by way of the C++ implementation
   of Bernhard Gruber (github.com/bernhardmgruber/bpa, Boost Software License)
   as extended by Claudio Silva (github.com/ctsilva/bpa). It produces the same
   triangles as BPA.jl on the same input:

   - every point within reach of the ball is a candidate for a pivot, the
     opposite vertex of the edge included, so the ball of every triangle that
     is built is empty (the property that defines the algorithm);
   - the first contact is found from the position of the ball centre in the
     pivot plane, without trigonometry; simultaneous hits (cospherical points,
     as on a lattice) are resolved deterministically;
   - the front is a FIFO queue, so the surface grows breadth-first from each
     seed;
   - a seed must face along the normals of its three vertices, and may reuse
     vertices already on the front (its edges are glued at once); a pivot
     triangle is tested against the normal of the point the ball lands on
     only, a zero dot product (a point without a normal) being accepted;
   - points are found through a uniform grid of cells of side 2*radius;
   - existing faces of the mesh are kept: their boundary edges resume pivoting
     when the triangle they belong to admits an empty ball of the given radius
     (section 4.6 of the paper), so running the algorithm again with a larger
     radius on its own output closes the gaps the smaller ball could not cross.

   A cloud without any vertex normal is accepted: its seeds are then oriented
   away from the barycentre of the cloud, as in the earlier implementation.

   All the geometry is computed in double precision whatever the scalar type
   of the mesh.

   Usage:
      tri::BallPivoting<MyMesh> pivot(m, radius);   // radius 0: a guess from the sampling
      pivot.BuildMesh(cb);
*/
namespace vcg {
namespace tri {

namespace ballpivoting {

typedef Point3<double> Vec;
typedef Point3<int>    IVec;

inline Vec unit(Vec v) { v.Normalize(); return v; } // a zero vector stays zero

struct MeshEdge;

struct MeshPoint {
  Vec pos;
  Vec normal;
  int index;                    // vertex index in the mesh, -1 for a deleted vertex
  bool used;
  std::vector<MeshEdge *> edges; // every front edge ever incident, live or not
};

enum EdgeStatus { ACTIVE, INNER, BOUNDARY };

// A directed front edge a -> b of the triangle (a, b, opposite), wound like the triangle,
// with the centre of the ball that built the triangle. Loops of the front are a doubly
// linked list through prev/next. An inner edge has its two triangles and only records that
// the undirected edge is closed.
struct MeshEdge {
  MeshPoint *a;
  MeshPoint *b;
  MeshPoint *opposite;
  Vec center;
  MeshEdge *prev;
  MeshEdge *next;
  EdgeStatus status;
  MeshEdge(MeshPoint *a_, MeshPoint *b_, MeshPoint *o_, const Vec &c, EdgeStatus s = ACTIVE):
    a(a_), b(b_), opposite(o_), center(c), prev(0), next(0), status(s) {}
};

struct MeshFace {
  MeshPoint *p[3];
  MeshFace() { p[0] = p[1] = p[2] = 0; }
  MeshFace(MeshPoint *a, MeshPoint *b, MeshPoint *c) { p[0] = a; p[1] = b; p[2] = c; }
  Vec normal() const { return unit((p[1]->pos - p[0]->pos) ^ (p[2]->pos - p[0]->pos)); }
};

// A cell keeps a copy of each position next to the point, so that scanning a cell for the
// points within reach reads memory in order and follows the pointer only for those.
struct CellEntry {
  Vec pos;
  MeshPoint *p;
};
typedef std::vector<CellEntry> Cell;

struct CellKey {
  int x, y, z;
  bool operator==(const CellKey &o) const { return x == o.x && y == o.y && z == o.z; }
};
struct CellKeyHash {
  size_t operator()(const CellKey &k) const {
    return size_t(k.x) * 73856093u ^ size_t(k.y) * 19349663u ^ size_t(k.z) * 83492791u;
  }
};

// Uniform grid of cells of side 2 * radius: every point a ball touching a point in a cell
// can touch lies in that cell or one of its 26 neighbours. Dense when the bounding box
// needs few cells per point, hashed otherwise (a radius that is tiny against the box).
// The points of a cell are kept in index order, and the cells are enumerated with x
// fastest, so that queries return points in the same order as BPA.jl.
class Grid {
public:
  Grid(std::vector<MeshPoint> &points, double radius): cellSize(radius * 2), dense(true) {
    bool first = true;
    for (size_t i = 0; i < points.size(); i++) {
      if (points[i].index < 0) continue;
      if (first) { lower = points[i].pos; first = false; }
      for (int k = 0; k < 3; k++) lower[k] = std::min(lower[k], points[i].pos[k]);
    }
    if (first) lower = Vec(0, 0, 0);
    // the box is sized from the cell coordinates themselves, so no point rounds outside
    dims = IVec(1, 1, 1);
    std::vector<IVec> coords(points.size());
    for (size_t i = 0; i < points.size(); i++) {
      if (points[i].index < 0) continue;
      coords[i] = cellIndex(points[i].pos);
      for (int k = 0; k < 3; k++) dims[k] = std::max(dims[k], coords[i][k] + 1);
    }
    const double total = double(dims[0]) * double(dims[1]) * double(dims[2]);
    dense = total <= std::max(8.0 * double(points.size()), 4096.0);
    if (dense) cells.resize(size_t(total));
    for (size_t i = 0; i < points.size(); i++) {
      if (points[i].index < 0) continue;
      CellEntry e;
      e.pos = points[i].pos;
      e.p = &points[i];
      get(coords[i]).push_back(e);
    }
    // the non-empty cells in linear order (x fastest), for the seed scan
    if (dense) {
      for (size_t i = 0; i < cells.size(); i++)
        if (!cells[i].empty()) order.push_back(&cells[i]);
    } else {
      std::vector<std::pair<CellKey, Cell *> > tmp;
      for (SparseMap::iterator it = sparse.begin(); it != sparse.end(); ++it)
        tmp.push_back(std::make_pair(it->first, &it->second));
      std::sort(tmp.begin(), tmp.end(), keyLess);
      for (size_t i = 0; i < tmp.size(); i++) order.push_back(tmp[i].second);
    }
  }

  IVec cellIndex(const Vec &p) const {
    IVec r;
    for (int k = 0; k < 3; k++) {
      double v = std::floor((p[k] - lower[k]) / cellSize);
      if (v < -2e9) v = -2e9;
      if (v > 2e9) v = 2e9;
      r[k] = int(v);
    }
    return r;
  }

  // The non-empty cells, in a fixed order.
  const std::vector<Cell *> &cellsInOrder() const { return order; }

  // The points within `r` (at most the cell size) of `point`, in grid order.
  void neighbors(const Vec &point, double r, std::vector<CellEntry> &result) const {
    result.clear();
    const IVec c = cellIndex(point);
    const double r2 = r * r;
    for (int zo = -1; zo <= 1; zo++)
      for (int yo = -1; yo <= 1; yo++)
        for (int xo = -1; xo <= 1; xo++) {
          const Cell *cell = find(IVec(c[0] + xo, c[1] + yo, c[2] + zo));
          if (!cell) continue;
          for (size_t i = 0; i < cell->size(); i++)
            if (((*cell)[i].pos - point).SquaredNorm() <= r2) result.push_back((*cell)[i]);
        }
  }

  // No point other than up to three given ones lies strictly inside the ball. Points within
  // a relative 1e-9 of the sphere count as outside, so that the points the ball touches, and
  // points exactly cospherical with them, do not fail the test.
  bool ballIsEmpty(const Vec &center, double radius, const MeshPoint *e0, const MeshPoint *e1, const MeshPoint *e2) const {
    const double r2 = radius * (1 - 1e-9) * radius * (1 - 1e-9);
    const IVec c = cellIndex(center);
    for (int zo = -1; zo <= 1; zo++)
      for (int yo = -1; yo <= 1; yo++)
        for (int xo = -1; xo <= 1; xo++) {
          const Cell *cell = find(IVec(c[0] + xo, c[1] + yo, c[2] + zo));
          if (!cell) continue;
          for (size_t i = 0; i < cell->size(); i++) {
            const CellEntry &e = (*cell)[i];
            if ((e.pos - center).SquaredNorm() < r2 && e.p != e0 && e.p != e1 && e.p != e2) return false;
          }
        }
    return true;
  }

  double cellSize;

private:
  typedef std::unordered_map<CellKey, Cell, CellKeyHash> SparseMap;

  static bool keyLess(const std::pair<CellKey, Cell *> &a, const std::pair<CellKey, Cell *> &b) {
    if (a.first.z != b.first.z) return a.first.z < b.first.z;
    if (a.first.y != b.first.y) return a.first.y < b.first.y;
    return a.first.x < b.first.x;
  }

  const Cell *find(const IVec &i) const {
    if (i[0] < 0 || i[1] < 0 || i[2] < 0) return 0;
    if (dense) {
      if (i[0] >= dims[0] || i[1] >= dims[1] || i[2] >= dims[2]) return 0;
      return &cells[(size_t(i[2]) * dims[1] + i[1]) * dims[0] + i[0]];
    }
    CellKey k = {i[0], i[1], i[2]};
    SparseMap::const_iterator it = sparse.find(k);
    return it == sparse.end() ? 0 : &it->second;
  }

  Cell &get(const IVec &i) {
    if (dense) return cells[(size_t(i[2]) * dims[1] + i[1]) * dims[0] + i[0]];
    CellKey k = {i[0], i[1], i[2]};
    return sparse[k];
  }

  Vec lower;
  bool dense;
  IVec dims;
  std::vector<Cell> cells;
  SparseMap sparse;
  std::vector<Cell *> order;
};

// Centre of the ball of the given radius through the three points of `f`, on the side of
// its normal; false if the circumradius is larger than the ball.
inline bool computeBallCenter(const MeshFace &f, double radius, Vec &center) {
  const Vec ac = f.p[2]->pos - f.p[0]->pos;
  const Vec ab = f.p[1]->pos - f.p[0]->pos;
  const Vec abXac = ab ^ ac;
  const double nn = abXac.dot(abXac);
  if (nn <= 1e-20 * ab.dot(ab) * ac.dot(ac)) // degenerate triangle
    return false;
  const Vec toCircumCircleCenter = ((abXac ^ ab) * ac.dot(ac) + (ac ^ abXac) * ab.dot(ab)) / (2 * nn);
  const Vec circumCircleCenter = f.p[0]->pos + toCircumCircleCenter;
  const double heightSquared = radius * radius - toCircumCircleCenter.dot(toCircumCircleCenter);
  if (heightSquared < 0)
    return false;
  center = circumCircleCenter + f.normal() * std::sqrt(heightSquared);
  return true;
}

// ---- the front predicates (section 4.4): an edge is live while it is active or boundary.

inline bool isLive(const MeshEdge *e) { return e->status != INNER; }

inline bool onFront(const MeshPoint *p) {
  for (size_t i = 0; i < p->edges.size(); i++)
    if (isLive(p->edges[i])) return true;
  return false;
}

// In the mesh with a complete fan of triangles: no more triangles may touch it.
inline bool isInterior(const MeshPoint *p) { return p->used && !onFront(p); }

// The front holds the directed edge i -> j.
inline bool hasEdge(const MeshPoint *i, const MeshPoint *j) {
  for (size_t k = 0; k < i->edges.size(); k++) {
    const MeshEdge *e = i->edges[k];
    if (isLive(e) && e->a == i && e->b == j) return true;
  }
  return false;
}

// The undirected edge {i, j} already has two triangles.
inline bool isClosed(const MeshPoint *i, const MeshPoint *j) {
  for (size_t k = 0; k < i->edges.size(); k++) {
    const MeshEdge *e = i->edges[k];
    if (e->status == INNER && ((e->a == i && e->b == j) || (e->a == j && e->b == i))) return true;
  }
  return false;
}

// ---- ball pivoting (section 4.3)

// Angle below which a candidate counts as touching the ball in its initial position.
const double TOUCH_TOLERANCE = 1e-6;
// Angle within which two hits count as simultaneous (cospherical points, as on a lattice).
const double TIE_TOLERANCE = 1e-7;

// The ball centre pivoting around edge (a, b) moves on the circle m + r (cos t u + sin t v)
// in the plane perpendicular to the edge through its midpoint m: u points from m to the
// current centre, and v = (b - a) x u is the direction in which the ball leaves the
// current triangle.
struct PivotFrame {
  Vec m, a, u, v;
  double r;
};

inline bool pivotFrame(const MeshEdge *e, PivotFrame &fr) {
  const Vec m = (e->a->pos + e->b->pos) / 2.0;
  Vec a = e->b->pos - e->a->pos;
  const double la = a.Norm();
  if (la == 0) return false;
  a /= la;
  Vec w = e->center - m;
  w -= a * w.dot(a); // numerical drift along the edge
  const double r = w.Norm();
  if (r <= 1e-12 * la) return false;
  fr.m = m;
  fr.a = a;
  fr.u = w / r;
  fr.v = a ^ fr.u;
  fr.r = r;
  return true;
}

// A position of the ball centre on the pivot circle, as coordinates (x, y) in the frame
// (u, v): the centre is m + x u + y v with x^2 + y^2 = r^2, and the pivot angle is the
// argument of (x, y) in [0, 2 pi). Contacts are compared by angle without evaluating it.
struct Contact {
  double x, y;
  Contact(): x(0), y(0) {}
  Contact(double x_, double y_): x(x_), y(y_) {}
};

inline Vec centerAt(const PivotFrame &fr, const Contact &c) { return fr.m + fr.u * c.x + fr.v * c.y; }

// The angle of c lies in [pi, 2 pi).
inline bool lowerHalf(const Contact &c) { return c.y < 0 || (c.y == 0 && c.x < 0); }

// The pivot angle of p is smaller than that of q: an angle in [0, pi) precedes any in
// [pi, 2 pi); within one half the two differ by less than pi, and the sign of the 2-D
// cross product decides.
inline bool angleLess(const Contact &p, const Contact &q) {
  const bool hp = lowerHalf(p);
  const bool hq = lowerHalf(q);
  return hp == hq ? p.x * q.y - p.y * q.x > 0 : hq;
}

// The pivot angles of p and q differ by at most the angle whose sine is sinTol, r2 being
// the squared radius of the pivot circle. Two contacts on either side of angle 0 are not a
// tie: the ball reaches one at once and the other after a full turn.
inline bool angleTie(const Contact &p, const Contact &q, double r2, double sinTol) {
  return p.x * q.x + p.y * q.y > 0 && std::fabs(p.x * q.y - p.y * q.x) <= r2 * sinTol &&
         (lowerHalf(p) == lowerHalf(q) || p.x < 0);
}

// The contact is within TOUCH_TOLERANCE of angle 0.
inline bool touching(const Contact &c, double r) { return c.x > 0 && std::fabs(c.y) < r * std::sin(TOUCH_TOLERANCE); }

// The ball centre at the first contact of the pivoting ball with x, or false if the ball
// never reaches x. With d = x - m and (d_u, d_v) its components in the pivot plane,
// |centre - x|^2 = radius^2 is the line d_u X + d_v Y = r K with K = (r^2 + |d|^2 -
// radius^2) / (2 r), whose intersections with the circle of radius r are F +- h (-d_v,
// d_u) / R: F = (r K / R^2)(d_u, d_v) is the foot of the perpendicular from the origin
// and h = r sqrt(1 - (K / R)^2) the half-chord. A contact at (nearly) zero means x
// touches the initial ball: if the ball is moving into x the hit is immediate, else that
// contact is ignored and the other one, where the ball comes back to x from the far
// side, is used. The opposite vertex of the edge is such a point.
inline bool pivotContact(const PivotFrame &fr, const Vec &x, double radius, Contact &out) {
  const Vec d = x - fr.m;
  const double du = d.dot(fr.u);
  const double dv = d.dot(fr.v);
  const double R = std::sqrt(du * du + dv * dv);
  if (R <= 1e-12 * radius) return false;
  const double r = fr.r;
  double ratio = (r * r + d.dot(d) - radius * radius) / (2 * r) / R;
  if (std::fabs(ratio) > 1 + 1e-9) return false;
  if (ratio > 1) ratio = 1;
  if (ratio < -1) ratio = -1;
  const double s = r * ratio / R;
  const double h = r * std::sqrt(1 - ratio * ratio) / R;
  const Contact ca(s * du - h * dv, s * dv + h * du); // phi + alpha
  const Contact cb(s * du + h * dv, s * dv - h * du); // phi - alpha
  const bool ta = touching(ca, r);
  const bool tb = touching(cb, r);
  if (ta || tb) {
    if (dv > 0) { out = Contact(r, 0); return true; } // the centre moves along v: into x
    if (ta && tb) return false;
    out = ta ? cb : ca;
    return true;
  }
  out = angleLess(cb, ca) ? cb : ca;
  return true;
}

} // namespace ballpivoting

template <class MESH> class BallPivoting {
public:
  typedef typename MESH::VertexType     VertexType;
  typedef typename MESH::FaceType       FaceType;
  typedef typename MESH::FaceIterator   FaceIterator;
  typedef typename MESH::ScalarType     ScalarType;
  typedef typename MESH::CoordType      CoordType;

  typedef ballpivoting::Vec       Vec;
  typedef ballpivoting::MeshPoint MeshPoint;
  typedef ballpivoting::MeshEdge  MeshEdge;
  typedef ballpivoting::MeshFace  MeshFace;
  typedef ballpivoting::CellEntry CellEntry;
  typedef ballpivoting::Cell      Cell;
  typedef ballpivoting::Grid      Grid;

  double radius;         // radius of the ball
  double min_edge;       // accepted and ignored: kept for source compatibility (was the clustering radius)
  double max_angle;      // cos of the largest dihedral angle allowed between a new triangle and
                         // the one it is pivoted from; -1 (the default) disables the test
  int seed_neighbors;    // only the nearest this many neighbours of a seed candidate are
                         // paired (0: all within 2 * radius, the paper's unbounded search)

  // stats of the last BuildMesh
  int seeds;             // seed triangles
  int pivots;            // pivots attempted
  int reactivated;       // boundary edges of the existing mesh that resumed pivoting

  // radius == 0: a guess from the bounding box and the number of points.
  // clustering: ignored, see min_edge.
  // angle: dihedral limit in radians; M_PI (the default) or more disables it.
  BallPivoting(MESH &_mesh, double _radius = 0, double clustering = 0.2, double angle = M_PI):
    radius(_radius), min_edge(clustering), max_angle(angle >= M_PI ? -1.0 : std::cos(angle)),
    seed_neighbors(100), seeds(0), pivots(0), reactivated(0),
    mesh(_mesh), qhead(0), grid(0), seedCursor(0), cloudHasNormals(false) {

    UpdateBounding<MESH>::Box(mesh);
    if (radius == 0) // radius ==0 means that an auto guess should be attempted.
      radius = std::sqrt(double(mesh.bbox.Diag()) * double(mesh.bbox.Diag()) / std::max(1, mesh.vn));

    // the points
    points.resize(mesh.vert.size());
    barycenter = Vec(0, 0, 0);
    int n = 0;
    for (size_t i = 0; i < mesh.vert.size(); i++) {
      MeshPoint &p = points[i];
      const VertexType &v = mesh.vert[i];
      p.index = v.IsD() ? -1 : int(i);
      p.used = v.IsD();
      p.pos = Vec(double(v.cP()[0]), double(v.cP()[1]), double(v.cP()[2]));
      p.normal = Vec(0, 0, 0);
      if (!v.IsD()) {
        if (HasPerVertexNormal(mesh))
          p.normal = ballpivoting::unit(Vec(double(v.cN()[0]), double(v.cN()[1]), double(v.cN()[2])));
        if (p.normal.SquaredNorm() > 0) cloudHasNormals = true;
        barycenter += p.pos;
        n++;
      }
    }
    if (n) barycenter /= n;

    grid = new Grid(points, radius);

    importExistingFaces();
  }

  ~BallPivoting() { delete grid; }

  // Pivot until every reachable edge has been tried and no seed is left, adding the
  // triangles to the mesh. `call` is told the progress every `interval` triangles.
  void BuildMesh(CallBackPos call = NULL, int interval = 512) {
    using namespace ballpivoting;
    if (call) (*call)(0, "Ball pivoting");
    const double expected = std::max(1.0, 2.0 * mesh.vn);
    int added = 0;

    // Fig. 5 of the paper: pivot until the front is exhausted, then seed again among the
    // points still unused, until no seed is left. The queue may start with the edges of
    // the existing mesh that importExistingFaces() reactivated.
    for (;;) {
      MeshEdge *e_ij;
      while ((e_ij = getActiveEdge()) != 0) {
        pivots++;
        MeshPoint *k = 0;
        Vec center;
        if (ballPivot(e_ij, k, center) && canAddTriangle(e_ij, k)) {
          addFace(e_ij->a->index, k->index, e_ij->b->index);
          join(e_ij, k, center);
          added++;
          if (call && interval > 0 && added % interval == 0) {
            const int perc = int(100.0 * added / expected);
            (*call)(std::min(perc, 99), "Adding Faces");
          }
        } else {
          e_ij->status = BOUNDARY;
        }
      }
      MeshFace seed;
      Vec ballCenter;
      if (!findSeedTriangle(seed, ballCenter)) break;
      seeds++;
      addFace(seed.p[0]->index, seed.p[1]->index, seed.p[2]->index);
      addSeed(seed, ballCenter);
      added++;
    }
    if (call) (*call)(100, "Ball pivoting");
  }

private:
  MESH &mesh;
  std::vector<MeshPoint> points;
  std::deque<MeshEdge> edges;     // stable addresses, in creation order
  std::vector<MeshEdge *> queue;  // FIFO of active edges; entries no longer active are skipped
  size_t qhead;
  Grid *grid;
  size_t seedCursor;
  bool cloudHasNormals;           // else seeds are oriented away from the barycentre
  Vec barycenter;
  std::vector<CellEntry> neighborhood; // reused buffers
  std::vector<CellEntry> seedPairs;

  void addFace(int v0, int v1, int v2) {
    FaceIterator fi = Allocator<MESH>::AddFace(mesh, size_t(v0), size_t(v1), size_t(v2));
    if (FaceType::HasNormal())
      fi->N() = TriangleNormal(*fi).Normalize();
    if (HasVFAdjacency(mesh)) {
      for (int j = 0; j < 3; ++j) {
        (*fi).VFp(j) = (*fi).V(j)->VFp();
        (*fi).VFi(j) = (*fi).V(j)->VFi();
        (*fi).V(j)->VFp() = &(*fi);
        (*fi).V(j)->VFi() = j;
      }
    }
  }

  // ---- the front

  MeshEdge *getActiveEdge() {
    using namespace ballpivoting;
    while (qhead < queue.size()) {
      MeshEdge *e = queue[qhead++];
      if (e->status == ACTIVE) return e;
    }
    queue.clear();
    qhead = 0;
    return 0;
  }

  // A new active edge a -> b of the triangle (a, b, o), registered at both endpoints and
  // queued. Loop links are left to the caller.
  MeshEdge *insertEdge(MeshPoint *a, MeshPoint *b, MeshPoint *o, const Vec &center) {
    edges.push_back(MeshEdge(a, b, o, center));
    MeshEdge *e = &edges.back();
    a->edges.push_back(e);
    b->edges.push_back(e);
    a->used = b->used = true;
    queue.push_back(e);
    return e;
  }

  static void link(MeshEdge *a, MeshEdge *b) { a->next = b; b->prev = a; }

  // The live front edge opposite to `edge`, if there is one.
  static MeshEdge *findReverseEdgeOnFront(const MeshEdge *edge) {
    using namespace ballpivoting;
    for (size_t i = 0; i < edge->a->edges.size(); i++) {
      MeshEdge *e = edge->a->edges[i];
      if (isLive(e) && e->a == edge->b && e->b == edge->a) return e;
    }
    return 0;
  }

  // Remove the coincident, oppositely oriented pair from the front and relink the loops:
  // the four cases of Fig. 7 (two-edge loop, consecutive edges, same loop, different loops).
  static void glue(MeshEdge *e1, MeshEdge *e2) {
    using namespace ballpivoting;
    if (e1->next == e2 && e2->next == e1) {
      // the two edges form a loop by themselves: nothing to relink
    } else if (e1->next == e2) {
      link(e1->prev, e2->next);
    } else if (e2->next == e1) {
      link(e2->prev, e1->next);
    } else {
      MeshEdge *p1 = e1->prev, *n1 = e1->next, *p2 = e2->prev, *n2 = e2->next;
      link(p1, n2);
      link(p2, n1);
    }
    e1->status = e2->status = INNER;
  }

  // Glue each of the freshly inserted edges against its opposite on the front, if any.
  static void glueOpposites(MeshEdge *const *es, int n) {
    using namespace ballpivoting;
    for (int i = 0; i < n; i++) {
      if (!isLive(es[i])) continue;
      if (MeshEdge *other = findReverseEdgeOnFront(es[i])) glue(es[i], other);
    }
  }

  // The seed triangle (a, b, c) as one loop of three edges, glued to the front where its
  // edges coincide with opposite front edges.
  void addSeed(const MeshFace &f, const Vec &center) {
    MeshEdge *es[3];
    es[0] = insertEdge(f.p[0], f.p[1], f.p[2], center);
    es[1] = insertEdge(f.p[1], f.p[2], f.p[0], center);
    es[2] = insertEdge(f.p[2], f.p[0], f.p[1], center);
    link(es[0], es[1]);
    link(es[1], es[2]);
    link(es[2], es[0]);
    glueOpposites(es, 3);
  }

  // The join operation (Fig. 6): the triangle (i, k, j) has been built by pivoting e(i, j);
  // e(i, j) leaves the front and e(i, k), e(k, j) take its place in the loop, then each new
  // edge is glued to its opposite on the front, if there is one (Fig. 7).
  void join(MeshEdge *e_ij, MeshPoint *k, const Vec &center) {
    using namespace ballpivoting;
    MeshEdge *prev = e_ij->prev, *next = e_ij->next;
    e_ij->status = INNER;
    MeshEdge *es[2];
    es[0] = insertEdge(e_ij->a, k, e_ij->b, center);
    es[1] = insertEdge(k, e_ij->b, e_ij->a, center);
    link(prev, es[0]);
    link(es[0], es[1]);
    link(es[1], next);
    glueOpposites(es, 2);
  }

  // ---- existing faces: their vertices are used, their border edges become boundary edges
  // of the front in the order the faces were created (which resumes the queue order of a
  // previous run), their inner edges between border vertices are recorded as closed.
  struct EdgeRec { int count; };
  typedef std::unordered_map<long long, EdgeRec> EdgeMap;

  long long edgeKey(int a, int b) const {
    return (long long)std::min(a, b) * (long long)points.size() + std::max(a, b);
  }

  void importExistingFaces() {
    using namespace ballpivoting;
    if (mesh.fn == 0) return;

    EdgeMap undirected;
    for (size_t i = 0; i < mesh.face.size(); i++) {
      const FaceType &f = mesh.face[i];
      if (f.IsD()) continue;
      for (int k = 0; k < 3; k++) {
        const int a = int(Index(mesh, f.cV(k))), b = int(Index(mesh, f.cV((k + 1) % 3)));
        points[a].used = true;
        undirected[edgeKey(a, b)].count++;
      }
    }

    // the border edges, wound like their triangle, in face order
    std::vector<MeshEdge *> border;
    std::vector<std::vector<MeshEdge *> > outgoing(points.size());
    for (size_t i = 0; i < mesh.face.size(); i++) {
      const FaceType &f = mesh.face[i];
      if (f.IsD()) continue;
      for (int k = 0; k < 3; k++) {
        const int a = int(Index(mesh, f.cV(k))), b = int(Index(mesh, f.cV((k + 1) % 3))), o = int(Index(mesh, f.cV((k + 2) % 3)));
        if (undirected[edgeKey(a, b)].count != 1) continue;
        edges.push_back(MeshEdge(&points[a], &points[b], &points[o], Vec(0, 0, 0), BOUNDARY));
        MeshEdge *e = &edges.back();
        border.push_back(e);
        outgoing[a].push_back(e);
      }
    }
    // link them into loops: after each edge a -> b comes an unclaimed border edge out of b
    for (size_t i = 0; i < border.size(); i++) {
      MeshEdge *e = border[i];
      std::vector<MeshEdge *> &out = outgoing[e->b->index];
      for (size_t j = 0; j < out.size(); j++)
        if (out[j]->prev == 0 && out[j] != e) { e->next = out[j]; out[j]->prev = e; break; }
    }
    // an edge that could not be linked both ways cannot be on the front; unlinking it may
    // orphan its neighbours, so repeat until stable
    for (bool changed = true; changed;) {
      changed = false;
      for (size_t i = 0; i < border.size(); i++) {
        MeshEdge *e = border[i];
        if (e->status == INNER) continue;
        if (e->prev == 0 || e->next == 0) {
          if (e->prev) e->prev->next = 0;
          if (e->next) e->next->prev = 0;
          e->prev = e->next = 0;
          e->status = INNER;
          changed = true;
        }
      }
    }
    std::vector<bool> isBorderVertex(points.size(), false);
    for (size_t i = 0; i < border.size(); i++) {
      MeshEdge *e = border[i];
      if (e->status == INNER) continue;
      e->a->edges.push_back(e);
      e->b->edges.push_back(e);
      isBorderVertex[e->a->index] = true;
      isBorderVertex[e->b->index] = true;
    }
    // the closed edges a pivot or a seed on the border could ask about
    for (size_t i = 0; i < mesh.face.size(); i++) {
      const FaceType &f = mesh.face[i];
      if (f.IsD()) continue;
      for (int k = 0; k < 3; k++) {
        const int a = int(Index(mesh, f.cV(k))), b = int(Index(mesh, f.cV((k + 1) % 3))), o = int(Index(mesh, f.cV((k + 2) % 3)));
        const int count = undirected[edgeKey(a, b)].count;
        if (count >= 2 && a > b) continue; // once per undirected edge
        if (!isBorderVertex[a] || !isBorderVertex[b]) continue;
        if (count == 1 && edgeIsLiveBorder(a, b)) continue;
        if (isClosed(&points[a], &points[b])) continue;
        edges.push_back(MeshEdge(&points[a], &points[b], &points[o], Vec(0, 0, 0), INNER));
        MeshEdge *e = &edges.back();
        e->a->edges.push_back(e);
        e->b->edges.push_back(e);
      }
    }
    // section 4.6: a boundary edge whose triangle admits an empty ball of this radius
    // resumes pivoting with that ball
    for (size_t i = 0; i < border.size(); i++) {
      MeshEdge *e = border[i];
      if (e->status != BOUNDARY) continue;
      Vec c;
      if (!computeBallCenter(MeshFace(e->a, e->b, e->opposite), radius, c)) continue;
      if (!grid->ballIsEmpty(c, radius, e->a, e->b, e->opposite)) continue;
      e->status = ACTIVE;
      e->center = c;
      queue.push_back(e);
      reactivated++;
    }
  }

  bool edgeIsLiveBorder(int a, int b) const {
    using namespace ballpivoting;
    const MeshPoint &p = points[a];
    for (size_t k = 0; k < p.edges.size(); k++) {
      const MeshEdge *e = p.edges[k];
      if (e->status == BOUNDARY && ((e->a == &p && e->b == &points[b]) || (e->b == &p && e->a == &points[b]))) return true;
    }
    return false;
  }

  // ---- seed search (section 4.2 of the paper). Cells are visited from a cursor that
  // persists between calls; a cell holding a used point is skipped (the paper's heuristic
  // against spawning small components next to the surface, fig. 4c), else one candidate is
  // tried, the point projecting furthest along the cell's average normal, paired with its
  // nearest neighbours: the first triangle facing along its vertex normals with an empty
  // ball on that side is the seed. A candidate that fails cannot succeed later (points only
  // become used), so the cursor never moves back.
  struct DistanceLess {
    Vec from;
    bool operator()(const CellEntry &a, const CellEntry &b) const { return (a.pos - from).SquaredNorm() < (b.pos - from).SquaredNorm(); }
  };

  bool findSeedTriangle(MeshFace &seed, Vec &ballCenter) {
    using namespace ballpivoting;
    const std::vector<Cell *> &cells = grid->cellsInOrder();
    for (; seedCursor < cells.size(); seedCursor++) {
      const Cell &cell = *cells[seedCursor];
      bool anyUsed = false;
      Vec navg(0, 0, 0), centroid(0, 0, 0);
      for (size_t i = 0; i < cell.size(); i++) {
        if (cell[i].p->used) { anyUsed = true; break; }
        navg += cell[i].p->normal;
        centroid += cell[i].pos;
      }
      if (anyUsed) continue;
      centroid /= double(cell.size());
      if (!cloudHasNormals) navg = unit(centroid - barycenter);
      MeshPoint *p1 = cell[0].p;
      double best = (cell[0].pos - centroid).dot(navg);
      for (size_t i = 1; i < cell.size(); i++) {
        const double d = (cell[i].pos - centroid).dot(navg);
        if (d > best) { best = d; p1 = cell[i].p; }
      }
      if (trySeed(p1, seed, ballCenter)) return true;
    }
    return false;
  }

  // Neighbour pairs (a, b) of s in order of distance, the interior vertices left out, the
  // candidate triangle oriented by the vertex normals, accepted when it keeps the mesh a
  // manifold and its ball, centred on the outward side, holds no other point.
  bool trySeed(MeshPoint *s, MeshFace &seed, Vec &ballCenter) {
    using namespace ballpivoting;
    grid->neighbors(s->pos, 2 * radius, neighborhood);
    seedPairs.clear();
    for (size_t i = 0; i < neighborhood.size(); i++)
      if (neighborhood[i].p != s && !isInterior(neighborhood[i].p)) seedPairs.push_back(neighborhood[i]);
    DistanceLess less = {s->pos};
    std::stable_sort(seedPairs.begin(), seedPairs.end(), less);
    if (seed_neighbors > 0 && seedPairs.size() > size_t(seed_neighbors)) seedPairs.resize(size_t(seed_neighbors));
    for (size_t ia = 0; ia < seedPairs.size(); ia++) {
      MeshPoint *a = seedPairs[ia].p;
      for (size_t ib = ia + 1; ib < seedPairs.size(); ib++) {
        MeshPoint *b = seedPairs[ib].p;
        MeshFace f(s, a, b);
        const Vec n = f.normal();
        if (cloudHasNormals) {
          // the seed must face along the normals of all three of its vertices (paper,
          // section 4.2): a vertex without a normal cannot fix a seed's winding
          if (agrees(n, f)) {}
          else if (agrees(-n, f)) std::swap(f.p[1], f.p[2]);
          else continue;
        } else {
          if (n.dot((s->pos + a->pos + b->pos) / 3.0 - barycenter) < 0) std::swap(f.p[1], f.p[2]);
        }
        if (!canAddSeed(f)) continue;
        Vec c;
        if (!computeBallCenter(f, radius, c)) continue;
        if (!grid->ballIsEmpty(c, radius, f.p[0], f.p[1], f.p[2])) continue;
        seed = f;
        ballCenter = c;
        return true;
      }
    }
    return false;
  }

  static bool agrees(const Vec &n, const MeshFace &f) {
    return n.dot(f.p[0]->normal) > 0 && n.dot(f.p[1]->normal) > 0 && n.dot(f.p[2]->normal) > 0;
  }

  // Manifoldness test for a seed, all of whose vertices may already be in the mesh.
  static bool canAddSeed(const MeshFace &f) {
    using namespace ballpivoting;
    for (int i = 0; i < 3; i++) {
      if (isInterior(f.p[i])) return false;
      const MeshPoint *p = f.p[i], *q = f.p[(i + 1) % 3];
      if (hasEdge(p, q) || isClosed(p, q)) return false;
    }
    return true;
  }

  // The tests the paper applies to the point k the ball lands on when pivoting e = (a, b)
  // (fig. 5, line 3, and the "edge orientation checks" it mentions): the triangle (a, k, b)
  // must not face against the normal of k (a zero dot product passes: the front edge
  // already fixes the winding), k must be unused or on the front, and the mesh must stay a
  // manifold: neither new half-edge a -> k, k -> b may already be on the front with the
  // same orientation, nor either undirected edge closed. Optionally (max_angle) the
  // dihedral angle with the triangle of e is bounded.
  bool canAddTriangle(const MeshEdge *e, const MeshPoint *k) const {
    using namespace ballpivoting;
    const Vec normal = (k->pos - e->a->pos) ^ (e->b->pos - e->a->pos);
    if (normal.dot(k->normal) < 0) return false;
    if (!(!k->used || onFront(k))) return false;
    if (hasEdge(e->a, k) || hasEdge(k, e->b) || isClosed(e->a, k) || isClosed(k, e->b)) return false;
    if (max_angle > -1) {
      const Vec old = (e->b->pos - e->a->pos) ^ (e->opposite->pos - e->a->pos);
      if (unit(old).dot(unit(normal)) < max_angle) return false;
    }
    return true;
  }

  // Preference among points hit simultaneously: 0 if the triangle would be rejected, else
  // 1 plus the number of its new edges that glue to an existing front edge. The first pivot
  // into a cospherical polygon fixes a diagonal; later pivots can only respect it.
  int tieScore(const MeshEdge *e, const MeshPoint *k) const {
    using namespace ballpivoting;
    if (!canAddTriangle(e, k)) return 0;
    return 1 + int(hasEdge(k, e->a)) + int(hasEdge(e->b, k));
  }

  // Roll the ball around edge e from its stored centre and return the point it touches
  // first with the centre at that moment. Every point within reach of the ball is a
  // candidate, the opposite vertex of e included: if the ball comes back to it before
  // touching anything else there is no triangle to build. Because the initial ball is
  // empty, the ball at the first contact is empty too.
  bool ballPivot(const MeshEdge *e, MeshPoint *&hit, Vec &center) {
    using namespace ballpivoting;
    PivotFrame fr;
    if (!pivotFrame(e, fr)) return false;
    // any point the ball can touch lies within r + radius of the midpoint
    grid->neighbors(fr.m, std::min(2 * radius, (fr.r + radius) * (1 + 1e-8)), neighborhood);
    const double r2 = fr.r * fr.r;
    const double tieSin = std::sin(TIE_TOLERANCE);

    MeshPoint *best = 0;
    Contact bestContact;
    int nTies = 1;
    for (size_t i = 0; i < neighborhood.size(); i++) {
      MeshPoint *p = neighborhood[i].p;
      if (p == e->a || p == e->b) continue;
      Contact c;
      if (!pivotContact(fr, neighborhood[i].pos, radius, c)) continue;
      if (!best) {
        best = p;
        bestContact = c;
      } else if (angleTie(c, bestContact, r2, tieSin)) {
        nTies++;
        if (angleLess(c, bestContact)) bestContact = c;
      } else if (angleLess(c, bestContact)) {
        best = p;
        bestContact = c;
        nTies = 1;
      }
    }
    if (!best) return false;
    if (best == e->opposite && nTies == 1) return false; // the ball came back to the opposite vertex first
    if (nTies == 1) {
      hit = best;
      center = centerAt(fr, bestContact);
      return true;
    }

    // Simultaneous hits: pick deterministically. The opposite vertex is never chosen; a
    // point hit at the same angle gives a valid triangle with it on the ball's surface.
    int bestScore = -1;
    for (size_t i = 0; i < neighborhood.size(); i++) {
      MeshPoint *p = neighborhood[i].p;
      if (p == e->a || p == e->b || p == e->opposite) continue;
      Contact c;
      if (!pivotContact(fr, neighborhood[i].pos, radius, c)) continue;
      if (!(angleTie(c, bestContact, r2, tieSin) || angleLess(c, bestContact))) continue;
      const int score = tieScore(e, p);
      if (score > bestScore || (score == bestScore && p->index < best->index)) {
        bestScore = score;
        best = p;
        bestContact = c;
      }
    }
    hit = best;
    center = centerAt(fr, bestContact);
    return true;
  }
};

} // namespace tri
} // namespace vcg
#endif
