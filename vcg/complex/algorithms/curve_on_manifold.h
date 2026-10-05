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
#ifndef __VCGLIB_CURVE_ON_SURF_H
#define __VCGLIB_CURVE_ON_SURF_H

#include<vcg/complex/complex.h>
#include<vcg/simplex/face/topology.h>
#include<vcg/complex/algorithms/update/topology.h>
#include<vcg/complex/algorithms/update/color.h>
#include<vcg/complex/algorithms/update/normal.h>
#include<vcg/complex/algorithms/update/quality.h>
#include<vcg/complex/algorithms/clean.h>
#include<vcg/complex/append.h>
#include<vcg/complex/algorithms/mesh_assert.h>
#include<vcg/complex/algorithms/update/bounding.h>
#include<vcg/complex/algorithms/vertex_interpolation.h>
#include<vcg/complex/algorithms/point_sampling.h>
#include <vcg/space/index/grid_static_ptr.h>
#include <vcg/math/histogram.h>
#include<vcg/space/distance3.h>
#include <wrap/callback.h>
#include <vcg/space/planar_polygon_tessellation.h>
#include <array>
#include <deque>
#include <map>
#include <set>

namespace vcg {
namespace tri {
/// \ingroup trimesh
/// \brief A class for managing curves on a 2-manifold (Curve on Manifold - CoM).
/**
 * This class is used to project, simplify, smooth, and snap polylines (represented as edge meshes) 
 * over a triangulated surface (the "base mesh").
 * 
 * \par Overview
 * The CoM class provides tools to:
 * - Project polylines onto a surface
 * - Snap polyline vertices to mesh vertices or edges
 * - Refine polylines to follow surface features
 * - Simplify polylines while maintaining surface fidelity
 * - Split the base mesh along polylines for mesh cutting operations (in CoMEmbed)
 * 
 * \par Terminology
 * - **Base mesh**: The triangulated surface mesh (stored in `base`)
 * - **Polyline/Curve**: An edge mesh passed to the various methods (typically named `poly` in parameters)
 * - **Snapping**: The process of aligning polyline vertices to mesh vertices or edges using barycentric thresholds
 *
 * \par Curves, strands and control points
 * A curve on the surface is a graph of **control points**, joined by **strands**. A strand
 * stands for a geodesic between its two control points: its other vertices are only a
 * discretization of it, which smoothing moves towards the geodesic and refinement and
 * simplification add or remove. Control points are fixed: they are projected onto the
 * surface and never smoothed, snapped, simplified or merged away. They are the selected
 * vertices of the polyline; SetControlPoints() chooses them: every vertex, the ends and
 * junctions (the default), or the selection as given. Since a geodesic on a triangle mesh
 * is straight inside each face, RefineCurveByBaseMesh() leaves a strand with only its
 * control points and its exact edge crossings.
 * 
 * \par Requirements
 * The base mesh should be:
 * - 2-manifold
 * - Have reasonable triangle quality
 * - Have up-to-date topology (FaceFace adjacency)
 * - Have bounding box correctly set
 * 
 * \par Usage Pattern
 * 1. Initialize the class with a base mesh: `CoM<MeshType> com(baseMesh);`
 * 2. Call `Init()` to build the spatial acceleration structures
 * 3. Pass polylines (as edge meshes) to various methods for processing
 * 4. Adjust parameters via the `par` member for fine control
 * 
 * \par Implementation Notes
 * - The class uses barycentric coordinates to determine if a point of the polyline should snap to vertices or edges
 * - All spatial queries use a uniform grid for acceleration
 * - Many operations are iterative and may require multiple passes
 * 
 * \par Errors
 * Violated preconditions throw vcg::MissingPreconditionException instead of asserting:
 * a base mesh with no faces, non-manifold edges or zero-area faces (Init()), an empty
 * curve (MoveAndProject()), a curve point farther than Param::gridBailout from the
 * surface (every closest-face query), and a curve that CoMEmbed::SplitMeshWithPolyline() cannot
 * embed. After an exception the base mesh and the curve may be partially modified, so
 * a caller that needs them intact should work on copies.
 *
 * \par Attributes
 * The methods that change the base mesh keep its per-vertex color, quality, normal and
 * texture coordinates and its per-face attributes: the vertices and faces created by
 * splitting interpolate or inherit them. Progress and diagnostic messages go to
 * Param::cb, when set. The methods that change the base mesh are in CoMEmbed.
 *
 * \note There is some naming inconsistency: methods use both "Curve" and "Polyline" 
 *       interchangeably to refer to the edge mesh being processed.
 * 
 */

template <class MeshType> class CoMEmbed;

template <class MeshType>
class CoM
{
public:
  typedef typename MeshType::ScalarType     ScalarType;
  typedef typename MeshType::CoordType      CoordType;
  typedef typename MeshType::VertexType     VertexType;
  typedef typename MeshType::VertexPointer  VertexPointer;
  typedef typename MeshType::VertexIterator VertexIterator;
  typedef typename MeshType::EdgeIterator   EdgeIterator;
  typedef typename MeshType::EdgeType       EdgeType;
  typedef typename MeshType::FaceType       FaceType;
  typedef typename MeshType::FacePointer    FacePointer;
  typedef typename MeshType::FaceIterator   FaceIterator;
  typedef Box3<ScalarType>                  Box3Type;
  typedef Segment3<ScalarType>              Segment3Type;  
  typedef typename vcg::GridStaticPtr<FaceType, ScalarType> MeshGrid;  
  typedef typename vcg::GridStaticPtr<EdgeType, ScalarType> EdgeGrid;
  typedef typename face::Pos<FaceType> PosType;
  typedef typename tri::UpdateTopology<MeshType>::PEdge PEdge;
  
  /**
   * \brief Parameter class controlling the behavior of CoM algorithms
   * 
   * This class contains all the thresholds and tolerances used by the various
   * curve-on-manifold operations. Default values are computed relative to the
   * bounding box diagonal of the base mesh.
   */
  class Param 
  {
  public:
    
    ScalarType surfDistThr;        ///< Max distance between surface and curve; used in simplify and refine
    ScalarType minRefEdgeLen;      ///< Minimal admitted edge length (used in refine: never make edges shorter than this) 
    ScalarType maxSimpEdgeLen;     ///< Maximal admitted edge length (used in simplify: never make edges longer than this) 
    ScalarType maxMoveDelta;       ///< The maximum movement admitted during MoveAndProject (before projection) 
    ScalarType maxSnapThr;         ///< The maximum distance a snap may move a polyline vertex (onto a mesh vertex or edge)
    ScalarType gridBailout;        ///< The maximum distance bailout used in grid-based spatial queries
    ScalarType barycentricSnapThr; ///< Threshold for snapping barycentric coords to 0 or 1 (controls vertex/edge snapping)
    vcg::CallBackPos *cb = nullptr; ///< Receives progress and diagnostic messages (see wrap/callback.h)
    
    /// Constructor with default parameter initialization based on mesh size
    Param(MeshType &m) { SetDefault(m);}
    
    /// Set all parameters to reasonable defaults based on the mesh bounding box
    void SetDefault(MeshType &m)
    {
      surfDistThr        = m.bbox.Diag()/1000.0;
      minRefEdgeLen      = m.bbox.Diag()/16000.0;
      maxSimpEdgeLen     = m.bbox.Diag()/100.0;
      maxMoveDelta       = m.bbox.Diag()/100.0;
      maxSnapThr         = m.bbox.Diag()/10000.0;
      gridBailout        = m.bbox.Diag()/20.0;
      barycentricSnapThr = 0.05;
    }
    
    /// Print current parameter values to stdout
    void Dump() const
    {
      printf("surfDistThr    = %6.4f\n",surfDistThr   );
      printf("minRefEdgeLen  = %6.4f\n",minRefEdgeLen    );
      printf("maxSimpEdgeLen = %6.4f\n",maxSimpEdgeLen    );
      printf("maxMoveDelta   = %6.4f\n",maxMoveDelta);
    }
  };
  
  
  
  // ============================================================================
  // Data Members
  // ============================================================================
  
  MeshType &base;       ///< Reference to the base triangulated surface mesh
  MeshGrid uniformGrid; ///< Spatial acceleration structure for closest point queries
  Param par;            ///< Parameters controlling algorithm behavior
  
  /// Constructor: initializes the CoM with a base mesh. It updates the bounding box of
  /// the base mesh first, because every default in Param is a fraction of its diagonal:
  /// a stale box would silently give tolerances of the wrong scale.
  CoM(MeshType &_m) :base(_m),par(BoxUpdated(_m)){}

private:
  static MeshType &BoxUpdated(MeshType &m) { tri::UpdateBounding<MeshType>::Box(m); return m; }

  template <class> friend class CoMEmbed;

  int progress = 0; ///< last progress reported through Param::cb, in [0,100]

  /// Report progress \a pos with a message through Param::cb.
  template <class... Args>
  void Progress(int pos, const char *fmt, Args... args)
  {
    progress = pos;
    Log(fmt, args...);
  }

  /// Send a diagnostic message through Param::cb, at the progress last reported: a fixed
  /// position would move a caller's progress bar back.
  template <class... Args>
  void Log(const char *fmt, Args... args)
  {
    if (!par.cb) return;
    char msg[256];
    snprintf(msg, sizeof(msg), fmt, args...);
    par.cb(progress, msg);
  }

  /// Every closest-face query goes through here. The grid returns no face beyond
  /// Param::gridBailout, and all callers need one, so a null result is an error to
  /// report, not a pointer to dereference.
  FaceType *Closest(const CoordType &p, CoordType &closestP, ScalarType &closestDist)
  {
    FaceType *f = vcg::tri::GetClosestFaceBase(base, uniformGrid, p, par.gridBailout, closestDist, closestP);
    if (f == nullptr) {
      char msg[256];
      snprintf(msg, sizeof(msg), "CoM: no face of the base mesh within %g (Param::gridBailout) of the point (%g, %g, %g); the curve is too far from the surface.",
               double(par.gridBailout), double(p[0]), double(p[1]), double(p[2]));
      throw vcg::MissingPreconditionException(msg);
    }
    return f;
  }

public:
 
  // ============================================================================
  // Spatial Query Methods
  // ============================================================================
  
  /**
   * \brief Get the closest face to a query point
   * \param p The query point
   * \return Pointer to the closest face
   * \throws vcg::MissingPreconditionException if no face is within Param::gridBailout,
   *         as all the closest-face queries below
   */
  FaceType *GetClosestFace(const CoordType &p)
  {
    ScalarType closestDist;
    CoordType closestP;
    return Closest(p, closestP, closestDist);
  }
  
  /**
   * \brief Get the closest face and its barycentric coordinates
   * \param p The query point
   * \param ip Output: barycentric coordinates of the closest point on the returned face
   * \return Pointer to the closest face
   */
  FaceType *GetClosestFaceIP(const CoordType &p, CoordType &ip)
    {
      ScalarType closestDist;
      CoordType closestP;
      FaceType *f = Closest(p, closestP, closestDist);
      InterpolationParameters(*f, f->N(), closestP, ip);
      return f;
    }

  /**
   * \brief Get the closest face, barycentric coordinates, and normal
   * \param p The query point
   * \param ip Output: barycentric coordinates of the closest point on the returned face
   * \param in Output: normal at the closest point
   * \return Pointer to the closest face
   */
  FaceType *GetClosestFaceIP(const CoordType &p, CoordType &ip, CoordType &in)
    {
      FaceType *f = GetClosestFaceIP(p, ip);
      in = f->V(0)->cN()*ip[0] + f->V(1)->cN()*ip[1] + f->V(2)->cN()*ip[2];
      return f;
    }

  /**
   * \brief Get the closest face and the closest point on it
   * \param p The query point
   * \param closestP Output: the 3D coordinates of the closest point on the surface
   * \return Pointer to the closest face
   */
  FaceType *GetClosestFacePoint(const CoordType &p, CoordType &closestP)
  {
    ScalarType closestDist;
    return Closest(p, closestP, closestDist);
  }
  
  // ============================================================================
  // Barycentric Coordinate and Snapping Methods
  // ============================================================================
  
  /**
   * \brief Check if a barycentric coordinate is snapped to an edge
   * \param ip The barycentric coordinate (must be snapped)
   * \param ei Output: index (0,1,2) of the edge opposite to the zero coordinate
   * \return true if snapped to an edge (exactly one coordinate is 0, other two are positive)
   * 
   * Edge snapping means the point lies on one of the three edges of the triangle.
   * Edge i is the edge from vertex i to vertex (i+1)%3, opposite to vertex (i+2)%3.
   */
  bool IsSnappedEdge(CoordType &ip, int &ei)
  {
    for(int i=0;i<3;++i)
      if(ip[i]>0.0 && ip[(i+1)%3]>0.0 && ip[(i+2)%3]==0.0 ) {
        ei=i;
        return true; 
      }
    ei=-1;
    return false;
  }

  /**
   * \brief Check if a barycentric coordinate is snapped to a vertex
   * \param ip The barycentric coordinate (must be snapped)
   * \param vi Output: index (0,1,2) of the vertex with coordinate == 1.0
   * \return true if snapped to a vertex (one coordinate is 1.0, others are 0.0)
   */
  bool IsSnappedVertex(CoordType &ip, int &vi)
  {
    for(int i=0;i<3;++i)
      if(ip[i]==1.0 && ip[(i+1)%3]==0.0 && ip[(i+2)%3]==0.0 ) {
        vi=i;
        return true; 
      }
    vi=-1;
    return false;
  }

  /**
   * \brief Find the vertex pointer for a vertex-snapped barycentric coordinate
   * \param fp The face containing the point
   * \param ip The barycentric coordinate (should be snapped to a vertex)
   * \return Pointer to the snapped vertex, or nullptr if not vertex-snapped
   */
  VertexPointer FindVertexSnap(FacePointer fp, CoordType &ip)
  {
    for(int i=0;i<3;++i)
      if(ip[i]==1.0 && ip[(i+1)%3]==0.0 && ip[(i+2)%3]==0.0 ) return fp->V(i);
    return 0;
  }
  

  
  // ============================================================================
  // Utility Functions
  // ============================================================================
  
  /**
   * @brief SnapPolyline snaps the vertexes of a polyline onto the base mesh
   * @param poly
   * @param newVertVec the vector of the indexes of the snapped vertices
   * @return true if it has modified the polyline
   * 
   * Polyline vertices can be snapped either on vertexes or on edges. 
   * Usually the only points that we should allow to not be snapped are the endpoints and non manifold points.
   * Vertexes are colored according to their snapping state 
   * 
   */  
    
  bool SnapPolyline(MeshType &poly)
  {
    tri::Allocator<MeshType>::CompactEveryVector(poly);     
    tri::UpdateTopology<MeshType>::VertexEdge(poly);
    // Where each vertex would snap.
    std::vector<FaceType *> face(poly.vert.size());
    std::vector<CoordType> ip(poly.vert.size()), raw(poly.vert.size());
    std::map<VertexType *, std::vector<size_t>> onVertex;  // mesh vertex -> polyline vertices
    for (size_t i = 0; i < poly.vert.size(); ++i)
    {
      face[i] = GetClosestFaceIP(poly.vert[i].cP(), raw[i]);
      ip[i] = raw[i];
      // A control point is not snapped; one already on a mesh vertex occupies it.
      if (poly.vert[i].IsS()) RoundingSnap(ip[i], *face[i]);
      else BarycentricSnap(ip[i], *face[i]);
      if (VertexPointer v = FindVertexSnap(face[i], ip[i])) onVertex[v].push_back(i);
    }
    // A mesh vertex takes one polyline vertex: the nearest, with the ones joined to it by
    // polyline edges (they will collapse onto it). The others stay off it, on an edge if
    // they may snap there, so two strands passing by the same vertex are not merged.
    int deniedCnt = 0;
    for (auto &ov : onVertex)
    {
      if (ov.second.size() < 2) continue;
      std::sort(ov.second.begin(), ov.second.end(), [&](size_t a, size_t b) {
        if (poly.vert[a].IsS() != poly.vert[b].IsS()) return poly.vert[a].IsS();  // a control point first
        return Distance(poly.vert[a].cP(), ov.first->cP()) < Distance(poly.vert[b].cP(), ov.first->cP()); });
      std::set<size_t> kept{ov.second[0]};
      for (bool grew = true; grew; )  // the nearest and what is joined to it along the polyline
      {
        grew = false;
        for (size_t c : ov.second)
        {
          if (kept.count(c)) continue;
          std::vector<VertexPointer> star;
          edge::VVStarVE(&poly.vert[c], star);
          for (VertexPointer w : star)
            if (kept.count(tri::Index(poly, w))) { kept.insert(c); grew = true; break; }
        }
      }
      for (size_t c : ov.second)
      {
        if (kept.count(c) || poly.vert[c].IsS()) continue;
        // Only onto the nearest edge, if that snap alone is allowed.
        ++deniedCnt;
        CoordType &q = ip[c];
        q = raw[c];
        const int k = (q[0] <= q[1] && q[0] <= q[2]) ? 0 : (q[1] <= q[2] ? 1 : 2);
        const FaceType &f = *face[c];
        const ScalarType len = Distance(f.cP((k + 1) % 3), f.cP((k + 2) % 3));
        if (q[k] <= par.barycentricSnapThr && len > 0 && q[k] * (DoubleArea(f) / len) <= par.maxSnapThr)
        {
          q[k] = 0;
          q[(k + 1) % 3] /= q[(k + 1) % 3] + q[(k + 2) % 3];
          q[(k + 2) % 3] = 1 - q[(k + 1) % 3];
        }
      }
    }
    int vertSnapCnt=0, edgeSnapCnt=0;
    for (size_t i = 0; i < poly.vert.size(); ++i)
    {
      if (poly.vert[i].IsS()) continue;  // control points stay where they are
      int zeros = 0;
      for (int k = 0; k < 3; ++k) zeros += ip[i][k] == 0;
      if (zeros == 0) continue;
      const FaceType &f = *face[i];
      poly.vert[i].P() = f.cP(0)*ip[i][0] + f.cP(1)*ip[i][1] + f.cP(2)*ip[i][2];
      if (zeros == 2) vertSnapCnt++; else edgeSnapCnt++;
    }
    const int dupCnt = CollapseNullEdges(poly);
    Log("SnapPolyline %i vertices: snapped %i onto vertices and %i onto edges, %i kept off an occupied vertex, %i null edges collapsed",
        poly.vn, vertSnapCnt, edgeSnapCnt, deniedCnt, dupCnt);
    return vertSnapCnt==0 && edgeSnapCnt==0 && dupCnt==0;
  }
  
  /// Which polyline vertices SetControlPoints() makes control points.
  enum ControlPoints
  {
    AllVertices,   ///< every vertex: the polyline is taken as it is, only projected and refined
    EndsAndNodes,  ///< the ends (degree 1) and the junctions (degree > 2): each strand is free
    Selected       ///< the vertices already selected, left as they are
  };

  /**
   * \brief Choose the control points of a polyline (stored as its vertex selection).
   *
   * Call it once on the input curve, before smoothing or refining it. With EndsAndNodes
   * the curve keeps its topology and its ends, and every strand between them is free to
   * become a geodesic; a helix on a cylinder with fixed ends stays a helix, since the
   * strand cannot unwind around the cylinder.
   */
  void SetControlPoints(MeshType &poly, ControlPoints mode = EndsAndNodes)
  {
    if (mode == Selected) return;
    tri::UpdateTopology<MeshType>::VertexEdge(poly);
    ForEachVertex(poly, [&](VertexType &v) {
      const int deg = edge::VEDegree<EdgeType>(&v);
      if (mode == AllVertices || deg != 2) v.SetS(); else v.ClearS();
    });
  }

  /**
   * \brief Make the points where strands meet into junctions (opt-in).
   * \return the number of junctions made
   *
   * After CoMEmbed::SplitMeshWithPolyline() two strands that cross or touch each have a
   * vertex at the same position. By default they stay separate strands; this merges every
   * set of coinciding vertices into one and merges the edges two strands share, so the
   * meeting points become part of the curve's topology; where the merged curve branches
   * the vertex is a junction, and a control point.
   */
  int ConnectCrossings(MeshType &poly)
  {
    std::map<CoordType, int> count;
    for (const VertexType &v : poly.vert) if (!v.IsD()) ++count[v.cP()];
    tri::Clean<MeshType>::RemoveDuplicateVertex(poly);
    tri::Clean<MeshType>::RemoveDuplicateEdge(poly);
    tri::Allocator<MeshType>::CompactEveryVector(poly);
    // Where the merged strands now branch it is a junction; where they only ran together
    // (an overlapping stretch) the vertex is an ordinary one.
    tri::UpdateTopology<MeshType>::VertexEdge(poly);
    int junctions = 0;
    for (VertexType &v : poly.vert)
      if (count[v.cP()] > 1 && edge::VEDegree<EdgeType>(&v) > 2) { v.SetS(); ++junctions; }
    return junctions;
  }

  
   void SelectUniformlyDistributed(MeshType &poly, int k)
   {
     tri::TrivialPointerSampler<MeshType> tps;
     ScalarType samplingRadius = tri::Stat<MeshType>::ComputeEdgeLengthSum(poly)/ScalarType(k);
     tri::SurfaceSampling<MeshType, typename tri::TrivialPointerSampler<MeshType> >::EdgeMeshUniform(poly,tps,samplingRadius);     
     for(int i=0;i<tps.sampleVec.size();++i)
       tps.sampleVec[i]->SetS();
   }
   
   
    
  /*
   * Make an edge mesh 1-manifold by splitting all the
   * vertexes that have more than two incident edges
   * 
   * It performs the split in three steps. 
   * - First it collects and counts the vertices to be splitten. 
   * - Then it adds the vertices to the mesh and 
   * - lastly it updates the poly with the newly added vertices. 
   * 
   * singSplitFlag allows to ubersplit each singularity in a number of vertex of the same order of its degree. 
   * This is not really necessary but helps the management of sharp turns in the poly mesh.
   * \todo add corner detection and split.
   */
  
  void DecomposeNonManifoldPolyline(MeshType &poly, bool singSplitFlag = true)
  {
    tri::Allocator<MeshType>::CompactEveryVector(poly);
    std::vector<int> degreeVec(poly.vn, 0);
    tri::UpdateTopology<MeshType>::VertexEdge(poly);
    int neededVert=0;
    int delta;
    if(singSplitFlag) delta = 1;
                 else delta = 2;
      
    for(VertexIterator vi=poly.vert.begin(); vi!=poly.vert.end();++vi)
    {
      std::vector<EdgeType *> starVec;
      edge::VEStarVE(&*vi,starVec);
      degreeVec[tri::Index(poly, *vi)] = starVec.size();
      if(starVec.size()>2)
        neededVert += starVec.size()-delta;
    }
    Log("DecomposeNonManifold Adding %i vert to a polyline of %i vert",neededVert,poly.vn);
    VertexIterator firstVi = tri::Allocator<MeshType>::AddVertices(poly,neededVert);
    
    for(size_t i=0;i<degreeVec.size();++i)
    {
      if(degreeVec[i]>2)
      {
        std::vector<EdgeType *> edgeStarVec;
        edge::VEStarVE(&(poly.vert[i]),edgeStarVec);
        assert(edgeStarVec.size() == degreeVec[i]);
        for(size_t j=delta;j<edgeStarVec.size();++j)
        {
          EdgeType *ep = edgeStarVec[j];
          int ind; // index of the vertex to be changed
          if(tri::Index(poly,ep->V(0)) == i) ind = 0;
              else ind = 1;
  
          ep->V(ind) = &*firstVi;
          ep->V(ind)->P() = poly.vert[i].P();
          ep->V(ind)->N() = poly.vert[i].N();
          ++firstVi;
        }
      }
    }
    assert(firstVi == poly.vert.end());
  }
  
  // ============================================================================
  // Initialization
  // ============================================================================
  
  /**
   * \brief Initialize the CoM data structures for processing
   * 
   * This must be called after construction and whenever the base mesh is modified.
   * It performs:
   * - Face normal computation
   * - Face-Face topology update
   * - Spatial acceleration structure (uniform grid) construction
   * 
   * \note Call this before using any polyline processing methods.
   * \warning If the base mesh topology changes, call Init() again.
   */
  void Init()
  {
    if (base.fn == 0)
      throw vcg::MissingPreconditionException("CoM: the base mesh has no faces.");
    UpdateNormal<MeshType>::PerFaceNormalized(base);
    UpdateTopology<MeshType>::FaceFace(base);    
    MeshAssert<MeshType>::FFTwoManifoldEdge(base);
    // A zero-area face has no barycentric coordinates: the closest-point query either
    // misses it or returns uninitialized ones, which the snapping then trusts.
    MeshAssert<MeshType>::NoZeroAreaFace(base);
    // Construction of the uniform grid
    uniformGrid.Set(base.face.begin(), base.face.end());    
  }
  
  // ============================================================================
  // Simplification Methods
  // ============================================================================
  
  /**
   * \brief Remove duplicate/zero-length edges from a polyline
   * \param poly The polyline to simplify
   * 
   * Removes vertices that have collapsed to the same position.
   */
  void SimplifyNullEdges(MeshType &poly)
  {
      int cnt=CollapseNullEdges(poly);
      if(cnt)
          Log("SimplifyNullEdges: Collapsed %i zero-length edges",cnt);
  }

  /**
   * \brief Merge the two ends of every zero-length polyline edge.
   * \return the number of edges collapsed
   *
   * Vertices are merged only along the polyline, never because they happen to coincide:
   * two strands passing through the same point stay two strands, and are not joined into a
   * junction that smoothing would then treat as one. (Clean::RemoveDuplicateVertex, used
   * here before, merged them.) Compacts \a poly and updates its VE adjacency.
   */
  int CollapseNullEdges(MeshType &poly)
  {
    tri::UpdateTopology<MeshType>::VertexEdge(poly);
    int cnt = 0;
    for (bool changed = true; changed; )
    {
      changed = false;
      for (size_t i = 0; i < poly.edge.size(); ++i)
      {
        EdgeType &e = poly.edge[i];
        if (e.IsD() || e.V(0)->cP() != e.V(1)->cP()) continue;
        if (e.V(0) == e.V(1)) { edge::VEDetach(e); tri::Allocator<MeshType>::DeleteEdge(poly, e); }
        else edge::VEEdgeCollapseToVertex(poly, &e, e.V(1)->IsS() && !e.V(0)->IsS() ? 1 : 0); // keep a control point
        ++cnt;
        changed = true;
      }
    }
    tri::Allocator<MeshType>::CompactEveryVector(poly);
    tri::UpdateTopology<MeshType>::VertexEdge(poly);
    return cnt;
  }
  void Simplify(MeshType &poly)
  {
    int startEn = poly.en;
    Distribution<ScalarType> hist;
    for(int i =0; i<poly.en;++i) 
      hist.Add(edge::Length(poly.edge[i]));
        
    UpdateTopology<MeshType>::VertexEdge(poly);
    
    for(int i =0; i<poly.vn;++i)
    {
      std::vector<VertexPointer> starVecVp;
      edge::VVStarVE(&(poly.vert[i]),starVecVp);      
      if ((starVecVp.size()==2) && (!poly.vert[i].IsS()))
      {
        ScalarType newSegLen = Distance(starVecVp[0]->P(), starVecVp[1]->P());
        Segment3Type seg(starVecVp[0]->P(),starVecVp[1]->P());
        ScalarType segDist;
        CoordType closestPSeg;
        SegmentPointDistance(seg,poly.vert[i].cP(),closestPSeg,segDist);
        CoordType fp,fn;
        ScalarType maxSurfDist = MaxSegDist(starVecVp[0], starVecVp[1],fp,fn);
        
        if((maxSurfDist < par.surfDistThr) && (newSegLen < par.maxSimpEdgeLen) )
        {
          edge::VEEdgeCollapse(poly,&(poly.vert[i]));          
        }
      }
    }
    tri::UpdateTopology<MeshType>::TestVertexEdge(poly);
    tri::Allocator<MeshType>::CompactEveryVector(poly);
    tri::UpdateTopology<MeshType>::TestVertexEdge(poly);
//    printf("Simplify %5i -> %5i (total len %5.2f)\n",startEn,poly.en,hist.Sum());
  }
  
  void EvaluateHausdorffDistance(MeshType &poly, Distribution<ScalarType> &dist)
  {
    dist.Clear();
    tri::UpdateTopology<MeshType>::VertexEdge(poly);
    tri::UpdateQuality<MeshType>::VertexConstant(poly,0);
    for(int i =0; i<poly.edge.size();++i)
    {      
      CoordType farthestP, farthestN;      
      ScalarType maxDist = MaxSegDist(poly.edge[i].V(0),poly.edge[i].V(1), farthestP, farthestN, &dist);      
      poly.edge[i].V(0)->Q()+= maxDist;
      poly.edge[i].V(1)->Q()+= maxDist;
    }
    for(int i=0;i<poly.vn;++i)
    {
      ScalarType deg = edge::VEDegree<EdgeType>(&poly.vert[i]);
      poly.vert[i].Q()/=deg;
    }
    tri::UpdateColor<MeshType>::PerVertexQualityRamp(poly,0,dist.Max());    
  }
  

  /**
   * \brief Snap the barycentric coordinates of a point of face \a f onto a vertex or an edge.
   * \param ip Input/Output: barycentric coordinates in \a f (summing to 1)
   * \param f  The face they refer to
   * \return true if the point is now on a vertex or an edge (at least one coordinate is 0)
   *
   * A coordinate is set to 0 only if it is within Param::barycentricSnapThr of 0 **and**
   * the snap moves the point by at most Param::maxSnapThr. The barycentric bound alone is
   * only meaningful on well-shaped, uniformly sized triangles: on a large or skinny one it
   * can move a point far, which the distance bound prevents; on tiny triangles the
   * barycentric bound keeps a point from being snapped across them. Two coordinates at 0
   * put the point on a vertex. The remaining coordinates are renormalized so they sum to
   * exactly 1. Setting either threshold to 0 disables snapping.
   *
   * \sa IsSnappedVertex, IsSnappedEdge
   */
  bool BarycentricSnap(CoordType &ip, const FaceType &f)
  {
    return BarycentricSnap(ip, f, par.barycentricSnapThr, par.maxSnapThr);
  }

  /// BarycentricSnap() with explicit thresholds, e.g. rounding-level ones.
  static bool BarycentricSnap(CoordType &ip, const FaceType &f, ScalarType barThr, ScalarType distThr)
  {
    // Setting coordinate i to 0 moves the point onto the opposite edge, by ip[i] times the
    // height of the triangle over that edge.
    const ScalarType area2 = DoubleArea(f);
    bool zero[3];
    for (int i = 0; i < 3; ++i)
    {
      const ScalarType len = Distance(f.cP((i + 1) % 3), f.cP((i + 2) % 3));
      zero[i] = ip[i] <= barThr && len > 0 && ip[i] * (area2 / len) <= distThr;
    }
    if (zero[0] && zero[1] && zero[2])  // cannot leave the triangle: keep the largest
      zero[ip[0] >= ip[1] && ip[0] >= ip[2] ? 0 : (ip[1] >= ip[2] ? 1 : 2)] = false;
    ScalarType sum = 0;
    for (int i = 0; i < 3; ++i) { if (zero[i]) ip[i] = 0; sum += ip[i]; }
    int last = -1;
    for (int i = 0; i < 3; ++i) if (!zero[i]) { ip[i] /= sum; last = i; }
    // Make the sum exactly 1: the last kept coordinate takes what the others leave.
    ip[last] = 1;
    for (int i = 0; i < 3; ++i) if (i != last) ip[last] -= ip[i];
    return zero[0] || zero[1] || zero[2];
  }
  
  
  // Given a segment find the maximum distance from it to the original surface. 
  // It is used to evaluate the Haustdorff distance of a Segment from the mesh.
  ScalarType MaxSegDist(VertexType *v0, VertexType *v1, CoordType &farthestPointOnSurf, CoordType &farthestN, Distribution<ScalarType> *distanceDistribution=0)
  {
    ScalarType maxSurfDist = 0;
    const ScalarType sampleNum = 10;
    for(ScalarType k = 1;k<sampleNum;++k)
    {
      ScalarType surfDist;
      CoordType closestPSurf;
      CoordType samplePnt = (v0->P()*k +v1->P()*(sampleNum-k))/sampleNum;          
      FaceType *f = Closest(samplePnt, closestPSurf, surfDist);
      if(distanceDistribution)
        distanceDistribution->Add(surfDist);
      if(surfDist > maxSurfDist)
      {
        maxSurfDist = surfDist;
        farthestPointOnSurf = closestPSurf;
        farthestN = f->N();
      }
    }
    return maxSurfDist;
  }
  
  
  /**
   * @brief RefineCurve
   * @param poly the curve to be refined
   * @param uniformFlag
   * 
   * Make one pass of refinement for all the edges of the curve that are distant from the basemesh
   * uses two parameters:
   * - par.minRefEdgeLen 
   * - par.surfDistThr
   */
    
  void RefineCurveByDistance(MeshType &poly)
  {
    tri::Allocator<MeshType>::CompactEveryVector(poly);    
    int startEdgeSize = poly.en;
    for(int i =0; i<startEdgeSize;++i)
    {
      EdgeType &ei = poly.edge[i];
      if(edge::Length(ei)>par.minRefEdgeLen)  
      {      
        CoordType farthestP, farthestN;
        ScalarType maxDist = MaxSegDist(ei.V(0),ei.V(1),farthestP, farthestN);
        if(maxDist > par.surfDistThr)  
        {
          edge::VEEdgeSplit(poly, &ei, farthestP, farthestN); 
        }
      }
    }
//    tri::Allocator<MeshType>::CompactEveryVector(poly);
//    printf("Refine %i -> %i\n",startEdgeSize,poly.en);fflush(stdout);
  }
  
  /**
   * \brief Make every polyline segment lie inside one face or along one edge of the base
   *        mesh, by inserting the exact points where it crosses the mesh edges.
   * \param poly The polyline; its own vertices are snapped first (SnapPolyline), then kept
   *
   * Each segment is traced across the mesh: in the current face it heads for its far end,
   * projected onto the face plane; the edge where it leaves is the first barycentric
   * coordinate to reach zero, and the crossing point is made with that coordinate exactly
   * zero, so it lies exactly on the edge. The trace then continues in the faces on the
   * other side, through a vertex if it leaves through one, until it reaches a face that
   * also holds the far end. Two strands traced through the same faces therefore keep their
   * order along every edge they cross: close parallel curves stay parallel, which
   * approximate crossing points (the bisection used here before) did not guarantee.
   *
   * \throws vcg::MissingPreconditionException if a segment cannot be traced, e.g. one that
   *         leaves the surface across a border
   */
  void RefineCurveByBaseMesh(MeshType &poly)
  {
    SnapPolyline(poly);  // compacts poly
    const int startEn = poly.en;
    MeshType out;
    tri::Allocator<MeshType>::AddVertices(out, poly.vert.size());
    for (size_t i = 0; i < poly.vert.size(); ++i) out.vert[i].ImportData(poly.vert[i]);
    for (size_t ei = 0; ei < poly.edge.size(); ++ei)
    {
      Progress(int(100 * ei / poly.edge.size()), "RefineCurveByBaseMesh: tracing segment %lu of %lu", (unsigned long)ei + 1, (unsigned long)poly.edge.size());
      const size_t i0 = tri::Index(poly, poly.edge[ei].cV(0)), i1 = tri::Index(poly, poly.edge[ei].cV(1));
      std::vector<CoordType> crossings;
      TraceSegment(poly.vert[i0].cP(), poly.vert[i1].cP(), crossings);
      size_t prev = i0;
      for (const CoordType &c : crossings)
      {
        const size_t vi = tri::Allocator<MeshType>::AddVertices(out, 1) - out.vert.begin();
        out.vert[vi].ImportData(poly.vert[i0]);
        out.vert[vi].P() = c;
        out.vert[vi].ClearS();  // a crossing point is not one of the locked vertices
        auto e = tri::Allocator<MeshType>::AddEdges(out, 1);
        e->V(0) = &out.vert[prev]; e->V(1) = &out.vert[vi];
        prev = vi;
      }
      auto e = tri::Allocator<MeshType>::AddEdges(out, 1);
      e->V(0) = &out.vert[prev]; e->V(1) = &out.vert[i1];
    }
    poly.Clear();
    tri::Append<MeshType, MeshType>::MeshCopy(poly, out);
    // A geodesic is straight inside a face and along an edge: drop the free vertices that
    // only bend a strand there, so it runs straight between its control points and edge
    // crossings. A vertex is dropped when every face it belongs to also holds both its
    // neighbours (inside a face: the neighbours are in it; on an edge: on the same edge),
    // so the shortcut stays where the strand was; strands meeting the same edges in the
    // same order cannot cross once straight.
    tri::UpdateTopology<MeshType>::VertexEdge(poly);
    int droppedCnt = 0;
    auto faces = [&](const CoordType &p) {
      std::vector<FaceType *> fs;
      for (const auto &l : Locate(p)) fs.push_back(l.first);
      std::sort(fs.begin(), fs.end());
      return fs;
    };
    for (VertexType &v : poly.vert)
    {
      if (v.IsD() || v.IsS() || edge::VEDegree<EdgeType>(&v) != 2) continue;
      std::vector<VertexPointer> nb;
      edge::VVStarVE(&v, nb);
      const std::vector<FaceType *> fv = faces(v.cP()), fa = faces(nb[0]->cP()), fb = faces(nb[1]->cP());
      if (!std::includes(fa.begin(), fa.end(), fv.begin(), fv.end()) ||
          !std::includes(fb.begin(), fb.end(), fv.begin(), fv.end())) continue;
      edge::VEEdgeCollapse(poly, &v);
      ++droppedCnt;
    }
    tri::Allocator<MeshType>::CompactEveryVector(poly);
    SimplifyNullEdges(poly);
    Log("RefineCurveByBaseMesh %i en -> %i en, %i vertices inside faces dropped", startEn, poly.en, droppedCnt);
  }

  /// Snap at rounding level only: a point this close to a vertex or an edge is on it.
  /// Shared by the trace and by CoMEmbed, so both agree on where a curve point is.
  bool RoundingSnap(CoordType &ip, const FaceType &f) const
  {
    return BarycentricSnap(ip, f, ScalarType(1e-4), base.bbox.Diag() * ScalarType(1e-6));
  }

  /// Where a point is on the mesh: every face it belongs to, with its barycentric
  /// coordinates there. One face for a point inside a face, the two faces of an edge, or
  /// the faces around a vertex.
  typedef std::vector<std::pair<FaceType *, CoordType>> Location;

  Location Locate(const CoordType &p)
  {
    CoordType ip;
    FaceType *f = GetClosestFaceIP(p, ip);
    RoundingSnap(ip, *f);
    return LocationFrom(f, ip);
  }

  Location LocationFrom(FaceType *f, const CoordType &ip) const
  {
    Location loc{{f, ip}};
    int zeros = 0, one = -1, zero = -1;
    for (int i = 0; i < 3; ++i) { if (ip[i] == 0) { ++zeros; zero = i; } else one = i; }
    if (zeros == 1)  // on the edge opposite to corner `zero`: add the face across it
    {
      const int e = (zero + 1) % 3;
      FaceType *g = f->FFp(e);
      if (g != f) {
        CoordType ig(0, 0, 0);
        for (int j = 0; j < 3; ++j) {
          if (g->V(j) == f->V(e)) ig[j] = ip[e];
          if (g->V(j) == f->V((e + 1) % 3)) ig[j] = ip[(e + 1) % 3];
        }
        loc.push_back({g, ig});
      }
    }
    else if (zeros == 2)  // on vertex V(one): every face around it, across edges incident to it
    {
      VertexPointer v = f->V(one);
      for (size_t k = 0; k < loc.size(); ++k)
        for (int i = 0; i < 3; ++i) {
          FaceType *h = loc[k].first;
          if (h->V(i) != v && h->V((i + 1) % 3) != v) continue;
          FaceType *g = h->FFp(i);
          bool seen = false;
          for (const auto &l : loc) seen |= l.first == g;
          if (seen) continue;
          CoordType ig(0, 0, 0);
          for (int j = 0; j < 3; ++j) if (g->V(j) == v) ig[j] = 1;
          loc.push_back({g, ig});
        }
    }
    return loc;
  }

  /// The points where the segment from \a p to \a q crosses mesh edges, in order.
  void TraceSegment(const CoordType &p, const CoordType &q, std::vector<CoordType> &crossings)
  {
    const Location target = Locate(q);
    Location cur = Locate(p);
    for (int step = 0; ; ++step)
    {
      for (const auto &c : cur) for (const auto &t : target)
        if (c.first == t.first) return;  // both ends in one face: the rest lies in it
      if (step > 4 * int(base.fn) + 16)
        throw vcg::MissingPreconditionException("CoM: could not trace a curve segment across the mesh.");
      // Among the faces of the current point, the one the segment heads into.
      FaceType *g = nullptr;
      CoordType bx, d;
      ScalarType bestInward = -1;
      for (const auto &c : cur)
      {
        CoordType bq;
        InterpolationParameters(*c.first, c.first->N(), q, bq);
        const CoordType dir = bq - c.second;
        ScalarType inward = std::numeric_limits<ScalarType>::max();
        for (int k = 0; k < 3; ++k) if (c.second[k] == 0) inward = std::min(inward, dir[k]);
        if (inward > 0 && inward > bestInward) { bestInward = inward; g = c.first; bx = c.second; d = dir; }
      }
      if (g == nullptr)
        throw vcg::MissingPreconditionException("CoM: a curve segment leaves the surface across a border.");
      // Leave g where the first coordinate reaches zero.
      int exitK = -1;
      ScalarType sExit = std::numeric_limits<ScalarType>::max();
      for (int k = 0; k < 3; ++k)
        if (d[k] < 0 && bx[k] > 0 && bx[k] / -d[k] < sExit) { sExit = bx[k] / -d[k]; exitK = k; }
      if (exitK < 0)
        throw vcg::MissingPreconditionException("CoM: could not trace a curve segment across the mesh.");
      CoordType b = bx + d * sExit;
      b[exitK] = 0;
      const ScalarType rest = b[(exitK + 1) % 3] + b[(exitK + 2) % 3];
      b[(exitK + 1) % 3] /= rest;
      b[(exitK + 2) % 3] = 1 - b[(exitK + 1) % 3];
      RoundingSnap(b, *g);  // through a vertex, if it passes that close to one
      crossings.push_back(g->cP(0) * b[0] + g->cP(1) * b[1] + g->cP(2) * b[2]);
      cur = LocationFrom(g, b);
    }
  }
  
  
  /**
   * @brief LaplacianFunctor basic Laplacian smoothing functor
   *
   * It computes the desired position for each vertex as the average of its
   * current position and the positions of its 1-ring neighbors. Used as the
   * position functor in the SmoothProject function.
   */
  struct LaplacianFunctor {
    std::vector<CoordType> operator()(const MeshType &poly) const {
      std::vector<CoordType> posVec(poly.vn, CoordType(0,0,0));
      std::vector<int>       cntVec(poly.vn, 0);
      for(int i=0; i<poly.en; ++i)
        for(int j=0; j<2; ++j) {
          int vi = tri::Index(poly, poly.edge[i].V0(j));
          posVec[vi] += poly.edge[i].V1(j)->P();
          cntVec[vi] += 1;
        }
      for(int i=0; i<poly.vn; ++i)
        posVec[i] = (poly.vert[i].P() + posVec[i]) / ScalarType(cntVec[i]+1);
      return posVec;
    }
  };

  /**
   * @brief QualityDistanceFieldFunctor quality field based smoothing functor
   *
   * @param poly the input curve mesh
   * @param com the CurveOnManifold class itself to quick access closest face 
   * @param scale a scaling factor to control the step size of the movement
   * along the quality, if 0 it will be automatically set to 1/2 of the CoM parameter `par.maxMoveDelta` (default is 0)
   * @param smoothBlend a blending factor to control the influence of the
   * smoothing (default is 0.5)
   *
   * It compute the new position using the quality field of the mesh assuming
   * that it is a distance field sampled per vertices and that we would like to move toward the zero of the distance field. 
   * We use gradient of the quality field for the direction and the sign of the quality for the versus of the direction. 
   * We move of a quantity proportional to the quality value at the vertex.
   *
   */

  struct QualityDistanceFieldFunctor
  {

    CoM<MeshType> &com;
    ScalarType scale;
    ScalarType smoothBlend = 0.5;
    QualityDistanceFieldFunctor(CoM<MeshType> &_com, ScalarType _scale=0, ScalarType _smoothBlend = 0.5) : com(_com), scale(_scale), smoothBlend(_smoothBlend) {};

    std::vector<CoordType> operator()(const MeshType &poly) const
    {
      // Step 1: Compute smoothed position using Laplacian smoothing
      std::vector<CoordType> smoothPosVec(poly.vn, CoordType(0, 0, 0));
      std::vector<int> cntVec(poly.vn, 0);
      for (int i = 0; i < poly.en; ++i)
        for (int j = 0; j < 2; ++j)
        {
          int vi = tri::Index(poly, poly.edge[i].V0(j));
          smoothPosVec[vi] += poly.edge[i].V1(j)->P();
          cntVec[vi] += 1;
        }
      for (int i = 0; i < poly.vn; ++i)
        smoothPosVec[i] = (poly.vert[i].P() + smoothPosVec[i]) / ScalarType(cntVec[i] + 1);

      // Step 2: Compute field-based position using the gradient of the quality field and moving toward zero
      std::vector<CoordType> fieldPosVec(poly.vn, CoordType(0, 0, 0));
      for (const VertexType &v : poly.vert)
      {
        CoordType ip; 
        FacePointer f = com.GetClosestFaceIP(v.P(),ip); // throws rather than return null
        ScalarType q = f->V(0)->Q() * ip[0] + f->V(1)->Q() * ip[1] + f->V(2)->Q() * ip[2];        
        CoordType fieldDir = GradientScalarField(*f, f->V(0)->Q(),f->V(1)->Q(),f->V(2)->Q());
        
        fieldPosVec[tri::Index(poly, v)] = v.P() + fieldDir * scale * q; // Move towards higher quality (lower distance)
      }
      std::vector<CoordType> PosVec(poly.vn, CoordType(0, 0, 0));
      for (int i = 0; i < poly.vn; ++i)
        PosVec[i] = (smoothPosVec[i] * smoothBlend + fieldPosVec[i] * (1.0 - smoothBlend));

      return PosVec;
    }
  };

  /**
   * @brief MoveAndProject
   * @param poly
   * @param iterNum
   * @param moveWeight    [0..1] blend toward the desired position returned by the functor
   * @param projectWeight [0..1] blend toward the closest point on the surface
   * @param desiredPos    functor: std::vector<CoordType>(const MeshType&)
   *
   * Generic version of SmoothProject: the per-vertex desired position is
   * supplied by a functor instead of being hard-coded as a Laplacian average.
   */
  template<typename PositionFunctor>
  void MoveAndProject(MeshType &poly, int iterNum, ScalarType moveWeight, ScalarType projectWeight,
                      PositionFunctor &desiredPos)
  {
    tri::RequireCompactness(poly);
    tri::UpdateTopology<MeshType>::VertexEdge(poly);
    if (poly.en == 0 || base.fn == 0)
      throw vcg::MissingPreconditionException("CoM: the curve and the base mesh must both be non-empty.");
    for(int k=0;k<iterNum;++k)
    {
      Progress(100*k/iterNum, "MoveAndProject iteration %i of %i, %i vertices", k+1, iterNum, poly.vn);
      if(k==iterNum-1) projectWeight=1;
      std::vector<CoordType> desired = desiredPos(poly);

      for(int i=0; i<poly.vn; ++i)
        if(!poly.vert[i].IsS())
        {
          // Clamp the movement towards the desired position by maxMoveDist
          CoordType delta = desired[i] - poly.vert[i].P();
          ScalarType deltaLen = delta.Norm();
          if(deltaLen > par.maxMoveDelta) {
            delta *= (par.maxMoveDelta / deltaLen);
            desired[i] = poly.vert[i].P() + delta;
          } 

          CoordType newP = poly.vert[i].P()*(1.0-moveWeight) + desired[i]*moveWeight;
          
          CoordType closestP;
          FaceType *f = GetClosestFacePoint(newP, closestP);
          poly.vert[i].P() = newP*(1.0-projectWeight) +closestP*projectWeight;
          poly.vert[i].N() = f->N();
        }
        else  // a control point is only projected
        {
          CoordType closestP;
          FaceType *f = GetClosestFacePoint(poly.vert[i].P(), closestP);
          poly.vert[i].P() = closestP;
          poly.vert[i].N() = f->N();
        }
      
      tri::UpdateTopology<MeshType>::TestVertexEdge(poly);
      RefineCurveByDistance(poly);      
      tri::UpdateTopology<MeshType>::TestVertexEdge(poly);
      Simplify(poly);
      tri::UpdateTopology<MeshType>::TestVertexEdge(poly);
      CollapseNullEdges(poly);
    }
  }

  void SmoothProject(MeshType &poly, int iterNum, ScalarType smoothWeight, ScalarType projectWeight)
  {
    LaplacianFunctor lapFunct;
    MoveAndProject(poly, iterNum, smoothWeight, projectWeight, lapFunct);
  }


};

/** \ingroup trimesh
 * \brief Embed a curve in the surface it lies on: the operations on a CoM that change the surface.
 *
 * CoM is a query structure over a fixed surface; what only reads the surface, or changes
 * the curve, lives there. These two change the base mesh itself, its connectivity and its
 * face-edge selection, so they are kept apart. Both take the CoM already built over the
 * surface, typically the one just used to smooth and refine the curve, so its spatial grid
 * is not built twice. SplitMeshWithPolyline() re-initializes it on the split surface, so it
 * stays valid afterwards.
 *
 * \code
 * tri::CoM<MyMesh> com(base);
 * com.Init();
 * com.SmoothProject(poly, 10, 0.5, 0.5);
 * com.RefineCurveByBaseMesh(poly);
 * tri::CoMEmbed<MyMesh>::SplitMeshWithPolyline(com, poly);
 * tri::CoMEmbed<MyMesh>::TagFaceEdgeSelWithPolyLine(com, poly);
 * tri::CutMeshAlongSelectedFaceEdges(base);
 * \endcode
 * See the sample trimesh_topological_cut.cpp.
 */
template <class MeshType>
class CoMEmbed
{
public:
  typedef CoM<MeshType>                     CoMType;
  typedef typename MeshType::ScalarType     ScalarType;
  typedef typename MeshType::CoordType      CoordType;
  typedef typename MeshType::VertexType     VertexType;
  typedef typename MeshType::VertexPointer  VertexPointer;
  typedef typename MeshType::VertexIterator VertexIterator;
  typedef typename MeshType::EdgeIterator   EdgeIterator;
  typedef typename MeshType::FaceType       FaceType;
  typedef typename MeshType::FacePointer    FacePointer;
  typedef typename MeshType::FaceIterator   FaceIterator;

  
  /**
   * \brief Tag face edges of the base mesh that coincide with polyline edges
   * \param com The CoM built over the base mesh
   * \param poly The polyline (as edge mesh) to use for tagging
   * \param markFlag If true, clears all FaceEdgeS flags before tagging (default: true)
   * \return true if ALL edges of the polyline are fully snapped onto mesh edges
   * 
   * This function marks edges in the base mesh (using FaceEdgeS flag) where they 
   * coincide with edges from the polyline. The polyline edges must be snapped to 
   * mesh vertices at both endpoints for this to work.
   * This function requires VertexFace and FaceFace adjacency.
   * 
   * \note This is typically used as a preparation step before cutting the mesh
   *       along the polyline using CutMeshAlongCrease or similar functions.
   *      
   * 
   * \warning Returns false if any polyline edge is not properly snapped, or if
   *          the snapped vertices don't form an edge in the base mesh.
   * 
   * \sa SplitMeshWithPolyline, CutMeshAlongCrease
   */
    
static bool TagFaceEdgeSelWithPolyLine(CoMType &com, MeshType &poly,bool markFlag=true)
{
	if (markFlag)
		tri::UpdateFlags<MeshType>::FaceClearFaceEdgeS(com.base);

	tri::UpdateTopology<MeshType>::VertexFace(com.base);
	tri::UpdateTopology<MeshType>::FaceFace(com.base);

	for(EdgeIterator ei=poly.edge.begin(); ei!=poly.edge.end();++ei)
	{
		CoordType ip0,ip1;
		FaceType *f0 = com.GetClosestFaceIP(ei->cP(0),ip0);
		FaceType *f1 = com.GetClosestFaceIP(ei->cP(1),ip1);

		if(com.BarycentricSnap(ip0, *f0) && com.BarycentricSnap(ip1, *f1))
		{
			VertexPointer v0 = com.FindVertexSnap(f0,ip0);
			VertexPointer v1 = com.FindVertexSnap(f1,ip1);

			if(v0==0 || v1==0)
				return false;
			if(v0==v1)
				return false;

			FacePointer ff0,ff1;
			int e0,e1;
			bool ret=face::FindSharedFaces<FaceType>(v0,v1,ff0,ff1,e0,e1);
			if(ret)
			{
				assert(ret);
				assert(ff0->V(e0)==v0 || ff0->V(e0)==v1);
				ff0->SetFaceEdgeS(e0);
				ff1->SetFaceEdgeS(e1);
			} else {
				return false;
			}
		}
		else {
			return false;
		}
	}
	return true;
}

  
  /**
   * \brief Make the base mesh conform to the polyline: afterwards every polyline edge is an
   *        edge of the base mesh, with its face-edge selection bit set on both sides.
   * \param com The CoM built over the base mesh, typically the one just used to refine the
   *            polyline; it is re-initialized on the split mesh, so it stays valid
   * \param poly The polyline; its vertices are moved onto the vertices created for them
   *
   * Precondition: RefineCurveByBaseMesh() has been called, so every polyline segment lies
   * inside one face or along one edge of the base mesh.
   *
   * Curves may cross and touch. The mesh realizes the union of the curves: it gets a vertex
   * at every crossing, and curve points or segments that coincide within rounding distance
   * share mesh vertices and edges. The polyline keeps its own topology: a segment cut at a
   * crossing or a touching point gets a vertex of its own there, so two strands stay two
   * strands, each passing through the shared mesh vertex. CoM::ConnectCrossings() turns
   * those meeting points into junctions, when the caller wants them.
   *
   * Any number of segments may cross the same triangle, from different curves or from the
   * same one passing several times, as a spiral on a cylinder does. Each crossed triangle is
   * rebuilt on its own, in its barycentric frame: the polyline points on its edges are
   * inserted by splitting the sub-triangle that owns that piece of the edge, points inside
   * it by splitting the sub-triangle that contains them, and every segment is then recovered
   * as an edge by flipping the edges that cross it (S. W. Sloan, "A fast algorithm for
   * generating constrained Delaunay triangulations", Computers & Structures 47(3), 1993; the
   * Delaunay part is not needed here). Sub-triangles keep the attributes of their triangle,
   * and their wedge texture coordinates are the triangle's at each vertex; new vertices
   * interpolate the attributes of the surface where they are (VertexInterpolator).
   *
   * \throws vcg::MissingPreconditionException if a segment does not lie in one face or
   *         along one edge, or if a segment cannot be recovered
   * \sa CoM::RefineCurveByBaseMesh, TagFaceEdgeSelWithPolyLine
   */
  static void SplitMeshWithPolyline(CoMType &com, MeshType &poly)
  {
    MeshType &m = com.base;
    tri::Allocator<MeshType>::CompactEveryVector(poly);
    if (poly.vn == 0) return;
    typedef std::pair<size_t, size_t> EdgeKey;  // mesh vertex indices, smaller first
    auto key = [](size_t a, size_t b) { return a < b ? EdgeKey(a, b) : EdgeKey(b, a); };

    // 1. Where each polyline vertex is on the mesh: on a vertex, on an edge, or in a face.
    enum { OnVertex, OnEdge, InFace };
    struct Loc { int kind; size_t vert; EdgeKey edge; ScalarType t; size_t face[2]; int faceNum; CoordType ip; };
    std::vector<Loc> loc(poly.vert.size());
    for (size_t ci = 0; ci < poly.vert.size(); ++ci)
    {
      Loc &l = loc[ci];
      FaceType *f = com.GetClosestFaceIP(poly.vert[ci].cP(), l.ip);
      // Only rounding-level snapping here: SnapPolyline already put the vertices that may
      // be snapped exactly on their vertex or edge, and kept the others off on purpose.
      com.RoundingSnap(l.ip, *f);
      int z = -1, one = -1;
      for (int i = 0; i < 3; ++i) { if (l.ip[i] == 0) z = i; if (l.ip[i] == 1) one = i; }
      l.face[0] = tri::Index(m, f); l.faceNum = 1;
      if (one >= 0) { l.kind = OnVertex; l.vert = tri::Index(m, f->V(one)); }
      else if (z >= 0) {
        // On the edge opposite to corner z: from V(z+1), at parameter ip[z+2].
        l.kind = OnEdge;
        const size_t a = tri::Index(m, f->V((z + 1) % 3)), b = tri::Index(m, f->V((z + 2) % 3));
        l.edge = key(a, b);
        l.t = (a < b) ? l.ip[(z + 2) % 3] : l.ip[(z + 1) % 3];
        const int e = (z + 1) % 3;  // index of that edge in f
        if (f->FFp(e) != f) { l.face[1] = tri::Index(m, f->FFp(e)); l.faceNum = 2; }
      }
      else l.kind = InFace;
    }

    // 2. A mesh vertex for each of them: the existing one, or a new one shared by every
    //    polyline vertex within rounding distance (the mesh realizes the union of the
    //    curves: two strands meeting at a point meet at one mesh vertex).
    const ScalarType tol = m.bbox.Diag() * ScalarType(1e-6);
    std::vector<size_t> meshV(poly.vert.size());
    std::map<EdgeKey, std::vector<std::pair<ScalarType, size_t>>> edgePoints;  // (t from first, vertex)
    std::map<size_t, std::vector<std::pair<size_t, CoordType>>> facePoints;   // (vertex, barycentric)
    size_t newVertNum = 0;
    const size_t firstNewVert = m.vert.size();
    for (size_t ci = 0; ci < poly.vert.size(); ++ci)
    {
      const Loc &l = loc[ci];
      if (l.kind == OnVertex) { meshV[ci] = l.vert; continue; }
      if (l.kind == OnEdge) {
        auto &pts = edgePoints[l.edge];
        const ScalarType len = Distance(m.vert[l.edge.first].cP(), m.vert[l.edge.second].cP());
        auto it = std::find_if(pts.begin(), pts.end(), [&](const std::pair<ScalarType, size_t> &p) { return std::abs(p.first - l.t) * len <= tol; });
        if (it != pts.end()) { meshV[ci] = it->second; continue; }
        meshV[ci] = firstNewVert + newVertNum++;
        pts.push_back({l.t, meshV[ci]});
      } else {
        auto &pts = facePoints[l.face[0]];
        const FaceType &f = m.face[l.face[0]];
        const CoordType p = f.cP(0) * l.ip[0] + f.cP(1) * l.ip[1] + f.cP(2) * l.ip[2];
        auto it = std::find_if(pts.begin(), pts.end(), [&](const std::pair<size_t, CoordType> &q) {
          return Distance(p, f.cP(0) * q.second[0] + f.cP(1) * q.second[1] + f.cP(2) * q.second[2]) <= tol; });
        if (it != pts.end()) { meshV[ci] = it->first; continue; }
        meshV[ci] = firstNewVert + newVertNum++;
        pts.push_back({meshV[ci], l.ip});
      }
    }
    tri::Allocator<MeshType>::AddVertices(m, newVertNum);
    for (size_t ci = 0; ci < poly.vert.size(); ++ci)
    {
      const Loc &l = loc[ci];
      VertexType &nv = m.vert[meshV[ci]];
      if (l.kind == OnEdge && meshV[ci] >= firstNewVert) {
        const VertexType &va = m.vert[l.edge.first], &vb = m.vert[l.edge.second];
        nv.P() = va.cP() * (1 - l.t) + vb.cP() * l.t;
        VertexInterpolator<MeshType>::Lerp(m, nv, va, vb, l.t);
      } else if (l.kind == InFace && meshV[ci] >= firstNewVert) {
        const FaceType &f = m.face[l.face[0]];
        nv.P() = f.cP(0) * l.ip[0] + f.cP(1) * l.ip[1] + f.cP(2) * l.ip[2];
        VertexInterpolator<MeshType>::Barycentric(m, nv, *f.cV(0), *f.cV(1), *f.cV(2), l.ip);
      }
    }
    for (size_t ci = 0; ci < poly.vert.size(); ++ci) poly.vert[ci].P() = m.vert[meshV[ci]].cP();

    // 3. The segments: across one face, or along one edge.
    struct Segment { size_t polyEdge, a, b; };  // mesh vertices, in the polyline edge's direction
    std::map<size_t, std::vector<Segment>> faceSegments;
    std::vector<std::pair<Segment, EdgeKey>> edgeSegments;
    std::map<size_t, std::vector<size_t>> vertFaces;  // faces around the polyline vertices on mesh vertices
    for (size_t ci = 0; ci < poly.vert.size(); ++ci)
      if (loc[ci].kind == OnVertex) vertFaces[loc[ci].vert];
    if (!vertFaces.empty())
      for (size_t fi = 0; fi < m.face.size(); ++fi) if (!m.face[fi].IsD())
        for (int i = 0; i < 3; ++i) {
          auto it = vertFaces.find(tri::Index(m, m.face[fi].V(i)));
          if (it != vertFaces.end()) it->second.push_back(fi);
        }
    auto candidateFaces = [&](size_t ci) {
      if (loc[ci].kind == OnVertex) return vertFaces[loc[ci].vert];
      return std::vector<size_t>(loc[ci].face, loc[ci].face + loc[ci].faceNum);
    };
    for (size_t ei = 0; ei < poly.edge.size(); ++ei)
    {
      const auto &e = poly.edge[ei];
      const size_t c0 = tri::Index(poly, e.cV(0)), c1 = tri::Index(poly, e.cV(1));
      if (meshV[c0] == meshV[c1]) continue;  // a zero-length segment
      // The segment lies in the faces both its ends belong to: one face for a segment
      // across it, the two faces of an edge for a segment along it.
      std::vector<size_t> f0 = candidateFaces(c0), f1 = candidateFaces(c1), common;
      std::sort(f0.begin(), f0.end()); std::sort(f1.begin(), f1.end());
      std::set_intersection(f0.begin(), f0.end(), f1.begin(), f1.end(), std::back_inserter(common));
      if (common.empty())
        throw vcg::MissingPreconditionException("CoMEmbed: a curve segment does not lie in one face; call RefineCurveByBaseMesh first.");
      const Segment s{ei, meshV[c0], meshV[c1]};
      if (common.size() == 1) faceSegments[common[0]].push_back(s);
      else {
        const FaceType &f = m.face[common[0]];  // the edge both faces share
        EdgeKey ek(0, 0);
        for (int i = 0; i < 3; ++i)
          if (f.cFFp(i) == &m.face[common[1]]) ek = key(tri::Index(m, f.cV(i)), tri::Index(m, f.cV1(i)));
        edgeSegments.push_back({s, ek});
      }
    }

    std::set<size_t> touched;
    for (const auto &fp : facePoints) touched.insert(fp.first);
    for (const auto &fs : faceSegments) touched.insert(fs.first);
    for (size_t fi = 0; fi < m.face.size(); ++fi) if (!m.face[fi].IsD())
      for (int i = 0; i < 3; ++i)
        if (edgePoints.count(key(tri::Index(m, m.face[fi].V(i)), tri::Index(m, m.face[fi].V1(i))))) touched.insert(fi);

    // 4. Rebuild every touched face on its own. Where segments cross, the crossing point
    //    becomes a mesh vertex; every segment is cut at the curve points lying on it
    //    (crossings, other curves touching or overlapping it), and the pieces recovered.
    std::map<size_t, std::vector<size_t>> chains;  // polyline edge -> mesh vertices inside it, in order
    std::vector<std::pair<size_t, CoordType>> crossingVerts;  // (face, barycentric) of the new vertices
    const size_t firstCrossingVert = m.vert.size();
    std::vector<std::pair<size_t, LocalTriangulation>> rebuilt;
    size_t newFaceNum = 0, done = 0;
    for (size_t fi : touched)
    {
      com.Progress(int(100 * done / touched.size()), "SplitMeshWithPolyline: rebuilding face %lu of %lu", (unsigned long)done + 1, (unsigned long)touched.size());
      ++done;
      LocalTriangulation lt;
      const FaceType &f = m.face[fi];
      size_t corner[3];
      for (int k = 0; k < 3; ++k) { corner[k] = tri::Index(m, f.cV(k)); lt.AddPoint(corner[k], k == 1, k == 2, (1 << k) | (1 << ((k + 2) % 3))); }
      lt.tris.push_back({{0, 1, 2}});
      for (int k = 0; k < 3; ++k)
      {
        auto it = edgePoints.find(key(corner[k], corner[(k + 1) % 3]));
        if (it == edgePoints.end()) continue;
        std::vector<std::pair<ScalarType, size_t>> pts;  // parameter from corner k
        for (const auto &p : it->second) pts.push_back({corner[k] < corner[(k + 1) % 3] ? p.first : 1 - p.first, p.second});
        std::sort(pts.begin(), pts.end());
        int prev = k;
        for (const auto &p : pts)
        {
          const ScalarType s = p.first;
          const ScalarType u = lt.pts[k].u * (1 - s) + lt.pts[(k + 1) % 3].u * s;
          const ScalarType v = lt.pts[k].v * (1 - s) + lt.pts[(k + 1) % 3].v * s;
          const int pi = lt.AddPoint(p.second, u, v, 1 << k);
          lt.SplitEdge(prev, (k + 1) % 3, pi);
          prev = pi;
        }
      }
      for (const auto &p : facePoints[fi]) lt.InsertInterior(lt.AddPoint(p.first, p.second[1], p.second[2], 0));

      const std::vector<Segment> &segs = faceSegments[fi];
      auto pos3 = [&](int i) { return f.cP(0) * (1 - lt.pts[i].u - lt.pts[i].v) + f.cP(1) * lt.pts[i].u + f.cP(2) * lt.pts[i].v; };
      auto nearPoint = [&](const CoordType &p) {
        for (size_t i = 0; i < lt.pts.size(); ++i) if (Distance(pos3(int(i)), p) <= tol) return int(i);
        return -1;
      };
      // Crossings: a new mesh vertex where two segments cross, unless a point is already there.
      for (size_t i = 0; i < segs.size(); ++i)
        for (size_t j = i + 1; j < segs.size(); ++j)
        {
          const int p = lt.Local(segs[i].a), q = lt.Local(segs[i].b), a = lt.Local(segs[j].a), b = lt.Local(segs[j].b);
          if (!lt.Crosses(p, q, a, b)) continue;
          const long double op = lt.Orient(a, b, p), oq = lt.Orient(a, b, q);
          const ScalarType s = ScalarType(op / (op - oq));
          const ScalarType u = lt.pts[p].u + (lt.pts[q].u - lt.pts[p].u) * s, v = lt.pts[p].v + (lt.pts[q].v - lt.pts[p].v) * s;
          const CoordType x = f.cP(0) * (1 - u - v) + f.cP(1) * u + f.cP(2) * v;
          if (nearPoint(x) >= 0) continue;
          crossingVerts.push_back({fi, CoordType(1 - u - v, u, v)});
          lt.InsertInterior(lt.AddPoint(firstCrossingVert + crossingVerts.size() - 1, u, v, 0));
        }
      // Each segment cut at the curve points on it, then its pieces made edges.
      for (const Segment &sg : segs)
      {
        const int p = lt.Local(sg.a), q = lt.Local(sg.b);
        const Segment3<ScalarType> seg(pos3(p), pos3(q));
        std::vector<std::pair<ScalarType, int>> on;  // (parameter along the segment, point)
        for (size_t r = 0; r < lt.pts.size(); ++r)
        {
          if (int(r) == p || int(r) == q) continue;
          CoordType closest; ScalarType dist;
          SegmentPointDistance(seg, pos3(int(r)), closest, dist);
          const ScalarType t = (pos3(int(r)) - seg.P0()) * (seg.P1() - seg.P0()) / seg.SquaredLength();
          if (dist <= tol && t > 0 && t < 1 && Distance(pos3(int(r)), seg.P0()) > tol && Distance(pos3(int(r)), seg.P1()) > tol)
            on.push_back({t, int(r)});
        }
        std::sort(on.begin(), on.end());
        int prev = p;
        for (const auto &o : on) { lt.Recover(prev, o.second); prev = o.second; chains[sg.polyEdge].push_back(lt.pts[o.second].vert); }
        lt.Recover(prev, q);
      }
      newFaceNum += lt.tris.size() - 1;
      rebuilt.push_back({fi, std::move(lt)});
    }
    // A segment along an edge is cut at the curve points on that edge between its ends.
    for (const auto &es : edgeSegments)
    {
      const Segment &sg = es.first;
      auto tOf = [&](size_t v) -> ScalarType {
        if (v == es.second.first) return 0;
        if (v == es.second.second) return 1;
        for (const auto &p : edgePoints[es.second]) if (p.second == v) return p.first;
        return -1;
      };
      const ScalarType ta = tOf(sg.a), tb = tOf(sg.b);
      std::vector<std::pair<ScalarType, size_t>> between;
      for (const auto &p : edgePoints[es.second])
        if (p.second != sg.a && p.second != sg.b && (p.first - ta) * (p.first - tb) < 0)
          between.push_back({std::abs(p.first - ta), p.second});
      std::sort(between.begin(), between.end());
      for (const auto &b : between) chains[sg.polyEdge].push_back(b.second);
    }
    tri::Allocator<MeshType>::AddVertices(m, crossingVerts.size());
    for (size_t k = 0; k < crossingVerts.size(); ++k)
    {
      const FaceType &f = m.face[crossingVerts[k].first];
      const CoordType &b = crossingVerts[k].second;
      VertexType &nv = m.vert[firstCrossingVert + k];
      nv.P() = f.cP(0) * b[0] + f.cP(1) * b[1] + f.cP(2) * b[2];
      VertexInterpolator<MeshType>::Barycentric(m, nv, *f.cV(0), *f.cV(1), *f.cV(2), b);
    }

    // 5. Write the sub-triangles back: the first in place of the face, the others new.
    const bool wedgeTex = tri::HasPerWedgeTexCoord(m);
    size_t nextFace = m.face.size();
    tri::Allocator<MeshType>::AddFaces(m, newFaceNum);
    for (auto &r : rebuilt)
    {
      FaceType &f = m.face[r.first];
      typename FaceType::TexCoordType wt[3];
      bool faux[3], border[3], edgeSel[3];
      for (int k = 0; k < 3; ++k) {
        if (wedgeTex) wt[k] = f.WT(k);
        faux[k] = f.IsF(k); border[k] = f.IsB(k); edgeSel[k] = f.IsFaceEdgeS(k);
      }
      const LocalTriangulation &lt = r.second;
      for (size_t ti = 0; ti < lt.tris.size(); ++ti)
      {
        FaceType &g = (ti == 0) ? f : m.face[nextFace++];
        if (ti > 0) g.ImportData(f);
        for (int i = 0; i < 3; ++i)
        {
          const auto &a = lt.pts[lt.tris[ti][i]], &b = lt.pts[lt.tris[ti][(i + 1) % 3]];
          g.V(i) = &m.vert[a.vert];
          if (wedgeTex) {
            g.WT(i) = wt[0];
            g.WT(i).P() = wt[0].P() * (1 - a.u - a.v) + wt[1].P() * a.u + wt[2].P() * a.v;
          }
          // An edge along an edge of the original face keeps that edge's bits; the
          // others are inside it.
          const int shared = a.edges & b.edges;
          const int k = shared ? (shared & 1 ? 0 : (shared & 2 ? 1 : 2)) : -1;
          if (k >= 0 && faux[k]) g.SetF(i); else g.ClearF(i);
          if (k >= 0 && border[k]) g.SetB(i); else g.ClearB(i);
          if (k >= 0 && edgeSel[k]) g.SetFaceEdgeS(i); else g.ClearFaceEdgeS(i);
        }
      }
    }
    com.Init();

    // 6. The polyline follows: a segment cut at some points gets its own vertices there, so
    //    two strands crossing or touching stay two strands, each with a vertex at the
    //    shared mesh vertex (CoM::ConnectCrossings merges them, if junctions are wanted).
    //    Every polyline edge is now a mesh edge: select it on both sides.
    std::set<EdgeKey> curveEdges;
    const size_t oldEdgeNum = poly.edge.size();
    for (size_t ei = 0; ei < oldEdgeNum; ++ei)
    {
      if (poly.edge[ei].IsD()) continue;
      const size_t c0 = tri::Index(poly, poly.edge[ei].cV(0)), c1 = tri::Index(poly, poly.edge[ei].cV(1));
      size_t prev = c0, prevMV = meshV[c0];
      auto ch = chains.find(ei);
      if (ch != chains.end())
        for (size_t mv : ch->second)
        {
          const size_t vi = tri::Allocator<MeshType>::AddVertices(poly, 1) - poly.vert.begin();
          poly.vert[vi].ImportData(poly.vert[c0]);
          poly.vert[vi].ClearS();  // a cut point is not a control point
          poly.vert[vi].P() = m.vert[mv].cP();
          auto ne = tri::Allocator<MeshType>::AddEdges(poly, 1);
          ne->V(0) = &poly.vert[prev]; ne->V(1) = &poly.vert[vi];
          curveEdges.insert(key(prevMV, mv));
          prev = vi; prevMV = mv;
        }
      poly.edge[ei].V(0) = &poly.vert[prev];  // the original edge is the last piece
      if (prevMV != meshV[c1]) curveEdges.insert(key(prevMV, meshV[c1]));
    }
    for (FaceType &f : m.face) if (!f.IsD())
      for (int i = 0; i < 3; ++i)
        if (curveEdges.count(key(tri::Index(m, f.V0(i)), tri::Index(m, f.V1(i))))) f.SetFaceEdgeS(i);
  }

private:
  /// The triangulation of one face of the base mesh being rebuilt, in its barycentric frame:
  /// a point at (u, v) is V0 + u (V1 - V0) + v (V2 - V0). Triangles are counterclockwise in
  /// that frame, as the face is.
  struct LocalTriangulation
  {
    struct Point { size_t vert; ScalarType u, v; int edges; };  // edges: bit k if on face edge k
    std::vector<Point> pts;
    std::vector<std::array<int, 3>> tris;
    std::set<std::pair<int, int>> recovered;  // segments already made edges, smaller index first
    static constexpr long double eps = 1e-13;  // orientation tolerance, in the unit frame

    int AddPoint(size_t vert, ScalarType u, ScalarType v, int edges)
    {
      pts.push_back({vert, u, v, edges});
      return int(pts.size()) - 1;
    }
    int Local(size_t vert) const
    {
      for (size_t i = 0; i < pts.size(); ++i) if (pts[i].vert == vert) return int(i);
      throw vcg::MissingPreconditionException("CoMEmbed: a curve segment ends outside the face it was assigned to.");
    }
    long double Orient(int a, int b, int c) const
    {
      return planar_polygon_detail::Orient2D(Point2d(pts[a].u, pts[a].v), Point2d(pts[b].u, pts[b].v), Point2d(pts[c].u, pts[c].v));
    }
    /// The triangle with the directed edge (a,b), and the position of a in it.
    bool FindEdge(int a, int b, size_t &ti, int &i) const
    {
      for (ti = 0; ti < tris.size(); ++ti)
        for (i = 0; i < 3; ++i)
          if (tris[ti][i] == a && tris[ti][(i + 1) % 3] == b) return true;
      return false;
    }
    /// Split the edge (a,b) at p, which lies on it, in the triangles on both its sides.
    void SplitEdge(int a, int b, int p)
    {
      for (int side = 0; side < 2; ++side)
      {
        size_t ti; int i;
        if (!FindEdge(side ? b : a, side ? a : b, ti, i)) continue;
        const std::array<int, 3> t = tris[ti];
        tris[ti] = {{t[i], p, t[(i + 2) % 3]}};
        tris.push_back({{p, t[(i + 1) % 3], t[(i + 2) % 3]}});
      }
    }
    /// Insert a point strictly inside the face.
    void InsertInterior(int p)
    {
      for (size_t ti = 0; ti < tris.size(); ++ti)
      {
        const std::array<int, 3> t = tris[ti];
        long double o[3];
        for (int i = 0; i < 3; ++i) o[i] = Orient(t[i], t[(i + 1) % 3], p);
        if (o[0] < -eps || o[1] < -eps || o[2] < -eps) continue;
        int onEdge = -1, onCount = 0;
        for (int i = 0; i < 3; ++i) if (o[i] <= eps) { onEdge = i; ++onCount; }
        if (onCount > 1)
          throw vcg::MissingPreconditionException("CoMEmbed: two curve points coincide inside a face.");
        if (onEdge >= 0) { SplitEdge(t[onEdge], t[(onEdge + 1) % 3], p); return; }
        tris[ti] = {{t[0], t[1], p}};
        tris.push_back({{t[1], t[2], p}});
        tris.push_back({{t[2], t[0], p}});
        return;
      }
      throw vcg::MissingPreconditionException("CoMEmbed: a curve point is outside the face it was assigned to.");
    }
    static int Sign(long double o) { return o > eps ? 1 : (o < -eps ? -1 : 0); }
    /// True if the segments (p,q) and (a,b) cross at a point inside both.
    bool Crosses(int p, int q, int a, int b) const
    {
      return Sign(Orient(p, q, a)) * Sign(Orient(p, q, b)) < 0 && Sign(Orient(a, b, p)) * Sign(Orient(a, b, q)) < 0;
    }
    /// Make (p,q) an edge by Sloan's flips: an edge crossing it is flipped when its two
    /// triangles form a strictly convex quad, and requeued otherwise; a new edge that still
    /// crosses is queued again. Curves that do not cross never queue a recovered segment;
    /// one that does is reported rather than flipped away.
    void Recover(int p, int q)
    {
      std::deque<std::pair<int, int>> queue;
      for (const auto &t : tris)
        for (int i = 0; i < 3; ++i)
        {
          const int a = t[i], b = t[(i + 1) % 3];
          if (a < b && Crosses(p, q, a, b)) {
            if (recovered.count({a, b}))
              throw vcg::MissingPreconditionException("CoMEmbed: two curve segments cross inside a face.");
            queue.push_back({a, b});
          }
        }
      size_t guard = 0;
      while (!queue.empty())
      {
        if (++guard > 100 * (tris.size() + 10) * (tris.size() + 10))
          throw vcg::MissingPreconditionException("CoMEmbed: could not recover a curve segment; curves may touch or cross.");
        const int a = queue.front().first, b = queue.front().second;
        queue.pop_front();
        size_t t1, t2; int i1, i2;
        if (!FindEdge(a, b, t1, i1) || !FindEdge(b, a, t2, i2))
          throw vcg::MissingPreconditionException("CoMEmbed: a curve segment crosses the boundary of its face.");
        const int c = tris[t1][(i1 + 2) % 3], d = tris[t2][(i2 + 2) % 3];
        if (Orient(c, a, d) > eps && Orient(d, b, c) > eps)
        {
          tris[t1] = {{c, a, d}};
          tris[t2] = {{d, b, c}};
          if (Crosses(p, q, c, d)) queue.push_back({c, d});
        }
        else queue.push_back({a, b});
      }
      size_t ti; int i;
      if (!FindEdge(p, q, ti, i) && !FindEdge(q, p, ti, i))
        throw vcg::MissingPreconditionException("CoMEmbed: could not recover a curve segment; it passes through another curve point.");
      recovered.insert({std::min(p, q), std::max(p, q)});
    }
  };
};

} // end namespace tri
} // end namespace vcg

#endif // __VCGLIB_CURVE_ON_SURF_H
