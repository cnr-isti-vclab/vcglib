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
#include<vcg/complex/algorithms/mesh_assert.h>
#include<vcg/complex/algorithms/update/bounding.h>
#include<vcg/complex/algorithms/refine.h>
#include<vcg/complex/algorithms/vertex_interpolation.h>
#include<vcg/complex/algorithms/create/platonic.h>
#include<vcg/complex/algorithms/point_sampling.h>
#include <vcg/space/index/grid_static_ptr.h>
#include <vcg/space/index/kdtree/kdtree.h>
#include <vcg/math/histogram.h>
#include<vcg/space/distance3.h>
#include <vcg/complex/algorithms/attribute_seam.h>
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
    ScalarType maxSnapThr;         ///< The maximum distance allowed when snapping a polyline vertex onto a mesh vertex (currently unused)
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
      maxSnapThr         = m.bbox.Diag()/1000.0;
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
   * \brief Test if a barycentric coordinate is well snapped 
   * \param ip The barycentric coordinate to test (must sum to 1.0)
   * \return true if no snapping is needed (all coords are either 0, 1, or far from boundaries)
   * 
   * A barycentric coordinate is "well snapped" if each component is either:
   * - Exactly 0.0 or 1.0, OR
   * - Far enough from 0 and 1 (outside the barycentricSnapThr threshold)
   * 
   * This indicates the point doesn't need further snapping adjustment.
   */
  bool IsWellSnapped(const CoordType &ip)
  {
      for(int i=0;i<3;++i)
          if( (ip[i]< par.barycentricSnapThr         && ip[i]!= 0.0) ||
              (ip[i]> (1.0 - par.barycentricSnapThr) && ip[i]!= 1.0))
              return false;
      assert(ip[0]+ip[1]+ip[2] == 1.0);
      return true;
  }
  
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
  

  /**
   * \brief Find the minimum distance from a sample point to the polyline
   * \param samplePnt The point to measure distance from
   * \param edgeGrid Spatial acceleration grid for the polyline edges
   * \param poly The polyline (as edge mesh)
   * \param closestPoint Output: the closest point on the polyline
   * \return The minimum distance from samplePnt to the polyline
   */
  ScalarType MinDistOnEdge(CoordType samplePnt, EdgeGrid &edgeGrid, MeshType &poly, CoordType &closestPoint)
  {
      ScalarType polyDist;
      EdgeType *cep = vcg::tri::GetClosestEdgeBase(poly,edgeGrid,samplePnt,par.gridBailout,polyDist,closestPoint);        
      return polyDist;    
  }
  
  /**
   * \brief Find the closest point on a mesh edge to the polyline (static version)
   * \param v0 First vertex of the mesh edge
   * \param v1 Second vertex of the mesh edge
   * \param edgeGrid Spatial acceleration grid for the polyline edges
   * \param poly The polyline (as edge mesh)
   * \param closestPoint Output: the point on the edge [v0,v1] closest to the polyline
   * \return The minimum distance from the edge to the polyline
   * 
   * This samples the edge [v0,v1] uniformly and finds which sample is closest to the polyline.
   */
  static ScalarType MinDistOnEdge(VertexType *v0,VertexType *v1, EdgeGrid &edgeGrid, MeshType &poly, CoordType &closestPoint)
  {
    ScalarType minPolyDist = std::numeric_limits<ScalarType>::max();
    const ScalarType sampleNum = 50;
    const ScalarType maxDist = poly.bbox.Diag()/10.0;
    for(ScalarType k = 0;k<sampleNum+1;++k)
    {
      ScalarType polyDist;
      CoordType closestPPoly;
      CoordType samplePnt = (v0->P()*k +v1->P()*(sampleNum-k))/sampleNum;          
      
      EdgeType *cep = vcg::tri::GetClosestEdgeBase(poly,edgeGrid,samplePnt,maxDist,polyDist,closestPPoly);        
      
      if(polyDist < minPolyDist)
      {
        minPolyDist = polyDist;
        closestPoint = samplePnt;
//        closestPoint = closestPPoly;
      }
    }
    return minPolyDist;    
  }
  
  // ============================================================================
  // Attribute Extraction and Comparison (for Seam Processing)
  // ============================================================================
  
  /**
   * \brief Extract vertex attributes for seam processing
   * \param srcMesh Source mesh (unused but required by interface)
   * \param f The face containing the vertex
   * \param whichWedge Which vertex (0,1,2) of the face to extract
   * \param dstMesh Destination mesh (unused but required by interface)
   * \param v Output: vertex with copied attributes
   * 
   * This is a callback function used by the attribute_seam system.
   * It copies all per-vertex properties and uses the face color.
   * 
   * \note This is used when splitting the mesh along seams/polylines.
   */
  static inline void ExtractVertex(const MeshType & srcMesh, const FaceType & f, int whichWedge, const MeshType & dstMesh, VertexType & v)
  {
      (void)srcMesh;
      (void)dstMesh;
      // This is done to preserve every single perVertex property
      // perVextex Texture Coordinate is instead obtained from perWedge one.
      v.ImportData(*f.cV(whichWedge));
      v.C() = f.cC();
  }
  
  /**
   * \brief Compare two vertices for seam compatibility
   * \param m The mesh (unused but required by interface)
   * \param vA First vertex
   * \param vB Second vertex
   * \return true if vertices are compatible across a seam
   * 
   * This callback is used by the attribute_seam system to determine if two
   * vertices can be considered the same across a seam boundary.
   * Current implementation: Red and Blue colored vertices are considered incompatible.
   * 
   * \note This is part of the mesh cutting/seam processing infrastructure.
   */
  static inline bool CompareVertex(const MeshType & m, const VertexType & vA, const VertexType & vB)
  {
      (void)m;
      
      if(vA.C() == Color4b(Color4b::Red) && vB.C() == Color4b(Color4b::Blue) ) return false;
      if(vA.C() == Color4b(Color4b::Blue) && vB.C() == Color4b(Color4b::Red) ) return false;
      return true;      
  }
  
  // ============================================================================
  // Utility Functions
  // ============================================================================
  
  /**
   * \brief Compute quality-weighted linear interpolation between two vertices
   * \param v0 First vertex
   * \param v1 Second vertex
   * \return Interpolated position weighted by inverse quality values
   * 
   * Points with higher quality (larger absolute value) contribute less to the result.
   * This is useful for adaptive refinement based on error metrics stored in quality.
   */
  static CoordType QLerp(VertexType *v0, VertexType *v1)
  {
    
    ScalarType qSum = fabs(v0->Q())+fabs(v1->Q());      
    ScalarType w0 = (qSum - fabs(v0->Q()))/qSum;
    ScalarType w1 = (qSum - fabs(v1->Q()))/qSum;      
    return v0->P()*w0 + v1->P()*w1;      
  }
  
  
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
    int vertSnapCnt=0;
    int edgeSnapCnt=0;
    int borderCnt=0,midCnt=0,nonmanifCnt=0;
    for(VertexIterator vi=poly.vert.begin(); vi!=poly.vert.end();++vi)
    {
      CoordType ip;
      FaceType *f = GetClosestFaceIP(vi->cP(),ip);
      if(BarycentricSnap(ip))
      {
        if(ip[0]>0 && ip[1]>0) { vi->P() = f->P(0)*ip[0]+f->P(1)*ip[1]; edgeSnapCnt++; assert(ip[2]==0); }
        if(ip[0]>0 && ip[2]>0) { vi->P() = f->P(0)*ip[0]+f->P(2)*ip[2]; edgeSnapCnt++; assert(ip[1]==0); }
        if(ip[1]>0 && ip[2]>0) { vi->P() = f->P(1)*ip[1]+f->P(2)*ip[2]; edgeSnapCnt++; assert(ip[0]==0); }
        
        if(ip[0]==1.0) { vi->P() = f->P(0); vertSnapCnt++; assert(ip[1]==0 && ip[2]==0); }
        if(ip[1]==1.0) { vi->P() = f->P(1); vertSnapCnt++; assert(ip[0]==0 && ip[2]==0); }
        if(ip[2]==1.0) { vi->P() = f->P(2); vertSnapCnt++; assert(ip[0]==0 && ip[1]==0); }
      }
      else
      {
        int deg = edge::VEDegree<EdgeType>(&*vi);
        if (deg > 2) nonmanifCnt++;
        if (deg < 2) borderCnt++;
        if (deg== 2) midCnt++;
      }
    }
    Log("SnapPolyline %i vertices:  snapped %i onto vert and %i onto edges %i nonmanif, %i border, %i mid",
        poly.vn, vertSnapCnt, edgeSnapCnt, nonmanifCnt,borderCnt,midCnt);
    int dupCnt=tri::Clean<MeshType>::RemoveDuplicateVertex(poly);
    tri::Allocator<MeshType>::CompactEveryVector(poly);     
    if(dupCnt) Log("SnapPolyline: Removed %i Duplicated vertices",dupCnt);
    
    return vertSnapCnt==0 && edgeSnapCnt==0 && dupCnt==0;
  }
  
   void SelectBoundaryVertex(MeshType &poly)
   {
     tri::UpdateSelection<MeshType>::VertexClear(poly);
     tri::UpdateTopology<MeshType>::VertexEdge(poly);
     ForEachVertex(poly, [&](VertexType &v){
       if(edge::VEDegree<EdgeType>(&v)==1) v.SetS();
     });
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
      int cnt=tri::Clean<MeshType>::RemoveDuplicateVertex(poly);
      if(cnt)
          Log("SimplifyNullEdges: Removed %i Duplicated vertices",cnt);
  }
  
  void SimplifyMidEdge(MeshType &poly)
  {
   int startVn;
   int midEdgeCollapseCnt=0;
   tri::Allocator<MeshType>::CompactEveryVector(poly); 
   do
   {
    startVn = poly.vn;
    for(int ei =0; ei<poly.en; ++ei)
    {
      VertexType *v0=poly.edge[ei].V(0);
      VertexType *v1=poly.edge[ei].V(1);
      CoordType ip0,ip1;    
      FaceType *f0=GetClosestFaceIP(v0->P(),ip0);
      FaceType *f1=GetClosestFaceIP(v1->P(),ip1);
      
      bool snap0=BarycentricSnap(ip0);
      bool snap1=BarycentricSnap(ip1);
      int e0i,e1i;
      bool e0 = IsSnappedEdge(ip0,e0i);
      bool e1 = IsSnappedEdge(ip1,e1i);
      if(e0 && e1)
        if( (          f0 == f1           &&          e0i == e1i) || 
            (          f0 == f1->FFp(e1i) &&          e0i == f1->FFi(e1i)) || 
            (f0->FFp(e0i) == f1           && f0->FFi(e0i) == e1i) || 
            (f0->FFp(e0i) == f1->FFp(e1i) && f0->FFi(e0i) == f1->FFi(e1i)) ) 
        {
          CoordType newp = (v0->P()+v1->P())/2.0;
          v0->P()=newp;
          v1->P()=newp;
          midEdgeCollapseCnt++;
        }
    }
    tri::Clean<MeshType>::RemoveDuplicateVertex(poly);
    tri::Allocator<MeshType>::CompactEveryVector(poly);     
//    printf("SimplifyMidEdge %5i -> %5i %i mid %i ve \n",startVn,poly.vn,midEdgeCollapseCnt);
   } while(startVn>poly.vn);
  } 
  
  /**
   * @brief SimplifyMidFace remove all the vertices that in the mid of a face 
   * and between two of the points snapped onto the edges of the same face
   * @param poly
   * 
   * It assumes that the mesh has been snapped and refined by the BaseMesh
   * 
   */
  void SimplifyMidFace(MeshType &poly)
  {
   int startVn= poly.vn;;
   int midFaceCollapseCnt=0;
   int vertexEdgeCollapseCnt=0;
   int curVn;
   do
   {
    tri::Allocator<MeshType>::CompactEveryVector(poly); 
    curVn = poly.vn;
    UpdateTopology<MeshType>::VertexEdge(poly);
    for(int i =0; i<poly.vn;++i)
    {
      std::vector<VertexPointer> starVecVp;
      edge::VVStarVE(&(poly.vert[i]),starVecVp);      
      if( (starVecVp.size()==2) )
      {
        CoordType ipP, ipN, ipI; 
        FacePointer fpP = GetClosestFaceIP(starVecVp[0]->P(),ipP);
        FacePointer fpN = GetClosestFaceIP(starVecVp[1]->P(),ipN);
        FacePointer fpI = GetClosestFaceIP(poly.vert[i].P(), ipI);
        
        bool snapP = (BarycentricSnap(ipP));
        bool snapN = (BarycentricSnap(ipN));
        bool snapI = (BarycentricSnap(ipI));
        VertexPointer vertexSnapP = 0;
        VertexPointer vertexSnapN = 0;
        VertexPointer vertexSnapI = 0;
        for(int j=0;j<3;++j)
        {
          if(ipP[j]==1.0) vertexSnapP=fpP->V(j);
          if(ipN[j]==1.0) vertexSnapN=fpN->V(j);
          if(ipI[j]==1.0) vertexSnapI=fpI->V(j);
        }
        
        bool collapseFlag=false;
        
        if((!snapI && snapP && snapN) ||              // First case a vertex that is not snapped between two snapped vertexes 
           (!snapI && !snapP && fpI==fpP) || // Or a two vertex not snapped but on the same face
           (!snapI && !snapN && fpI==fpN) )
        {
          collapseFlag=true;
          midFaceCollapseCnt++;
        } 
        
        else  // case 2) a vertex snap and edge snap we have to check that the edge do not share the same vertex of the vertex snap
          if(snapI && snapP && snapN && vertexSnapI==0 && (vertexSnapP!=0 || vertexSnapN!=0) )
          {
            for(int j=0;j<3;++j) {
              if(ipI[j]!=0 && (fpI->V(j)==vertexSnapP || fpI->V(j)==vertexSnapN)) {
                collapseFlag=true;                                          
                vertexEdgeCollapseCnt++;
              }
            }
          }            
        
        if(collapseFlag)  
          edge::VEEdgeCollapse(poly,&(poly.vert[i]));
      }
    }  
   } while(curVn>poly.vn);
   Log("SimplifyMidFace %5i -> %5i %i mid %i ve",startVn,poly.vn,midFaceCollapseCnt,vertexEdgeCollapseCnt);
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
   * \brief Snap barycentric coordinates to 0 or 1 if within threshold
   * \param ip Input/Output: barycentric coordinates (must sum to 1.0)
   * \return true if the point was snapped to a vertex or edge (at least one coord became 0)
   * 
   * **This is one of the MOST IMPORTANT functions in the class** - it's used throughout!
   * 
   * Given barycentric coordinates of a point in a triangle, this function decides 
   * whether it should be "snapped" to a vertex or edge based on the 
   * `par.barycentricSnapThr` threshold.
   * 
   * **Algorithm:**
   * 1. If any coordinate is within `barycentricSnapThr` of 0, snap it to 0
   * 2. If any coordinate is within `barycentricSnapThr` of 1, snap it to 1
   * 3. Renormalize to ensure sum = 1.0
   * 4. If sum is still not exactly 1.0 (due to floating point), adjust the non-snapped coordinate
   * 
   * **Snapping Cases:**
   * - One coord = 1.0, others = 0 → Snapped to a vertex
   * - One coord = 0, others > 0 → Snapped to an edge
   * - All coords > 0 and < 1 → Interior point, NOT snapped
   * 
   * **Return Value:**
   * - `true`: Point is on a vertex or edge (at least one coordinate is 0)
   * - `false`: Point is in the interior of the triangle
   * 
   * \note This function MODIFIES the input coordinates in-place!
   * \note The threshold `par.barycentricSnapThr` (default 0.05) controls snapping sensitivity
   * 
   * \warning Side effect: modifies ip parameter! Consider renaming to BarycentricSnapInPlace()
   * 
   * \sa IsWellSnapped, IsSnappedVertex, IsSnappedEdge
   */
  bool BarycentricSnap(CoordType &ip)
  {
    for(int i=0;i<3;++i)
    {
      if(ip[i] <= par.barycentricSnapThr) ip[i]=0;
      if(ip[i] >= 1.0-par.barycentricSnapThr) ip[i]=1;
    }
    ScalarType sum = ip[0]+ip[1]+ip[2];
    
    for(int i=0;i<3;++i) 
      if(ip[i]!=1.0) ip[i]/=sum;
    
    sum = ip[0]+ip[1]+ip[2];
    
    if(sum!=1.0){
        for(int i=0;i<3;++i)
            if(ip[i]>0.0 && ip[i]<1.0) // if it is non snapped
                ip[i]=1.0-(ip[(i+1)%3]+ip[(i+2)%3]);
    }
    
    sum = ip[0]+ip[1]+ip[2];     
    assert(sum ==1.0);
    assert(IsWellSnapped(ip));
    if(ip[0]==0 || ip[1]==0 || ip[2]==0) return true;
    return false;
  }
  
  
  /**
   * @brief TestSplitSegWithMesh  Given a poly segment decide if it should be split along elements of base mesh. 
   * @param v0
   * @param v1
   * @param splitPt
   * @return true if it should be split
   * 
   * We make a few samples onto the edge and if some of them snaps onto a an edge we use it.
   * In case there are more than one candidate we choose the sample closeset to its snapping point.
   * We explicitly avoid snapping twice on the same edge by checking the starting and ending edges.
   * 
   * Two cases:
   * - poly edge pass near a vertex of the mesh
   * - poly edge cross one or more edges
   * 
   * Note that we have to check the case where 
   */
  bool TestSplitSegWithMesh(VertexType *v0, VertexType *v1, CoordType &splitPt)
  {
    Segment3Type segPoly(v0->P(),v1->P());
    const ScalarType sampleNum = 40;    
    CoordType ip0,ip1;
    
    FaceType *f0=GetClosestFaceIP(v0->P(),ip0);
    FaceType *f1=GetClosestFaceIP(v1->P(),ip1);
    if(f0==f1) return false;
    
    bool snap0=false,snap1=false; // true if the segment start/end on a edge/vert
    
    Segment3Type seg0; // The two segments to be avoided 
    Segment3Type seg1; // from which the current poly segment can start
    VertexPointer vertexSnap0 = 0;
    VertexPointer vertexSnap1 = 0;
    if(BarycentricSnap(ip0)) { 
      snap0=true; 
      for(int i=0;i<3;++i) {
        if(ip0[i]==1.0) vertexSnap0=f0->V(i);
        if(ip0[i]==0.0) seg0=Segment3Type(f0->P1(i),f0->P2(i)); 
      }        
    } 
    if(BarycentricSnap(ip1)) { 
      snap1=true; 
      for(int i=0;i<3;++i){
        if(ip1[i]==1.0) vertexSnap1=f1->V(i);
        if(ip1[i]==0.0) seg1=Segment3Type(f1->P1(i),f1->P2(i)); 
      }        
    } 
    
    CoordType bestSplitPt(0,0,0);
    ScalarType bestDist = std::numeric_limits<ScalarType>::max();
    for(ScalarType k = 1;k<sampleNum;++k)
    {
      CoordType samplePnt = segPoly.Lerp(k/sampleNum);    
      CoordType ip;
      FaceType *f=GetClosestFaceIP(samplePnt,ip);
//      BarycentricEdgeSnap(ip);
      if(BarycentricSnap(ip))
      {
        VertexPointer vertexSnapI = 0;        
        for(int i=0;i<3;++i)
          if(ip[i]==1.0) vertexSnapI=f->V(i);
        CoordType closestPt = f->P(0)*ip[0]+f->P(1)*ip[1]+f->P(2)*ip[2];
        if(Distance(samplePnt,closestPt) < bestDist )  
        {
          ScalarType dist0=std::numeric_limits<ScalarType>::max();
          ScalarType dist1=std::numeric_limits<ScalarType>::max();
          CoordType closestSegPt;
          if(snap0) SegmentPointDistance(seg0,closestPt,closestSegPt,dist0);
          if(snap1) SegmentPointDistance(seg1,closestPt,closestSegPt,dist1);
          if( (!vertexSnapI && (dist0 > par.surfDistThr/1000 && dist1>par.surfDistThr/1000) ) ||
              ( vertexSnapI!=vertexSnap0 && vertexSnapI!=vertexSnap1)  )
          {
            bestDist = Distance(samplePnt,closestPt);
            bestSplitPt = closestPt;            
          }
        }      
      }
    }
    if(bestDist < par.surfDistThr*100)
    {
      splitPt = bestSplitPt;
      return true;
    }
    
    return false;
  }
  /**
   * @brief SnappedOnSameFace Return true if the two points are snapped to a common face;
   * @param f0
   * @param i0
   * @param f1
   * @param i0
   * @return 
   * 
   * Require FFAdj. se assume that both SNAPPED. Three cases:
   * - Edge Edge - true iff the two edges belongs to a common face. 
   * - Vert Edge - true iff there is one of the two snapped edge faces has the vert as non-edge face;  
   * - Vert Vert 
   * 
   */
  bool SnappedOnSameFace(FacePointer f0, CoordType i0, FacePointer f1, CoordType i1)
  {
   if(f0==f1) return true;
   int e0,e1;
   int v0,v1;
   bool e0Snap = IsSnappedEdge(i0,e0);
   bool e1Snap = IsSnappedEdge(i1,e1);
   bool v0Snap = IsSnappedVertex(i0,v0);
   bool v1Snap = IsSnappedVertex(i1,v1);
   FacePointer f0p=0; int e0p=-1;  // When Edge snap the other face and the index of the snapped edge on the other face
   FacePointer f1p=0; int e1p=-1;
   assert((e0Snap != v0Snap) && (e1Snap != v1Snap));
   // For EdgeSnap compute the 'other' face stuff 
   if(e0Snap){
     f0p = f0->FFp(e0); e0p=f0->FFi(e0); assert(f0p->FFp(e0p)==f0);
   }
   if(e1Snap){
     f1p = f1->FFp(e1); e1p=f1->FFi(e1); assert(f1p->FFp(e1p)==f1);
   }
   
   if(e0Snap && e1Snap) {
    if(f0==f1p || f0p==f1p || f0p==f1 || f0==f1) return true;
   }
   
   if(e0Snap && v1Snap)  {
     assert(v1>=0 && v1<3 && v0==-1 && e1==-1);
     if(f0->V2(e0)  ==f1->V(v1)) return true;
     if(f0p->V2(e0p)==f1->V(v1)) return true;
   }
     
   if(e1Snap && v0Snap)  {
     assert(v0>=0 && v0<3 && v1==-1 && e0==-1);
     if(f1->V2(e1)  ==f0->V(v0)) return true;
     if(f1p->V2(e1p)==f0->V(v0)) return true;
   }
     
   if(v1Snap && v0Snap)  {
     PosType startPos(f0,f0->V(v0));
     PosType curPos=startPos;
     do
     {
       assert(curPos.V()==f0->V(v0));
       if(curPos.VFlip()==f1->V(v1)) return true;
       curPos.FlipE();
       curPos.FlipF();       
     }
     while(curPos!=startPos);   
   }   
   return false;    
  }
  
  /**
   * @brief TestSplitSegWithMesh  Given a poly segment decide if it should be split along elements of base mesh. 
   * @param v0
   * @param v1
   * @param splitPt
   * @return true if it should be split
   * 
   * We make a few samples onto the edge and if some of them snaps onto a an edge we use it.
   * In case there are more than one candidate we choose the sample closeset to its snapping point.
   * We explicitly avoid snapping twice on the same edge by checking the starting and ending edges.
   * 
   * Two cases:
   * - poly edge pass near a vertex of the mesh
   * - poly edge cross one or more edges
   * 
   * Note that we have to check the case where 
   */
  bool TestSplitSegWithMeshAdapt(VertexType *v0, VertexType *v1, CoordType &splitPt)
  {
    splitPt=(v0->P()+v1->P())/2.0;
      
    CoordType ip0,ip1,ipm;    
    FaceType *f0=GetClosestFaceIP(v0->P(),ip0);
    FaceType *f1=GetClosestFaceIP(v1->P(),ip1);
    FaceType *fm=GetClosestFaceIP(splitPt,ipm);
    
    if(f0==f1) return false;
    
    bool snap0=BarycentricSnap(ip0);
    bool snap1=BarycentricSnap(ip1);
    bool snapm=BarycentricSnap(ipm);
    
    splitPt = fm->P(0)*ipm[0]+fm->P(1)*ipm[1]+fm->P(2)*ipm[2];
    
    if(!snap0 && !snap1) {
      assert(f0!=f1);
      return true;
    }
    if(snap0 && snap1) 
    {
      if(SnappedOnSameFace(f0,ip0,f1,ip1)) 
        return false;            
    }
    
    if(snap0) {
      int e0,v0;
      if (IsSnappedEdge(ip0,e0)) {
        if(f0->FFp(e0) == f1) return false;
      }
      if(IsSnappedVertex(ip0,v0)) {
        for(int i=0;i<3;++i) 
          if(f1->V(i)==f0->V(v0)) return false;
      }
    }
    if(snap1) {
      int e1,v1;
      if (IsSnappedEdge(ip1,e1)) {
        if(f1->FFp(e1) == f0) return false;
      }
      if(IsSnappedVertex(ip1,v1)) {
        for(int i=0;i<3;++i) 
          if(f0->V(i)==f1->V(v1)) return false;
      }
    }
    
    return true;
  }
  
  
  bool TestSplitSegWithMeshAdaptOld(VertexType *v0, VertexType *v1, CoordType &splitPt)
  {
    Segment3Type segPoly(v0->P(),v1->P());
    const ScalarType sampleNum = 40;    
    CoordType ip0,ip1;    
    FaceType *f0=GetClosestFaceIP(v0->P(),ip0);
    FaceType *f1=GetClosestFaceIP(v1->P(),ip1);
    if(f0==f1) return false;
    
    bool snap0=BarycentricSnap(ip0);
    bool snap1=BarycentricSnap(ip1);
    
    if(!snap0 && !snap1) {
      assert(f0!=f1);
      splitPt=(v0->P()+v1->P())/2.0;
      return true;
    }
    if(snap0 && snap1) 
    {
      if(SnappedOnSameFace(f0,ip0,f1,ip1)) 
        return false;      
    }
    
    if(snap0) {
      int e0,v0;
      if (IsSnappedEdge(ip0,e0)) {
        if(f0->FFp(e0) == f1) return false;
      }
      if(IsSnappedVertex(ip0,v0)) {
        for(int i=0;i<3;++i) 
          if(f1->V(i)==f0->V(v0)) return false;
      }
    }
    splitPt=(v0->P()+v1->P())/2.0;
    return true;
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
   * @brief RefineCurveByBaseMesh
   * @param poly
   */
  
  void RefineCurveByBaseMesh(MeshType &poly)
  {
    tri::Allocator<MeshType>::CompactEveryVector(poly);    
    tri::UpdateTopology<MeshType>::VertexEdge(poly); // the edge splits below walk VE adjacency
    std::vector<int> edgeToRefineVec;
    for(int i=0; i<poly.en;++i) 
      edgeToRefineVec.push_back(i);
    int startEn=poly.en;  
    int iterCnt=0;
    while (!edgeToRefineVec.empty() && iterCnt<100) {
      iterCnt++;
      std::vector<int> edgeToRefineVecNext;
      for(int i=0; i<edgeToRefineVec.size();++i)
      {
        EdgeType &e = poly.edge[edgeToRefineVec[i]];
        CoordType splitPt;
        if(TestSplitSegWithMeshAdapt(e.V(0),e.V(1),splitPt))  
        {
          edge::VEEdgeSplit(poly, &e, splitPt); 
          edgeToRefineVecNext.push_back(edgeToRefineVec[i]);
          edgeToRefineVecNext.push_back(poly.en-1);
        } 
      }
      tri::Allocator<MeshType>::CompactEveryVector(poly);
      swap(edgeToRefineVecNext,edgeToRefineVec);
      Progress(iterCnt, "RefineCurveByBaseMesh %i en -> %i en",startEn,poly.en); // at most 100 iterations
    }
//
    SimplifyNullEdges(poly);
    SimplifyMidFace(poly);
    SimplifyMidEdge(poly);
    SnapPolyline(poly);    
    Log("RefineCurveByBaseMesh %i en -> %i en",startEn,poly.en);
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
      
      tri::UpdateTopology<MeshType>::TestVertexEdge(poly);
      RefineCurveByDistance(poly);      
      tri::UpdateTopology<MeshType>::TestVertexEdge(poly);
      Simplify(poly);
      tri::UpdateTopology<MeshType>::TestVertexEdge(poly);
      int dupVertNum = Clean<MeshType>::RemoveDuplicateVertex(poly);
      if(dupVertNum) {
        tri::Allocator<MeshType>::CompactEveryVector(poly);
        tri::UpdateTopology<MeshType>::VertexEdge(poly);
      }
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

		if(com.BarycentricSnap(ip0) && com.BarycentricSnap(ip1))
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
   * Preconditions: RefineCurveByBaseMesh() has been called, so every polyline segment lies
   * inside one face or along one edge of the base mesh; and the polyline has no contacts:
   * its vertices are distinct and its segments do not cross each other.
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
   *         along one edge, or if a segment cannot be recovered (curves touching or crossing)
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
      com.BarycentricSnap(l.ip);
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
    //    polyline vertex at exactly the same place.
    std::vector<size_t> meshV(poly.vert.size());
    std::map<std::pair<EdgeKey, ScalarType>, size_t> edgePointVert;
    size_t newVertNum = 0;
    const size_t firstNewVert = m.vert.size();
    for (size_t ci = 0; ci < poly.vert.size(); ++ci)
    {
      const Loc &l = loc[ci];
      if (l.kind == OnVertex) meshV[ci] = l.vert;
      else if (l.kind == OnEdge) {
        auto ins = edgePointVert.insert({{l.edge, l.t}, firstNewVert + newVertNum});
        if (ins.second) ++newVertNum;
        meshV[ci] = ins.first->second;
      }
      else meshV[ci] = firstNewVert + newVertNum++;
    }
    tri::Allocator<MeshType>::AddVertices(m, newVertNum);
    for (size_t ci = 0; ci < poly.vert.size(); ++ci)
    {
      const Loc &l = loc[ci];
      VertexType &nv = m.vert[meshV[ci]];
      if (l.kind == OnEdge) {
        const VertexType &va = m.vert[l.edge.first], &vb = m.vert[l.edge.second];
        nv.P() = va.cP() * (1 - l.t) + vb.cP() * l.t;
        VertexInterpolator<MeshType>::Lerp(m, nv, va, vb, l.t);
      } else if (l.kind == InFace) {
        const FaceType &f = m.face[l.face[0]];
        nv.P() = f.cP(0) * l.ip[0] + f.cP(1) * l.ip[1] + f.cP(2) * l.ip[2];
        VertexInterpolator<MeshType>::Barycentric(m, nv, *f.cV(0), *f.cV(1), *f.cV(2), l.ip);
      }
      poly.vert[ci].P() = m.vert[meshV[ci]].cP();
    }

    // 3. What each face has to take in: points on its edges, points inside it, segments.
    std::map<EdgeKey, std::vector<std::pair<ScalarType, size_t>>> edgePoints;  // (t from first, vertex)
    for (const auto &ep : edgePointVert) edgePoints[ep.first.first].push_back({ep.first.second, ep.second});
    std::map<size_t, std::vector<std::pair<size_t, CoordType>>> facePoints;   // (vertex, barycentric)
    std::map<size_t, std::vector<EdgeKey>> faceSegments;
    std::set<EdgeKey> curveEdges;
    std::map<size_t, std::vector<size_t>> vertFaces;  // faces around the polyline vertices on mesh vertices
    for (size_t ci = 0; ci < poly.vert.size(); ++ci)
    {
      if (loc[ci].kind == InFace) facePoints[loc[ci].face[0]].push_back({meshV[ci], loc[ci].ip});
      if (loc[ci].kind == OnVertex) vertFaces[loc[ci].vert];
    }
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
    for (const auto &e : poly.edge) if (!e.IsD())
    {
      const size_t c0 = tri::Index(poly, e.cV(0)), c1 = tri::Index(poly, e.cV(1));
      if (meshV[c0] == meshV[c1]) continue;  // a zero-length segment
      curveEdges.insert(key(meshV[c0], meshV[c1]));
      // The segment lies in the faces both its ends belong to: one face for a segment
      // across it, the two faces of an edge for a segment along it.
      std::vector<size_t> f0 = candidateFaces(c0), f1 = candidateFaces(c1), common;
      std::sort(f0.begin(), f0.end()); std::sort(f1.begin(), f1.end());
      std::set_intersection(f0.begin(), f0.end(), f1.begin(), f1.end(), std::back_inserter(common));
      if (common.empty())
        throw vcg::MissingPreconditionException("CoMEmbed: a curve segment does not lie in one face; call RefineCurveByBaseMesh first.");
      if (common.size() == 1) faceSegments[common[0]].push_back(key(meshV[c0], meshV[c1]));
    }

    std::set<size_t> touched;
    for (const auto &fp : facePoints) touched.insert(fp.first);
    for (const auto &fs : faceSegments) touched.insert(fs.first);
    for (size_t fi = 0; fi < m.face.size(); ++fi) if (!m.face[fi].IsD())
      for (int i = 0; i < 3; ++i)
        if (edgePoints.count(key(tri::Index(m, m.face[fi].V(i)), tri::Index(m, m.face[fi].V1(i))))) touched.insert(fi);

    // 4. Rebuild every touched face on its own.
    std::vector<std::pair<size_t, LocalTriangulation>> rebuilt;
    size_t newFaceNum = 0, done = 0;
    for (size_t fi : touched)
    {
      com.Progress(int(100 * done++ / touched.size()), "SplitMeshWithPolyline: rebuilding face %lu of %lu", (unsigned long)done, (unsigned long)touched.size());
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
      for (const EdgeKey &s : faceSegments[fi]) lt.Recover(lt.Local(s.first), lt.Local(s.second));
      newFaceNum += lt.tris.size() - 1;
      rebuilt.push_back({fi, std::move(lt)});
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

    // 6. Every polyline edge is now a mesh edge: select it on both sides.
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
    /// crosses is queued again. Recovered segments never cross, so none of them is flipped.
    void Recover(int p, int q)
    {
      std::deque<std::pair<int, int>> queue;
      for (const auto &t : tris)
        for (int i = 0; i < 3; ++i)
        {
          const int a = t[i], b = t[(i + 1) % 3];
          if (a < b && Crosses(p, q, a, b)) queue.push_back({a, b});
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
    }
  };
};

} // end namespace tri
} // end namespace vcg

#endif // __VCGLIB_CURVE_ON_SURF_H
