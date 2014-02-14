/****************************************************************************
* VCGLib                                                            o o     *
* Visual and Computer Graphics Library                            o     o   *
*                                                                _   O  _   *
* Copyright(C) 2004                                                \/)\/    *
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
#ifndef POLYGON_POLYCOORD_COLLAPSE_H
#define POLYGON_POLYCOORD_COLLAPSE_H

#include <vector>
#include <list>
#include <set>
#include <map>
#include <queue>
#include <utility>
#include <vcg/complex/complex.h>
#include <vcg/simplex/face/jumping_pos.h>

namespace vcg {
namespace tri {
/** \addtogroup trimesh */

/**
* @brief The PolycoordCollapse class provides methods to semplify a quad mesh, by collapsing the polycoords.
*
* This class is an implementation of a method very similar to that for mesh semplification proposed
* by Daniels et al. in "Quadrilateral mesh simplification", see http://www.cs.utah.edu/~jdaniels/research/asia2008_qms.htm
* The main function is PolycoordCollapse::CollapsePolycoord() which deletes all the quadrilateral faces in a polycoord.
* The polycoords that can be collapsed in this case are those forming a closed loop (a ring) or that start and end to
* mesh borders. A way to preserve the structure of the singularities is also provided.
* The convenient method PolycoordCollapse::CollapseAllPolycoords() finds and collapses all the polycoords on a mesh.
* The input mesh should be polygonal, i.e. it should have the vcg::face::PolyInfo component. Even though a generic
* triangle mesh can be given, actually the class does not perform any collapsing operation since it sees only triangles,
* in fact it does not consider faux edges.
*/
template < typename PolyMeshType >
class PolycoordCollapse {
public:
  typedef typename PolyMeshType::FaceType   FaceType;
  typedef typename PolyMeshType::VertexType VertexType;
  typedef typename PolyMeshType::CoordType  CoordType;

  /**
  * @brief The PC_ResultCode enum codifies the result type of a polycoord collapse operation.
  */
  enum PC_ResultCode {
    PC_SUCCESS = 0,
    PC_NOTMANIF = 1,
    PC_NOTQUAD = 2,
    PC_NOLINKCOND = 4,
    PC_SINGBOTH = 8,
    PC_SELFINTERSECT = 16,
    PC_VOID = 32
  };

  /**
  * @brief The PC_Coord struct identifies a coord of a polycoord passing through a quad.
  */
  struct PC_Coord {
    unsigned long mark;
    PC_ResultCode q;
    PC_Coord * prev;
    PC_Coord * next;
    PC_Coord() : mark(std::numeric_limits<unsigned long>::max()), q(PC_VOID), prev(NULL), next(NULL) { }
    inline void Reset() {
      mark = std::numeric_limits<unsigned long>::max();
      q = PC_VOID;
      prev = next = NULL;
    }
  };

  /**
  * @brief The PC_Coords class gives efficient access to each coord (relative to a face).
  */
  class PC_Coords {
  public:
    /**
     * @brief PC_Coords constructor.
     * @note Since each face corresponds to two coords, the actual size of the vector of coords is 2*mesh.face.size().
     * @param mesh
     */
    PC_Coords (const PolyMeshType &mesh) : _coords(2*mesh.face.size()), _currentCoord(NULL) {
      Reset(mesh);
    }

    /**
     * @brief ResetMarks
     */
    void ResetMarks() {
      typename std::vector<PC_Coord>::iterator it = _coords.begin();
      for (; it != _coords.end(); it++)
        (*it).mark = std::numeric_limits<unsigned long>::max();
    }

    /**
     * @brief Reset rearrages the container.
     * @note Since each face corresponds to two coords, the actual size of the vector of coords is 2*mesh.face.size().
     * @param mesh
     */
    void Reset(const PolyMeshType &mesh) {
      _coords.resize(2*mesh.face.size());
      for (size_t j = 0; j < _coords.size(); j++)
        _coords[j].Reset();
      _currentCoord = NULL;

      PC_Coord *coord = NULL;
      long long j = 0;
      for (size_t i = 0; i < _coords.size(); i++) {
        // set the prev
        coord = NULL;
        if ((long long)i-1 >= 0) {
          coord = &_coords[i-1];
          if (vcg::tri::HasPerFaceFlags(mesh)) {
            j = i-1;
            while (j >= 0 && mesh.face[j/2].IsD())
              j--;
            if (j >= 0)
              coord = &_coords[j];
            else
              coord = NULL;
          }
        }
        _coords[i].prev = coord;

        // set the next
        coord = NULL;
        if (i+1 < _coords.size()) {
          coord = &_coords[i+1];
          if (vcg::tri::HasPerFaceFlags(mesh)) {
            j = i+1;
            while (j < (long long)_coords.size() && mesh.face[j/2].IsD())
              j++;
            if (j < (long long)_coords.size())
              coord = &_coords[j];
            else
              coord = NULL;
          }
        }
        _coords[i].next = coord;
      }
      if (mesh.face.size() > 0) {
        // set the current coord (first - not deleted - face)
        _currentCoord = &_coords[0];
        if (vcg::tri::HasPerFaceFlags(mesh) && mesh.face[0].IsD())
          _currentCoord = _currentCoord->next;
      }
    }

    /**
     * @brief operator [], given a face index and an offset, it returns (a reference to) its corresponding PC_Coord.
     * @param face_edge A std::pair<size_t, unsigned char>(face_index, offset). The offset should be 0 or 1.
     * @return A reference to the corresponding PC_Coord.
     */
    inline PC_Coord & operator[] (const std::pair<size_t, unsigned char> &face_edge) {
      assert(face_edge.first >= 0 && 2*face_edge.first+face_edge.second < _coords.size());
      return _coords[2*face_edge.first + face_edge.second];
    }
    /**
     * @brief operator [], given a face index and an offset, it returns (a const reference to) its corresponding PC_Coord.
     * @param face_edge A std::pair<size_t, unsigned char>(face_index, offset). The offset should be 0 or 1.
     * @return A reference to the corresponding PC_Coord.
     */
    inline const PC_Coord & operator[] (const std::pair<size_t, unsigned char> &face_edge) const {
      assert(face_edge.first >= 0 && 2*face_edge.first+face_edge.second < _coords.size());
      return _coords[2*face_edge.first + face_edge.second];
    }

    /**
     * @brief operator [], given a coord, it returns its corresponding face index and edge.
     * @param coord The coord pointer.
     * @return A std::pair <size_t, unsigned char>(face_index, offset) with offset being 0 or 1.
     */
    inline std::pair<size_t, unsigned char> operator[] (PC_Coord const * const coord) {
      assert(coord >= &_coords[0] && coord < &_coords[0]+_coords.size());
      return std::pair<size_t, unsigned char>((coord - &_coords[0])/2, (coord - &_coords[0])%2);
    }

    /**
     * @brief UpdateCoord updates the coord information and links.
     * @param coord The coord to update.
     * @param mark The mark of the polycoord.
     * @param resultCode The code for the type of the polycoord.
     */
    inline void UpdateCoord (PC_Coord &coord, const unsigned long mark, const PC_ResultCode resultCode) {
      // update prev and next
      if (coord.q == PC_VOID) {
        if (coord.prev != NULL && &coord != _currentCoord)
          coord.prev->next = coord.next;
        if (coord.next != NULL && &coord != _currentCoord)
          coord.next->prev = coord.prev;
        coord.mark = mark;
      }
      coord.q = resultCode;
    }

    /**
     * @brief Next, if it's not at the end, it goes to the next coord.
     */
    inline void Next () {
      if (_currentCoord != NULL)
        _currentCoord = _currentCoord->next;
    }

    /**
     * @brief GetCurrent returns the current FaceType pointer and edge.
     * @param face_edge A std::pair where to store the FaceType pointer and the edge index.
     */
    inline void GetCurrent (std::pair<size_t, unsigned char> &face_edge) {
      if (_currentCoord != NULL) {
        face_edge.first = (_currentCoord - &_coords[0])/2;
        face_edge.second = (_currentCoord - &_coords[0])%2;
      } else {
        face_edge.first = std::numeric_limits<size_t>::max();
        face_edge.second = 0;
      }
    }

    /**
     * @brief End says if an end has been reached.
     * @return true if an end has been reached, false otherwise.
     */
    inline bool End () {
      return _currentCoord == NULL;
    }

  private:
    std::vector<PC_Coord>   _coords;
    PC_Coord                *_currentCoord;
  };

  /**
   * @brief The LinkCondition class provides a tool to check if a polycoord satisfies the link conditions.
   */
  class LinkConditions {
  private:
    struct LCEdge;
    struct LCVertex;
    typedef std::set<LCVertex *> LCVertexStar;  // define the star of a vertex
    typedef std::set<LCEdge *> LCEdgeStar;      // define the set of edges whose star involves a vertex

  public:
    /**
     * @brief LinkCondition constructor.
     * @param size The number of vertices of the mesh.
     */
    LinkConditions (const size_t size) : _lcVertices(size) { }

    /**
     * @brief Resize just resets the size of the container.
     * @param size
     */
    inline void Resize(const size_t size) {
      _lcVertices.resize(size);
    }

    /**
     * @brief CheckLinkConditions checks if collapsing the polycoord starting from startPos
     * satisfies the link conditions.
     * @warning The polycoord starts from startPos and ends to itself (if it's a loop) or to a border. In the latter case,
     * call this method starting from the opposite border of the strip of quads.
     * @param mesh The mesh for getting the vertex index.
     * @param startPos The starting position of the polycoord.
     * @return true if satisfied, false otherwise.
     */
    bool CheckLinkConditions (const PolyMeshType &mesh, const vcg::face::Pos<FaceType> &startPos) {
      assert(!startPos.IsNull());
      assert(mesh.vert.size() == _lcVertices.size());
      std::list<LCEdge> lcEdges;
      LCEdge *e = NULL;
      LCVertexStar intersection;

      // reset the stars
      LC_ResetStars(mesh, startPos);

      // compute the stars
      LC_computeStars(mesh, startPos, lcEdges);

      // for each edge e = (v1,v2)
      // if intersection( star(v1) , star(v2) ) == star(e)
      //      then collapse e
      // else
      //      return false (i.e. link conditions not satisfied)
      for (typename std::list<LCEdge>::iterator eIt = lcEdges.begin(); eIt != lcEdges.end(); eIt++) {
        e = &*eIt;
        // compute the intersetion
        SetIntersection(e->v1->star, e->v2->star, intersection);
        // if intersection( star(v1) , star(v2) ) != star(e) then return false
        if (intersection != e->star)
            return false;
        // else simulate the collapse
        LC_SimulateEdgeCollapse(*e);
      }
      // at this point all collapses are possible, thus return true
      return true;
    }

  private:
    /**
     * @brief SetIntersection computes the set intersection between two sets.
     * @param set1
     * @param set2
     * @param result The set resulting from the intersection.
     */
    static void SetIntersection (const LCVertexStar &set1, const LCVertexStar &set2, LCVertexStar &result) {
      typename LCVertexStar::const_iterator set1It = set1.begin();
      typename LCVertexStar::const_iterator set2It = set2.begin();
      result.clear();
      while (set1It != set1.end() && set2It != set2.end()) {
        if (*set1It < *set2It) ++set1It;
        else if (*set2It < *set1It) ++set2It;
        else {
          result.insert(*set1It);
          ++set1It;
          ++set2It;
        }
      }
    }

    /**
     * @brief LC_ResetStars resets the stars on a polycoord.
     * @param mesh The mesh for getting the vertex index.
     * @param startPos
     */
    void LC_ResetStars (const PolyMeshType &mesh, const vcg::face::Pos<FaceType> &startPos) {
      assert(!startPos.IsNull());
      assert(mesh.vert.size() == _lcVertices.size());
      vcg::face::Pos<FaceType> runPos = startPos;
      vcg::face::JumpingPos<FaceType> vStarPos;
      // reset the stars
      do {
        // reset the star of this edge endpoints
        _lcVertices[vcg::tri::Index(mesh, runPos.V())].edges.clear();
        _lcVertices[vcg::tri::Index(mesh, runPos.V())].star.clear();
        _lcVertices[vcg::tri::Index(mesh, runPos.VFlip())].edges.clear();
        _lcVertices[vcg::tri::Index(mesh, runPos.VFlip())].star.clear();
        // reset the stars of the vertices in the star of the second vertex
        runPos.FlipV();
        vStarPos.Set(runPos.F(), runPos.E(), runPos.V());
        do {
          vStarPos.FlipV();
          vStarPos.FlipE();
          while (vStarPos.V() != runPos.V()) {
            _lcVertices[vcg::tri::Index(mesh, vStarPos.V())].edges.clear();
            _lcVertices[vcg::tri::Index(mesh, vStarPos.V())].star.clear();
            vStarPos.FlipV();
            vStarPos.FlipE();
          }
          vStarPos.NextFE();
        } while (vStarPos != runPos);
        // reset the stars of the vertices in the star of the first vertex
        runPos.FlipV();
        vStarPos.Set(runPos.F(), runPos.E(), runPos.V());
        do {
          vStarPos.FlipV();
          vStarPos.FlipE();
          while (vStarPos.V() != runPos.V()) {
            _lcVertices[vcg::tri::Index(mesh, vStarPos.V())].edges.clear();
            _lcVertices[vcg::tri::Index(mesh, vStarPos.V())].star.clear();
            vStarPos.FlipV();
            vStarPos.FlipE();
          }
          vStarPos.NextFE();
        } while (vStarPos != runPos);
        // when arrive to a border, return
        if (runPos != startPos && runPos.IsBorder())
          break;
        // go on the next edge
        runPos.FlipE();
        runPos.FlipV();
        runPos.FlipE();
        runPos.FlipF();
      } while (runPos != startPos);
    }

    /**
     * @brief LC_computeStars computes the stars of edges and vertices of the polycoord from the starting pos
     * either to itself (if it's a loop) or to the border edge.
     * @param mesh The mesh for getting the vertex index.
     * @param startPos Starting position.
     * @param lcEdges List of edge stars.
     */
    void LC_computeStars (const PolyMeshType &mesh, const vcg::face::Pos<FaceType> &startPos, std::list<LCEdge> &lcEdges)
    {
      assert(!startPos.IsNull());
      assert(mesh.vert.size() == _lcVertices.size());
      LCEdge *lcedgeP = NULL;
      vcg::face::Pos<FaceType> runPos = startPos;
      vcg::face::JumpingPos<FaceType> vStarPos;
      vcg::face::Pos<FaceType> eStarPos;

      lcEdges.clear();
      /// compute the star of all the vertices and edges seen from the polycoord
      runPos = startPos;
      do {
        // create a lcedge
        lcEdges.push_back(LCEdge());
        lcedgeP = &lcEdges.back();
        // set lcvertices references
        lcedgeP->v1 = &_lcVertices[vcg::tri::Index(mesh, runPos.V())];
        lcedgeP->v2 = &_lcVertices[vcg::tri::Index(mesh, runPos.VFlip())];
        // add this edge to its vertices edge-stars
        lcedgeP->v1->edges.insert(lcedgeP);
        lcedgeP->v2->edges.insert(lcedgeP);
        // compute the star of this edge
        lcedgeP->star.insert(lcedgeP->v1);  // its endpoints, clearly
        lcedgeP->star.insert(lcedgeP->v2);  // its endpoints, clearly
        // navigate over the other vertices of this facet
        eStarPos = runPos;
        eStarPos.FlipE();
        eStarPos.FlipV();
        while (eStarPos.V() != runPos.VFlip()) {
          // add current vertex to the star of this edge
          lcedgeP->star.insert(&_lcVertices[vcg::tri::Index(mesh, eStarPos.V())]);
          // add this edge to the edge-star of the current vertex
          _lcVertices[vcg::tri::Index(mesh, eStarPos.V())].edges.insert(lcedgeP);
          // go on
          eStarPos.FlipE();
          eStarPos.FlipV();
        }
        // go on the opposite facet
        if (!runPos.IsBorder()) {
          eStarPos = runPos;
          eStarPos.FlipF();
          eStarPos.FlipE();
          eStarPos.FlipV();
          while (eStarPos.V() != runPos.VFlip()) {
            // add current vertex to the star of this edge
            lcedgeP->star.insert(&_lcVertices[vcg::tri::Index(mesh, eStarPos.V())]);
            // add this edge to the edge-star of the current vertex
            _lcVertices[vcg::tri::Index(mesh, eStarPos.V())].edges.insert(lcedgeP);
            // go on
            eStarPos.FlipE();
            eStarPos.FlipV();
          }
        }

        // compute the star of vertex v2
        runPos.FlipV();
        vStarPos.Set(runPos.F(), runPos.E(), runPos.V());
        // v2 is in its star
        _lcVertices[vcg::tri::Index(mesh, vStarPos.V())].star.insert(&_lcVertices[vcg::tri::Index(mesh, vStarPos.V())]);
        do {
          vStarPos.FlipV();
          vStarPos.FlipE();
          while (vStarPos.V() != runPos.V()) {
            // add the current vertex to the v2 star
            _lcVertices[vcg::tri::Index(mesh, runPos.V())].star.insert(&_lcVertices[vcg::tri::Index(mesh, vStarPos.V())]);
            // add v2 to the star of the current vertex
            _lcVertices[vcg::tri::Index(mesh, vStarPos.V())].star.insert(&_lcVertices[vcg::tri::Index(mesh, runPos.V())]);
            vStarPos.FlipV();
            vStarPos.FlipE();
          }
          vStarPos.NextFE();
        } while (vStarPos != runPos);

        // compute the star of vertex v1
        runPos.FlipV();
        vStarPos.Set(runPos.F(), runPos.E(), runPos.V());
        // v1 is in its star
        _lcVertices[vcg::tri::Index(mesh, vStarPos.V())].star.insert(&_lcVertices[vcg::tri::Index(mesh, vStarPos.V())]);
        do {
          vStarPos.FlipV();
          vStarPos.FlipE();
          while (vStarPos.V() != runPos.V()) {
            // add the current vertex to the v2 star
            _lcVertices[vcg::tri::Index(mesh, runPos.V())].star.insert(&_lcVertices[vcg::tri::Index(mesh, vStarPos.V())]);
            // add v2 to the star of the current vertex
            _lcVertices[vcg::tri::Index(mesh, vStarPos.V())].star.insert(&_lcVertices[vcg::tri::Index(mesh, runPos.V())]);
            vStarPos.FlipV();
            vStarPos.FlipE();
          }
          vStarPos.NextFE();
        } while (vStarPos != runPos);

        // when arrive to a border, stop
        if (runPos != startPos && runPos.IsBorder())
          break;

        // go on the next edge
        runPos.FlipE();
        runPos.FlipV();
        runPos.FlipE();
        runPos.FlipF();
      } while (runPos != startPos);

      // check if the starting pos or the border has been reached
      assert(runPos == startPos || runPos.IsBorder());
    }

    /**
     * @brief LC_SimulateEdgeCollapse simulates an edge collapse by updating the stars involved.
     * @param edge The edge to collapse.
     */
    void LC_SimulateEdgeCollapse (LCEdge &edge) {
      // let v1 and v2 be the two end points
      LCVertex *v1 = edge.v1;
      LCVertex *v2 = edge.v2;
      assert(v1 && v2);
      LCVertex *v = NULL;
      LCEdge *e = NULL;

      /// v2 merges into v1:
      // star(v1) = star(v1) U star(v2)
      v1->star.insert(v2->star.begin(), v2->star.end());
      v1->star.erase(v2);     // remove v2 from v1-star
      v2->star.erase(v1);     // remove v1 from v2-star
      // foreach v | v2 \in star(v) [i.e. v \in star(v2)]
      //      star(v) = star(v) U {v1} \ {v2}
      for (typename LCVertexStar::iterator vIt = v2->star.begin(); vIt != v2->star.end(); vIt++) {
        v = *vIt;
        v->star.insert(v1);
        v->star.erase(v2);
      }
      /// update the star of the edges which include v1 and v2 in their star
      // foreach e | v1 \in star(e) ^ v2 \in star(e)
      //      star(e) = star(e) \ {v1,v2} U {v1}
      for (typename LCEdgeStar::iterator eIt = v1->edges.begin(); eIt != v1->edges.end(); eIt++) {
        e = *eIt;
        e->star.erase(v2);
      }
      for (typename LCEdgeStar::iterator eIt = v2->edges.begin(); eIt != v2->edges.end(); eIt++) {
        e = *eIt;
        e->star.erase(v2);
        e->star.insert(v1);
      }
    }

    /**
     * @brief The LCVertex struct represents a vertex for the Link Conditions.
     */
    struct LCVertex {
      LCVertexStar star;  // vertex star
      LCEdgeStar edges;   // list of edges whose star involves this vertex
    };

    /**
     * @brief The LCEdge struct represents an edge for the Link Conditions.
     */
    struct LCEdge {
      LCVertex *v1, *v2;          // endpoints
      LCVertexStar star;          // edge star
      LCEdge() {v1 = v2 = NULL;}  // default contructor
    };

    /**
     * @brief _lcVertices is a vector of vertex stars for the link conditions.
     */
    std::vector<LCVertex> _lcVertices;
  };


  // PolyCoordCollapse's methods begin here::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  /**
   * @brief CollapsePolycoord performs all checks and then collapses the polycoord.
   *
   * @warning This function deletes faces and vertices by calling
   * vcg::tri::Allocator<PolyMeshType>::DeleteFace() and
   * vcg::tri::Allocator<PolyMeshType>::DeleteVertex().
   * The object PC_Coords coords is used to track the polycoords, and it has got
   * a size proportional to that of the mesh face container. If you actually
   * delete faces and vertices by calling vcg::tri::Allocator<PolyMeshType>::CompactFaceVector()
   * and vcg::tri::Allocator<PolyMeshType>::CompactVertexVector() after this function,
   * object PC_Coords coords then is not valid any more, so you MUST rearrange it
   * by calling PC_Coords.Reset(). For the same reason, you MUST rearrange LinkConditions linkConditions
   * by calling LinkConditions.Resize().
   * However, for efficiency, you SHOULD compact vertex and face containers at the end of all your
   * polycoord collapsing operations, without having to rearrange coords and linkConditions.
   * The function CollapseAllPolycoords() does this for you.
   *
   * @note Vertex flags, face flags, FF adjacency and FV adjacency are required. Not anything else.
   * Such components are automatically updated here. If the mesh has other components that may be
   * affected by this editing, you should update them later by yourself.
   *
   * @param mesh The polygonal mesh used for getting the face index and deleting the faces
   * (it SHOULD have the vcg::face::PolyInfo component).
   * @param pos Position of the polycoord.
   * @param mark Mark for the current polycoord.
   * @param coords Vector of coords.
   * @param linkConditions Link conditions checker.
   * @param checkSing true if singularities on both sides are not allowed.
   * @return A PC_ResultCode resulting from checks or PC_SUCCESS if the collapse has been performed.
   */
  static PC_ResultCode CollapsePolycoord (PolyMeshType &mesh,
                                          const vcg::face::Pos<FaceType> &pos,
                                          const unsigned long mark,
                                          PC_Coords &coords,
                                          LinkConditions &linkConditions,
                                          const bool checkSing = true) {
    vcg::tri::RequirePerVertexFlags(mesh);
    vcg::tri::RequirePerFaceFlags(mesh);

    if (mesh.face.size() == 0)
      return PC_VOID;

    if (pos.IsNull())
      return PC_VOID;

    vcg::face::Pos<FaceType> tempPos, startPos;

    // check if the sequence of facets is a polycoord and find the starting coord
    PC_ResultCode resultCode = CheckPolycoordFindStartPosition(pos, startPos, checkSing);
    // if not successful, visit the sequence for marking it and return
    if (resultCode != PC_SUCCESS) {
      // if not manifold, visit the entire polycoord ending on the non-manifold edge
      if (resultCode == PC_NOTMANIF) {
        tempPos = pos;
        VisitPolycoord(mesh, tempPos, coords, mark, resultCode);
        if (tempPos.IsManifold() && !tempPos.IsBorder()) {
          tempPos.FlipF();
          VisitPolycoord(mesh, tempPos, coords, mark, resultCode);
        }
        return resultCode;
      }
      // if not quad, visit all the polycoords passing through this coord
      if (resultCode == PC_NOTQUAD) {
        tempPos = startPos;
        do {
          if (!tempPos.IsBorder()) {
            tempPos.FlipF();
            VisitPolycoord(mesh, tempPos, coords, mark, resultCode);
            tempPos.FlipF();
          }
          tempPos.FlipV();
          tempPos.FlipE();
        } while (tempPos != startPos);
      }
      VisitPolycoord(mesh, startPos, coords, mark, resultCode);
      return resultCode;
    }
    // check if the link conditions are satisfied
    bool lc = linkConditions.CheckLinkConditions(mesh, startPos);
    // if not satisfied, visit the sequence for marking it and return
    if (!lc) {
      VisitPolycoord(mesh, startPos, coords, mark, PC_NOLINKCOND);
      return PC_NOLINKCOND;
    }
    // check if the polycoord does not intersect itself
    bool si = IsPolycoordSelfIntersecting(mesh, startPos, coords, mark);
    // if it self-intersects, visit the polycoord for marking it and return
    if (si) {
      VisitPolycoord(mesh, startPos, coords, mark, PC_SELFINTERSECT);
      return PC_SELFINTERSECT;
    }
    // at this point the polycoord is collapsable, visit it for marking
    VisitPolycoord(mesh, startPos, coords, mark, PC_SUCCESS);

    // now collapse
    CoordType point;
    int valenceA = 0, valenceB = 0;
    vcg::face::Pos<FaceType> runPos = startPos;
    vcg::face::JumpingPos<FaceType> tmpPos;
    bool onSideA = false, onSideB = false;
    vcg::face::Pos<FaceType> sideA, sideB;
    typedef std::queue<VertexType **> FacesVertex;
    typedef std::pair<VertexType *, FacesVertex> FacesVertexPair;
    typedef std::queue<FacesVertexPair> FacesVertexPairQueue;
    FacesVertexPairQueue vQueue;
    typedef std::pair<FaceType **, FaceType *> FFpPair;
    typedef std::pair<char *, char> FFiPair;
    typedef std::pair<FFpPair, FFiPair> FFPair;
    typedef std::queue<FFPair> FFQueue;
    FFQueue ffQueue;
    std::queue<VertexType *> verticesToDeleteQueue;
    std::queue<FaceType *> facesToDeleteQueue;

    if (checkSing) {
      do {
        runPos.FlipV();
        valenceB = runPos.NumberOfIncidentVertices();
        tmpPos.Set(runPos.F(), runPos.E(), runPos.V());
        if (tmpPos.FindBorder())
          valenceB++;
        runPos.FlipV();
        valenceA = runPos.NumberOfIncidentVertices();
        tmpPos.Set(runPos.F(), runPos.E(), runPos.V());
        if (tmpPos.FindBorder())
          valenceA++;
        if (valenceA != 4)
          onSideA = true;
        if (valenceB != 4)
          onSideB = true;
        assert(!onSideA || !onSideB);

        if (runPos != startPos && runPos.IsBorder())
          break;

        // go on next edge/face
        runPos.FlipE();
        runPos.FlipV();
        runPos.FlipE();
        runPos.FlipF();
      } while (runPos != startPos);
    }

    runPos = startPos;
    do {
      // compute new vertex
      point = (runPos.V()->P() + runPos.VFlip()->P()) / 2.f;
      if (checkSing) {
        if (onSideA)
          point = runPos.V()->P();
        if (onSideB)
          point = runPos.VFlip()->P();
      }
      runPos.V()->P() = point;
      // list the vertex pointer of the faces on the other side to be updated
      vQueue.push(FacesVertexPair());
      vQueue.back().first = runPos.V();
      tmpPos.Set(runPos.F(), runPos.E(), runPos.V());
      tmpPos.FlipV();
      tmpPos.NextFE();    // go to next face
      while (tmpPos.F() != runPos.F()) {
        if (tmpPos.F() != runPos.FFlip())
          vQueue.back().second.push(&tmpPos.F()->V(tmpPos.VInd()));
        tmpPos.NextFE();    // go to next face
      }

      // enqueue to delete the other vertex
      verticesToDeleteQueue.push(runPos.VFlip());

      // list the adjacencies
      sideA = runPos;
      sideA.FlipE();
      sideA.FlipF();
      sideB = runPos;
      sideB.FlipV();
      sideB.FlipE();
      sideB.FlipF();
      // first side
      if (!sideA.IsBorder()) {
        ffQueue.push(FFPair(FFpPair(),FFiPair()));
        ffQueue.back().first.first = &sideA.F()->FFp(sideA.E());
        ffQueue.back().second.first = &sideA.F()->FFi(sideA.E());
        if (!sideB.IsBorder()) {
          ffQueue.back().first.second = sideB.F();
          ffQueue.back().second.second = sideB.E();
        } else {
          ffQueue.back().first.second = sideA.F();
          ffQueue.back().second.second = sideA.E();
        }
      }
      // second side
      if (!sideB.IsBorder()) {
        ffQueue.push(FFPair(FFpPair(),FFiPair()));
        ffQueue.back().first.first = &sideB.F()->FFp(sideB.E());
        ffQueue.back().second.first = &sideB.F()->FFi(sideB.E());
        if (!sideA.IsBorder()) {
          ffQueue.back().first.second = sideA.F();
          ffQueue.back().second.second = sideA.E();
        } else {
          ffQueue.back().first.second = sideB.F();
          ffQueue.back().second.second = sideB.E();
        }
      }

      // enqueue to delete the face
      facesToDeleteQueue.push(runPos.F());

      // go on next edge/face
      runPos.FlipE();
      runPos.FlipV();
      runPos.FlipE();
      runPos.FlipF();
    } while (runPos != startPos && !runPos.IsBorder());
    assert(runPos == startPos || vcg::face::IsBorder(*startPos.F(),startPos.E()));
    if (runPos.IsBorder()) {
      // compute new vertex on the last (border) edge
      point = (runPos.V()->P() + runPos.VFlip()->P()) / 2.f;
      if (checkSing) {
        if (onSideA)
          point = runPos.V()->P();
        if (onSideB)
          point = runPos.VFlip()->P();
      }
      runPos.V()->P() = point;
      // list the vertex pointer of the faces on the other side to be updated
      vQueue.push(FacesVertexPair());
      vQueue.back().first = runPos.V();
      tmpPos.Set(runPos.F(), runPos.E(), runPos.V());
      tmpPos.FlipV();
      tmpPos.NextFE();    // go to next face
      while (tmpPos.F() != runPos.F()) {
        vQueue.back().second.push(&tmpPos.F()->V(tmpPos.VInd()));
        tmpPos.NextFE();
      }

      // enqueue to delete the other vertex
      verticesToDeleteQueue.push(runPos.VFlip());
    }

    // update vertices
    while (!vQueue.empty()) {
      while (!vQueue.front().second.empty()) {
        *vQueue.front().second.front() = vQueue.front().first;
        vQueue.front().second.pop();
      }
      vQueue.pop();
    }

    // update adjacencies
    while (!ffQueue.empty()) {
      *ffQueue.front().first.first = ffQueue.front().first.second;
      *ffQueue.front().second.first = ffQueue.front().second.second;
      ffQueue.pop();
    }

    // delete faces
    while (!facesToDeleteQueue.empty()) {
      vcg::tri::Allocator<PolyMeshType>::DeleteFace(mesh, *facesToDeleteQueue.front());
      facesToDeleteQueue.pop();
    }

    // delete vertices
    while (!verticesToDeleteQueue.empty()) {
      vcg::tri::Allocator<PolyMeshType>::DeleteVertex(mesh, *verticesToDeleteQueue.front());
      verticesToDeleteQueue.pop();
    }

    return PC_SUCCESS;
  }

  /**
   * @brief CollapseAllPolycoords finds and collapses all the polycoords.
   * @param mesh The input polygonal mesh (it SHOULD have the vcg::face::PolyInfo component).
   * @param checkSing true if singularities on both sides of a polycoord are not allowed.
   */
  static void CollapseAllPolycoords (PolyMeshType &mesh, const bool checkSing = true) {
    vcg::tri::RequireFFAdjacency(mesh);

    if (mesh.FN() == 0)
      return;

    vcg::face::Pos<FaceType> pos;
    PC_ResultCode resultCode;
    std::pair<size_t, unsigned char> face_edge;
    // construct the link conditions checker
    LinkConditions linkConditions(mesh.vert.size());
    // construct the vector of coords
    PC_Coords coords(mesh);
    unsigned long mark = 0;

    // iterate over all the coords
    while (!coords.End()) {
      // get the current coord
      coords.GetCurrent(face_edge);
      // construct a pos on the face and edge of the current coord
      pos.Set(&mesh.face[face_edge.first], face_edge.second, mesh.face[face_edge.first].V(face_edge.second));
      // (try to) collapse the polycoord
      resultCode = CollapsePolycoord(mesh, pos, mark, coords, linkConditions, checkSing);
      // go to the next coord
      coords.Next();

      std::cout << resultCode << std::endl;

      // increment the mark
      mark++;
      if (mark == std::numeric_limits<unsigned long>::max()) {
        coords.ResetMarks();
        mark = 0;
      }
    }
  }

private:
  /**
   * @brief IsVertexAdjacentToAnyNonManifoldEdge checks if a vertex is adjacent to any non-manifold edge.
   * @param pos The starting position.
   * @return true if adjacent to non-manifold edges, false otherwise.
   */
  static bool IsVertexAdjacentToAnyNonManifoldEdge (const vcg::face::Pos<FaceType> &pos) {
    assert(!pos.IsNull());
    vcg::face::JumpingPos<FaceType> jmpPos;
    jmpPos.Set(pos.F(), pos.E(), pos.V());
    do {
      if (!jmpPos.IsManifold())
        return true;
      jmpPos.NextFE();
    } while (jmpPos != pos);
    return false;
  }

  /**
   * @brief CheckPolycoordFindStartPosition checks if it's a collapsable polycoord.
   * @param pos Input The starting position.
   * @param startPos Output the new starting position (in case of borders).
   * @param checkSing true if singularities on both sides are not allowed.
   * @return PC_SUCCESS if it's a collapsable polycoord, otherwise the code for the cause (startPos is on it).
   */
  static PC_ResultCode CheckPolycoordFindStartPosition (const vcg::face::Pos<FaceType> &pos,
                                                        vcg::face::Pos<FaceType> &startPos,
                                                        const bool checkSing = true) {
    assert(!pos.IsNull());
    int valence = 0;
    bool singSideA = false, singSideB = false;
    bool borderA = false, borderB = false;
    bool polyBorderFound = false;
    vcg::face::JumpingPos<FaceType> jmpPos;

    startPos = pos;
    // check if it is a quad
    if (startPos.F()->VN() != 4)
      return PC_NOTQUAD;

    do {
      // check manifoldness
      if (IsVertexAdjacentToAnyNonManifoldEdge(startPos))
        return PC_NOTMANIF;
      startPos.FlipV();
      if (IsVertexAdjacentToAnyNonManifoldEdge(startPos))
        return PC_NOTMANIF;
      startPos.FlipV();

      // check if singularities are not in both sides
      if (checkSing) {
        // compute the valence of the vertex on side B
        startPos.FlipV();
        valence = startPos.NumberOfIncidentVertices();
        // if the vertex is on border increment its valence by 1 (virtually connect it to a dummy vertex)
        jmpPos.Set(startPos.F(), startPos.E(), startPos.V());
        if (jmpPos.FindBorder()) {
          borderB = true;
          valence++;
        }
        if (valence != 4)
          singSideB = true;
        // a 2-valence internl vertex cause a polycoord to touch itself, producing non-2manifoldness
        // in that case, a 2-valence vertex is dealt as 2 singularities in both sides
        if (valence == 2 && !borderB)
          singSideA = true;
        // compute the valence of the vertex on side A
        startPos.FlipV();
        valence = startPos.NumberOfIncidentVertices();
        // if the vertex is on border increment its valence by 1 (virtually connect it to a dummy vertex)
        jmpPos.Set(startPos.F(), startPos.E(), startPos.V());
        if (jmpPos.FindBorder()) {
          borderA = true;
          valence++;
        }
        if (valence != 4)
          singSideA = true;
        // a 2-valence internal vertex cause a polycoord to touch itself, producing non-2manifoldness
        // in that case, a 2-valence vertex is dealt as 2 singularities in both sides
        if (valence == 2 && !borderA)
          singSideB = true;
      }

      // if the first border has been reached, go on the other direction to find the other border
      if (startPos != pos && startPos.IsBorder() && !polyBorderFound) {
        startPos = pos;
        startPos.FlipF();
        polyBorderFound = true;
      }

      // if the other border has been reached, return
      if (polyBorderFound && startPos.IsBorder())
        break;

      // go to the next edge
      startPos.FlipE();
      startPos.FlipV();
      startPos.FlipE();
      // check manifoldness
      if (IsVertexAdjacentToAnyNonManifoldEdge(startPos))
        return PC_NOTMANIF;
      startPos.FlipV();
      if (IsVertexAdjacentToAnyNonManifoldEdge(startPos))
        return PC_NOTMANIF;
      startPos.FlipV();
      // go to the next face
      startPos.FlipF();
    } while (startPos != pos);

    // polycoord with singularities on both sides can not collapse
    if (singSideA && singSideB)
      return PC_SINGBOTH;

    // polycoords that are rings and have borders on both sides can not collapse
    if (!polyBorderFound && borderA && borderB)
      return PC_SINGBOTH;

    return PC_SUCCESS;
  }

  /**
   * @brief IsPolycoordSelfIntersecting checks if the input polycoord intersects itself.
   * @warning Don't call this function without being sure that it's a polycoord
   * (i.e. call CheckPolycoordFindStartPoint() before calling IsPolycoordSelfIntersecting().
   * @param mesh The mesh used for getting the face index.
   * @param startPos The starting position.
   * @param coords The vector of coords.
   * @param mark The current mark, used to identify quads already visited.
   * @return true if it intersects itself, false otherwise.
   */
  static bool IsPolycoordSelfIntersecting (const PolyMeshType &mesh,
                                           const vcg::face::Pos<FaceType> &startPos,
                                           const PC_Coords &coords,
                                           const unsigned long mark) {
    assert(!startPos.IsNull());
    vcg::face::Pos<FaceType> runPos = startPos;
    vcg::face::Pos<FaceType> tmpPos;
    std::pair<size_t, unsigned char> face_edge(std::numeric_limits<size_t>::max(), 0);
    do {
      assert(runPos.F()->VN() == 4);
      // check if we've already crossed this face
      face_edge.first = vcg::tri::Index(mesh, runPos.F());
      face_edge.second = (runPos.E()+1)%2;
      if (coords[face_edge].mark == mark)
        return true;
      // if this coord is adjacent to another coord of the same polycoord
      // i.e., this polycoord touches itself without intersecting
      // it might cause a wrong collapse, producing holes and non-2manifoldness
      tmpPos = runPos;
      tmpPos.FlipE();
      if (!tmpPos.IsBorder()) {
        tmpPos.FlipF();
        face_edge.first = vcg::tri::Index(mesh, tmpPos.F());
        face_edge.second = (tmpPos.E()+1)%2;
        if (coords[face_edge].mark == mark)
          return true;
      }
      tmpPos = runPos;
      tmpPos.FlipV();
      tmpPos.FlipE();
      if (!tmpPos.IsBorder()) {
        tmpPos.FlipF();
        face_edge.first = vcg::tri::Index(mesh, tmpPos.F());
        face_edge.second = (tmpPos.E()+1)%2;
        if (coords[face_edge].mark == mark)
          return true;
      }
      runPos.FlipE();
      runPos.FlipV();
      runPos.FlipE();
      runPos.FlipF();
    } while (runPos != startPos && !runPos.IsBorder());

    return false;
  }

  /**
   * @brief VisitPolycoord updates the information of a polycoord.
   * @param mesh The mesh used for getting the face index.
   * @param startPos The starting position.
   * @param coords The vector of coords.
   * @param mark The mark.
   * @param q The visiting type.
   */
  static void VisitPolycoord (const PolyMeshType &mesh,
                              const vcg::face::Pos<FaceType> &startPos,
                              PC_Coords &coords,
                              const unsigned long mark,
                              const PC_ResultCode q) {
    assert(!startPos.IsNull());
    vcg::face::Pos<FaceType> tmpPos, runPos = startPos;
    std::pair<size_t, unsigned char> face_edge(std::numeric_limits<size_t>::max(), 0);

    if (runPos.F()->VN() != 4)  // non-quads are not visited
      return;

    // follow the sequence of quads
    do {
      // check manifoldness
      tmpPos = runPos;
      do {
        if (!tmpPos.IsManifold()) {
          // update current coord
          face_edge.first = vcg::tri::Index(mesh, tmpPos.F());
          face_edge.second = tmpPos.E()%2;
          coords.UpdateCoord(coords[face_edge], mark, q);
          face_edge.second = (tmpPos.E()+1)%2;
          coords.UpdateCoord(coords[face_edge], mark, q);
          return;
        }
        tmpPos.FlipV();
        tmpPos.FlipE();
      } while (tmpPos != runPos);

      // update current coord
      face_edge.first = vcg::tri::Index(mesh, runPos.F());
      face_edge.second = runPos.E()%2;
      coords.UpdateCoord(coords[face_edge], mark, q);
      // if the polycoord has to collapse, i.e. q == PC_SUCCESS, also visit the orthogonal coord
      if (q == PC_SUCCESS) {
        face_edge.second = (runPos.E()+1)%2;
        coords.UpdateCoord(coords[face_edge], mark, q);
      }

      runPos.FlipE();
      runPos.FlipV();
      runPos.FlipE();
      runPos.FlipF();
    } while (runPos != startPos && !runPos.IsBorder() && runPos.F()->VN() == 4);
    assert(runPos == startPos || vcg::face::IsBorder(*startPos.F(),startPos.E())
           || runPos.F()->VN() != 4 || startPos.FFlip()->VN() != 4);
  }
};

}
}

#endif // POLYGON_POLYCOORD_COLLAPSE_H
