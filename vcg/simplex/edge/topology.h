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

#ifndef _VCG_EDGE_TOPOLOGY
#define _VCG_EDGE_TOPOLOGY

#include <vector>
#include <algorithm>
#include <vcg/simplex/edge/pos.h>

namespace vcg {
namespace edge {
/** \addtogroup edge */
/*@{*/template <class EdgeType>
inline bool IsEdgeManifoldFF( EdgeType const & e, const int j )
{
  assert(e.cFFp(j) != 0); // never try to use this on uncomputed topology

  if(EdgeType::HasFFAdjacency())
    return ( e.cFFp(j) == &e || &e == e.cFFp(j)->cFFp(e.cFFi(j)) );
  else
    return true;
}

/** Return a boolean that indicate if the j-th edge of the face is a border.
  @param j Index of the edge
  @return true if j is an edge of border, false otherwise
*/
template <class EdgeType>
inline bool IsEdgeBorder(EdgeType const & e,  const int j )
{
  if(EdgeType::HasEEAdjacency())
    return e.cEEp(j)==&e;

  assert(0);
  return true;
}

template <class EdgeType>
void VVStarVE(typename EdgeType::VertexType* vp, std::vector<typename EdgeType::VertexType *> &starVec)
{
  starVec.clear();
  edge::VEIterator<EdgeType> vei(vp);
  while(!vei.End())
      {
        starVec.push_back(vei.V1());
        ++vei;
      }
}

template <class EdgeType>
void VEStarVE(const typename EdgeType::VertexType* vp, std::vector<EdgeType *> &starVec)
{
  starVec.clear();
  edge::VEIterator<EdgeType> vei(vp);
  while(!vei.End())
      {
        starVec.push_back(vei.E());
        ++vei;
      }
}

/// Completely detach an edge from the VE adjacency. Useful before deleting it 
template <class EdgeType>
void VEDetach(EdgeType & e)
{
  VEDetach(e,0);
  VEDetach(e,1);
}

/// It detaches the given edge e from the VE adjacency on the vertex z
/// It is used for careful hand stictching of topologies.
template <class EdgeType>
void VEDetach(EdgeType & e, int z)
{
  typename EdgeType::VertexType *vz = e.V(z); // the vertex from which the edge must be detached.
  
  if(vz->VEp()==&e )  //if it is the first edge in the VE chain it detaches it from the begin
    {
        int fz = vz->VEi();
        vz->VEp() = e.VEp(fz);
        vz->VEi() = e.VEi(fz);
    }
    else  // scan the list of edges to find the current edge e to be detached
    {
    VEIterator<EdgeType> x(vz->VEp(),vz->VEi());
    VEIterator<EdgeType> y;

        for(;;)
        {
            y = x;
            ++x;
            assert(x.e!=0);
            if(x.e==&e) // found!
            {
                y.e->VEp(y.z) = e.VEp(z);
                y.e->VEi(y.z) = e.VEi(z);
                break;
            }
        }
    }
}

/// Append an edge in the VE list of vertex e->V(z)
template <class EdgeType>
void VEAppend(EdgeType* & e, int z)
{
    typename EdgeType::VertexType *v = e->V(z);
    if (v->VEp()!=0)
    {
        EdgeType *e0=v->VEp();
        int       z0=v->VEi();
        //append
        e->VEp(z)=e0;
        e->VEi(z)=z0;
    }
    v->VEp()=e;
    v->VEi()=z;
}


/*! Perform a simple edge collapse using VE adjacency
 * 
 * It collapses the two edges incidnent on the indicated vertex so that the passed edge survives, 
 * the indicated vertex is deleted, and the edge ajacent to e0 along z is deleted too.
 * It assumes that the edge mesh is 1-Manifold.
 * If the indicated vertex <vd> is boundary or non manifold the function do nothing.
 * 
 *    v0     vd       v1
 * ---O-------O-------O---
 *   z0   e0  z   e1  z1 
 * 
 *    v0              v1
 * ---O---------------O---
 *           e0
 * 
 * 
 */
template <class MeshType>
void VEEdgeCollapse(MeshType &poly, typename MeshType::EdgeType *e0, const int z)
{
  typedef typename MeshType::EdgeType EdgeType;
  typedef typename MeshType::VertexType VertexType;
  
  VertexType *vd = e0->V(z);
  
  std::vector<EdgeType *> starVecEp;
  edge::VEStarVE(vd,starVecEp);   
  if(starVecEp.size()!=2) return;
  
  EdgeType *e1=0; // this edge will be deleted
  if( starVecEp[0] == e0 ) e1 = starVecEp[1];
  if( starVecEp[1] == e0 ) e1 = starVecEp[0];
  assert(e1 && (e1!=e0) );

  //int z0 = (z+1)%2;
  int z1 = -1;
  if(e1->V(0) == vd) z1=1;
  if(e1->V(1) == vd) z1=0;
  assert(z1!=-1);
  
  VertexType *v1 = e1->V(z1); 
      
  edge::VEDetach(*e1);
  edge::VEDetach(*e0,z);
  e0->V(z) = v1;
  edge::VEAppend(e0, z);
  tri::Allocator<MeshType>::DeleteEdge(poly,*e1);
  tri::Allocator<MeshType>::DeleteVertex(poly,*vd);
}

template <class MeshType>
void VEEdgeCollapse(MeshType &poly, typename MeshType::VertexType *v)
{
  VEEdgeCollapse(poly,v->VEp(),v->VEi());
}



/*! Returns the number of incident edges over a vertex vp; Using the VE adjacency.  
 * 
 * It just follows the chain of incident edges of the VE adjacency.
*/

template <class EdgeType>
int VEDegree(const typename EdgeType::VertexType* vp)
{
  int cnt=0;
  edge::VEIterator<EdgeType> vei(vp);
  while(!vei.End()) 
  {
    ++cnt;
    ++vei;
  }    
  return cnt;
}


} // end namespace edge
} // end namespace vcg


#endif
