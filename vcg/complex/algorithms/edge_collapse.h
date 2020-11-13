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
#ifndef _VCG_EDGE_COLLAPSE_
#define _VCG_EDGE_COLLAPSE_


#include<vcg/simplex/face/pos.h>
#include<vcg/simplex/face/topology.h>

namespace vcg{
namespace tri{

template < class VERTEX_TYPE>
class BasicVertexPair {
public:
  inline BasicVertexPair() {}
  inline BasicVertexPair( VERTEX_TYPE * v0, VERTEX_TYPE * v1){V(0) = v0; V(1) = v1; }
  void Sort() {if(V(0)<V(0)) std::swap(V(0),V(0)); }
  VERTEX_TYPE *&V(int i) { return v[i]; }
  VERTEX_TYPE *cV(int i) const { return v[i]; }
private:
  VERTEX_TYPE *v[2]; // remember that v[0] will be deleted and v[1] will survive (eventually with a new position)
};


/** \addtogroup trimesh */
/*@{*/
/** This a static utility class for the edge collapse.
    It provides a common set of useful function for actually making an edge collapse over a trimesh.
    See also the corresponding class in the local optimization framework called TriEdgeCollapse
**/

template <class TRI_MESH_TYPE, class VertexPair>
class EdgeCollapser
{
public:
  typedef	TRI_MESH_TYPE TriMeshType;
  typedef	typename TriMeshType::FaceType FaceType;
  typedef	typename FaceType::VertexType VertexType;
  typedef	typename FaceType::VertexPointer VertexPointer;
  typedef	typename FaceType::VertexType::CoordType CoordType;
  typedef	typename TriMeshType::VertexType::ScalarType ScalarType;
  typedef typename vcg::face::VFIterator<FaceType>  VFIterator;
  typedef typename std::vector<vcg::face::VFIterator<FaceType> > VFIVec;

private:
  struct EdgeSet
  {
     VFIVec av0, av1, av01;
     VFIVec & AV0() { return av0;}  // Faces incident only on v0
//     VFIVec & AV1() { return av1;}  // Faces incident only on v1
     VFIVec & AV01(){ return av01;} // Faces incident only on both v0 and v1
  };

  static void FindSets(VertexPair &p, EdgeSet &es)
  {
    VertexType * v0 = p.V(0);
    VertexType * v1 = p.V(1);
    
    es.AV0().clear(); 
//    es.AV1().clear(); 
    es.AV01().clear();
    
    for(VFIterator x = VFIterator(v0); !x.End(); ++x)
    {
      bool foundV1=false;      
      for(int j=0;j<3;++j)
        if( x.f->V(j)==v1 )	{
          foundV1 = true;
          break;
        }
      if(!foundV1)  es.AV0().push_back( x ); // v1 not found -> so the face is incident only on v0
      else         es.AV01().push_back( x );
    }
    
//    for( VFIterator x = VFIterator(v1); !x.End(); ++x)
//    {
//      bool foundV0=false;      
//      for(int j=0;j<3;++j)
//        if( x.f->V(j)==v0 )	{
//          foundV0=true;
//          break;
//        }
//      if(!foundV0)	es.AV1().push_back( x ); // v0 not found -> so the face is incident only on v0
//    }
  }
  
  /*
    Link Conditions test, as described in

    Topology Preserving Edge Contraction
    T. Dey, H. Edelsbrunner,
    Pub. Inst. Math. 1999

    Lk (sigma) is the set of all the faces of the cofaces of sigma that are disjoint from sigma

    Lk(v0) inters Lk(v1) == Lk(v0-v1)

    To perform these tests using only the VF adjacency we resort to some virtual counters over
    the vertices and the edges, we implement them as std::maps, and we increase these counters
    by running over all the faces around each vertex of the collapsing edge.

    At the end (after adding dummy stuff) we should have
       2 for vertices not shared
       4 for vertices shared
       2 for edges shared
       1 for edges not shared.


*/

public:
  static bool LinkConditions(VertexPair &pos)
  {
    // at the end of the loop each vertex must be counted twice
    // except for boundary vertex.
    std::map<VertexPointer,int> VertCnt;
    std::map<std::pair<VertexPointer,VertexPointer>,int> EdgeCnt;

    // the list of the boundary vertexes for the two endpoints
    std::vector<VertexPointer> BoundaryVertexVec[2];

    // Collect vertexes and edges of V0 and V1
    VFIterator vfi;
    for(int i=0;i<2;++i)
    {
      vfi = VFIterator(pos.V(i));
      for( ;!vfi.End();++vfi)
      {
        ++ VertCnt[vfi.V1()];
        ++ VertCnt[vfi.V2()];
        if(vfi.V1()<vfi.V2()) ++EdgeCnt[std::make_pair(vfi.V1(),vfi.V2())];
                         else ++EdgeCnt[std::make_pair(vfi.V2(),vfi.V1())];
      }
      // Now a loop to add dummy stuff: add the dummy vertex and two dummy edges
      // (and remember to increase the counters for the two boundary vertexes involved)
      typename std::map<VertexPointer,int>::iterator vcmit;
      for(vcmit=VertCnt.begin();vcmit!=VertCnt.end();++vcmit)
      {
        if((*vcmit).second==1) // boundary vertexes are counted only once
          BoundaryVertexVec[i].push_back((*vcmit).first);
      }
      if(BoundaryVertexVec[i].size()==2)
      { // aha! one of the two vertex of the collapse is on the boundary
        // so add dummy vertex and two dummy edges
        VertCnt[0]+=2;
        ++ EdgeCnt[std::make_pair(VertexPointer(0),BoundaryVertexVec[i][0]) ] ;
        ++ EdgeCnt[std::make_pair(VertexPointer(0),BoundaryVertexVec[i][1]) ] ;
        // remember to hide the boundaryness of the two boundary vertexes
        ++VertCnt[BoundaryVertexVec[i][0]];
        ++VertCnt[BoundaryVertexVec[i][1]];
      }
    }

    // Final loop to find cardinality of Lk( V0-V1 )
    // Note that Lk(edge) is only a set of vertices.
    std::vector<VertexPointer> LkEdge;

    for( vfi = VFIterator(pos.V(0)); !vfi.End(); ++vfi)
    {
      if(vfi.V1() == pos.V(1) ) LkEdge.push_back(vfi.V2());
      if(vfi.V2() == pos.V(1) ) LkEdge.push_back(vfi.V1());
    }

    // if the collapsing edge was a boundary edge, we must add the dummy vertex.
    // Note that this implies that Lk(edge) >=2;
    if(LkEdge.size()==1)
    {
      LkEdge.push_back(0);
    }

    // NOW COUNT!!!
    size_t SharedEdgeCnt=0;
    typename std::map<std::pair<VertexPointer,VertexPointer>, int>::iterator eci;
    for(eci=EdgeCnt.begin();eci!=EdgeCnt.end();++eci)
      if((*eci).second == 2) SharedEdgeCnt ++;

    if(SharedEdgeCnt>0) return false;
    size_t SharedVertCnt=0;
    typename std::map<VertexPointer,int>::iterator vci;
    for(vci=VertCnt.begin();vci!=VertCnt.end();++vci)
      if((*vci).second == 4) SharedVertCnt++;

    if(SharedVertCnt != LkEdge.size() ) return false;

    return true;
  }

  // Main Collapsing Function: the one that actually performs the collapse of the edge denoted by the VertexPair c
  // Remember that v[0] will be deleted and v[1] will survive with the position indicated by p
  // To do a collapse onto a vertex simply pass p as the position of the surviving vertex
  static int Do(TriMeshType &m, VertexPair & c, const Point3<ScalarType> &p, const bool preserveFaceEdgeS = false)
  {
      EdgeSet es, es1;
      FindSets(c,es);

      if (preserveFaceEdgeS)
      {
          VertexPair c1(c.V(1), c.V(0));
          FindSets(c1, es1);
      }

      int n_face_del=0 ;

      static int VtoE[3][3] = { -1,  0,  2,
                                0, -1,  1,
                                2,  1, -1 };

      std::vector<VertexPointer> topVertices; topVertices.reserve(2);
      std::vector<VertexPointer> fan1V2S; fan1V2S.reserve(2);
      std::vector<VertexPointer> v2s; v2s.reserve(2);
      std::map <VertexPointer, bool> toSel;


      for(auto i=es.AV01().begin();i!=es.AV01().end();++i)
      {
          FaceType  & f = *((*i).f);
          assert(f.V((*i).z) == c.V(0));

          if (preserveFaceEdgeS)
          {
              VertexPointer top;
              size_t topIdx;
              if (f.V(((*i).z+1)%3) == c.V(1))
              {
                  top = f.V(((*i).z+2)%3);
                  topIdx = ((*i).z+2)%3;
              }
              else
              {
                  top = f.V(((*i).z+1)%3);
                  topIdx = ((*i).z+1)%3;
              }

              topVertices.push_back(top);

              if (f.IsFaceEdgeS(VtoE[((*i).z)][topIdx]))
                  fan1V2S.push_back(top);

              if (f.IsFaceEdgeS(VtoE[((*i).z+1)%3][((*i).z+2)%3]))
                  v2s.push_back(top);
          }

          vcg::face::VFDetach(f,((*i).z+1)%3);
          vcg::face::VFDetach(f,((*i).z+2)%3);
          Allocator<TriMeshType>::DeleteFace(m,f);
          n_face_del++;
      }

      // Very LOW LEVEL update of VF Adjacency;
      // for all the faces incident in v[0]
      // - v[0] will be deleted so we substitute v[0] with v[1]
      // - we prepend that face to the list of the faces incident on v[1]
      for(auto i=es.AV0().begin();i!=es.AV0().end();++i)
      {
          FaceType  & f = *((*i).f);

          if (preserveFaceEdgeS)
          {
              for (size_t j = 0; j < v2s.size(); ++j)
              {
                  if ((*i).f->V(((*i).z+1)%3) == v2s[j])
                  {
                      (*i).f->SetFaceEdgeS(VtoE[((*i).z)%3][((*i).z+1)%3]);
                      break;
                  }
                  if ((*i).f->V(((*i).z+2)%3) == v2s[j])
                  {
                      (*i).f->SetFaceEdgeS(VtoE[((*i).z)%3][((*i).z+2)%3]);
                      break;
                  }
              }
          }
          (*i).f->V((*i).z) = c.V(1);	// For each face in v0 we substitute v0 with v1
          (*i).f->VFp((*i).z) = c.V(1)->VFp();
          (*i).f->VFi((*i).z) = c.V(1)->VFi();
          c.V(1)->VFp() = (*i).f;
          c.V(1)->VFi() = (*i).z;
      }

      if (preserveFaceEdgeS)
      {
          for (auto i = es1.AV0().begin(); i != es1.AV0().end(); ++i)
          {
              FaceType  & f = *((*i).f);
              for (size_t j = 0; j < fan1V2S.size(); ++j)
              {
                  if ((*i).f->V(((*i).z+1)%3) == fan1V2S[j])
                  {
                      (*i).f->SetFaceEdgeS(VtoE[((*i).z)%3][((*i).z+1)%3]);
                      break;
                  }
                  if ((*i).f->V(((*i).z+2)%3) == fan1V2S[j])
                  {
                      (*i).f->SetFaceEdgeS(VtoE[((*i).z)%3][((*i).z+2)%3]);
                      break;
                  }
              }
          }
      }

      Allocator<TriMeshType>::DeleteVertex(m,*(c.V(0)));
      c.V(1)->P()=p;
      return n_face_del;
  }
  
};

} // end namespace tri
} // end namespace vcg
#endif
