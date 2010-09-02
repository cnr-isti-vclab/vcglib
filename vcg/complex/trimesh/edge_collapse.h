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
/****************************************************************************
  History
    $Log: not supported by cvs2svn $
    Revision 1.16  2006/10/07 15:04:25  cignoni
    removed a useless include

    Revision 1.15  2005/10/12 10:36:26  cignoni
    Removed unused local type Edge. Now it use the standard simplex edge.

    Revision 1.14  2004/12/10 01:04:42  cignoni
    better comments

    Revision 1.13  2004/11/23 10:34:45  cignoni
    passed parameters by reference in many funcs and gcc cleaning


****************************************************************************/

#ifndef __VCG_TETRA_TRI_COLLAPSE
#define __VCG_TETRA_TRI_COLLAPSE


#include<vcg/simplex/face/pos.h>
#include<vcg/simplex/face/topology.h>
#include<vcg/complex/trimesh/allocate.h>

namespace vcg{
namespace tri{	

/** \addtogroup trimesh */
/*@{*/
/** This a static utility class for the edge collapse. 
    It provides a common set of useful function for actually making an edge collapse over a trimesh. 
    See also the corresponding class in the local optimization framework called TriEdgeCollapse
**/

template <class TRI_MESH_TYPE> 
class EdgeCollapse
{
	public:
  /// The tetrahedral mesh type
  typedef	 TRI_MESH_TYPE TriMeshType;
  /// The tetrahedron type
  typedef	typename TriMeshType::FaceType FaceType;
	/// The vertex type
	typedef	typename FaceType::VertexType VertexType;
  typedef	typename FaceType::VertexPointer VertexPointer;
  /// The vertex iterator type
  typedef	typename TriMeshType::VertexIterator VertexIterator;
  /// The tetra iterator type
  typedef	typename TriMeshType::FaceIterator FaceIterator;
  /// The coordinate type
	typedef	typename FaceType::VertexType::CoordType CoordType;
  /// The scalar type
  typedef	typename TriMeshType::VertexType::ScalarType ScalarType;
  ///the container of tetrahedron type
  typedef typename TriMeshType::FaceContainer FaceContainer;
  ///the container of vertex type
  typedef typename TriMeshType::VertContainer VertContainer;
  ///half edge type
	typedef typename TriMeshType::FaceType::EdgeType EdgeType;
	/// vector of pos
	typedef typename std::vector<EdgeType> EdgeVec;
	///of VFIterator
	typedef typename vcg::face::VFIterator<FaceType>  VFI;
	/// vector of VFIterator
	typedef typename std::vector<vcg::face::VFIterator<FaceType> > VFIVec;



	/// Default Constructor
	EdgeCollapse()
		{
		};

  ~EdgeCollapse()
		{
		};

	static VFIVec & AV0(){static VFIVec av0; return av0;}
	static VFIVec & AV1(){static VFIVec av1; return av1;}
	static VFIVec & AV01(){static VFIVec av01; return av01;}


	void FindSets(EdgeType &p)
	{
		VertexType * v0 = p.V(0);
		VertexType * v1 = p.V(1);

		AV0().clear();  // Facce incidenti in v0 
		AV1().clear();  // Facce incidenti in v1 
		AV01().clear(); // Facce incidenti in v0 e v1
		
		VFI x;

		for( x.f = v0->VFp(), x.z = v0->VFi(); x.f!=0; ++x)
		{
			int zv1 = -1;

			for(int j=0;j<3;++j) 
				if( x.f->V(j)==&*v1 )	{
					zv1 = j;
					break;
				}
			if(zv1==-1) 	AV0().push_back( x ); // la faccia x.f non ha il vertice v1 => e' incidente solo in v0
			else    			AV01().push_back( x );
		}

		for( x.f = v1->VFp(), x.z = v1->VFi(); x.f!=0; ++x )
		{
			int zv0 = -1;

			for(int j=0;j<3;++j)
				if( x.f->V(j)==&*v0 )	{
					zv0 = j;
					break;
				}
			if(zv0==-1)	AV1().push_back( x ); // la faccia x.f non ha il vertice v1 => e' incidente solo in v0
		}
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
  bool LinkConditions(EdgeType pos)
  {
    typedef typename vcg::face::VFIterator<FaceType> VFIterator;
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



  bool LinkConditionsOld(EdgeType  pos){
		
		const int ADJ_1 = TriMeshType::VertexType::NewBitFlag();
		const int ADJ_E = TriMeshType::VertexType::NewBitFlag();
		//enum {ADJ_1= MeshType::VertexType::USER0,
		//	    ADJ_E= MeshType::VertexType::USER0<<1} ;
		// const int ALLADJ = ADJ_1|ADJ_E;
		const int NOTALLADJ = ~(ADJ_1 | ADJ_E | TriMeshType::VertexType::VISITED);
		const int NOTALLADJ1 = ~(ADJ_E | TriMeshType::VertexType::VISITED);

		//EdgePosB<MeshType::face_type::face_base> x;
		typename vcg::face::VFIterator<FaceType> x;
		// Clear visited and adj flag for all vertices adj to v0;
		for(x.f = pos.V(0)->VFp(), x.z = pos.V(0)->VFi(); x.f!=0; ++x )	 	{
			x.f->V1(x.z)->Flags() &= NOTALLADJ;
			x.f->V2(x.z)->Flags() &= NOTALLADJ;
		}
		// Clear visited flag for all vertices adj to v1 and set them adj1 to v1;
		for(x.f = pos.V(1)->VFp(), x.z = pos.V(1)->VFi(); x.f!=0; ++x )	 	{
			x.f->V1(x.z)->Flags() &= NOTALLADJ1;
			x.f->V2(x.z)->Flags() &= NOTALLADJ1;
		}
		// Mark vertices adj to v1 as  ADJ_1 and  adj1 to v1;
		for(x.f = pos.V(1)->VFp(), x.z = pos.V(1)->VFi(); x.f!=0; ++x )	 	{
			if(x.f->V1(x.z)==pos.V(0)) x.f->V2(x.z)->Flags() |= ADJ_E | ADJ_1;
			else x.f->V2(x.z)->Flags() |= ADJ_1;
			if(x.f->V2(x.z)==pos.V(0)) x.f->V1(x.z)->Flags() |= ADJ_E | ADJ_1;
			else x.f->V1(x.z)->Flags() |= ADJ_1;
		}
			
		// compute the number of:
		int adj01=0;  // vertices adjacents to both v0 and v1 
		int adje=0;   // vertices adjacents to an egde (usually 2)
		for(x.f = pos.V(0)->VFp(), x.z = pos.V(0)->VFi(); x.f!=0; ++x )	 	{
			if(!x.f->V1(x.z)->IsV()) {
				x.f->V1(x.z)->SetV();
				if(x.f->V1(x.z)->Flags()&ADJ_1) ++adj01;
				if(x.f->V1(x.z)->Flags()&ADJ_E) ++adje;
			}
			if(!x.f->V2(x.z)->IsV()) {
				x.f->V2(x.z)->SetV();
				if(x.f->V2(x.z)->Flags()&ADJ_1) ++adj01;
				if(x.f->V2(x.z)->Flags()&ADJ_E) ++adje;
			}
		}

		//bool val=TopoCheck2();
		//if(val != (adj01==adje)) 		printf("Wrong topo %i %i\n",adj01,adje);
		TriMeshType::VertexType::DeleteBitFlag(ADJ_E);
		TriMeshType::VertexType::DeleteBitFlag(ADJ_1);
		
    return (adj01==adje);
	}



	int DoCollapse(TriMeshType &m, EdgeType & c, const Point3<ScalarType> &p)
	{
		FindSets(c);
		typename VFIVec::iterator i;
		int n_face_del =0 ;
		
		//set Face Face topology
		if (TriMeshType::HasFFTopology())
		{
			//int e0=c.z;
			//int e1=c.f->FFi(c.z);	//opposite edge

			//FaceType *f0=c.f;
			//FaceType *f1=f0->FFp(c.z);
			//
			////take right indexes
			//FaceType *f00=f0->FFp((e0+1)%3);
			//FaceType *f01=f0->FFp((e0+2)%3);
			//int If00=f0->FFi((e0+1)%3);
			//int If01=f0->FFi((e0+2)%3);
			//
			////then attach faces
			//f00->FFp(If00)=f01;
			//f00->FFi(If00)=If01;
			//f01->FFp(If01)=f00;
			//f01->FFi(If01)=If00;

			////and the ones of face f1

			//f00=f1->FFp((e1+1)%3);
			//f01=f1->FFp((e1+2)%3);
			//If00=f1->FFi((e1+1)%3);
			//If01=f1->FFi((e1+2)%3);
			//
			////and attach faces
			//f00->FFp(If00)=f01;
			//f00->FFi(If00)=If01;
			//f01->FFp(If01)=f00;
			//f01->FFi(If01)=If00;
		}

		for(i=AV01().begin();i!=AV01().end();++i)
		{	
			FaceType  & f = *((*i).f);
			assert(f.V((*i).z) == c.V(0));
			vcg::face::VFDetach(f,((*i).z+1)%3);
			vcg::face::VFDetach(f,((*i).z+2)%3);
			Allocator<TriMeshType>::DeleteFace(m,f);
			//n_face_del++;
		}

		//set Vertex Face topology
		for(i=AV0().begin();i!=AV0().end();++i)
		{
			(*i).f->V((*i).z) = c.V(1);									 // In tutte le facce incidenti in v0, si sostituisce v0 con v1
			(*i).f->VFp((*i).z) = (*i).f->V((*i).z)->VFp(); // e appendo la lista di facce incidenti in v1 a questa faccia
			(*i).f->VFi((*i).z) = (*i).f->V((*i).z)->VFi();
			(*i).f->V((*i).z)->VFp() = (*i).f;
			(*i).f->V((*i).z)->VFi() = (*i).z;
		}
		
		Allocator<TriMeshType>::DeleteVertex(m,*(c.V(0)));
		//c.V(0)->SetD();
		c.V(1)->P()=p;
		return n_face_del;
	}

};

}
}
#endif 
