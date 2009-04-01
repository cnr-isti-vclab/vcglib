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

#ifndef __VCGLIB_EDGE_SUPPORT
#define __VCGLIB_EDGE_SUPPORT

#include <vector>
#include <vcg/complex/trimesh/allocate.h>
#include <vcg/complex/trimesh/clean.h>
#include <vcg/complex/trimesh/update/topology.h>
#include <vcg/complex/trimesh/base.h>

namespace vcg
{
	namespace tri{
		/// \ingroup trimesh 

		/// \headerfile edge_support.h vcg/complex/trimesh/edge_support.h

		/// \brief This class is used to build edge based data structure from indexed data structure and viceversa

		/**
		*/

		template <class MeshType  >
		class EdgeSupport{
		public:
			typedef typename MeshType::VertexType VertexType;
			typedef typename MeshType::VertexPointer VertexPointer;
			typedef typename MeshType::EdgePointer EdgePointer;
			typedef typename MeshType::EdgeType EdgeType;
			typedef typename MeshType::EdgeIterator EdgeIterator;
			typedef typename MeshType::FaceIterator FaceIterator;
			typedef typename MeshType::FaceType FaceType;

			struct VertexPairEdgePtr{
				VertexPairEdgePtr(VertexPointer _v0,VertexPointer _v1,EdgePointer _ep):v0(_v0),v1(_v1),ep(_ep){if(v0>v1) std::swap(v0,v1);}
                                const bool operator <(const VertexPairEdgePtr &o) const {return (v0 == o.v0)? (v1<o.v1):(v0<o.v0);}
                                const bool operator ==(const VertexPairEdgePtr &o) const {return (v0 == o.v0)&& (v1==o.v1);}

				VertexPointer v0,v1;
				EdgePointer ep;
			};
			struct FacePtrInt{
				FacePtrInt ( FaceType * _f,int _i):f(_f),i(_i){}
				FaceType * f;
				int i;
			};

			typedef std::vector<bool> BitVector;

			/**
			build a half-edge data structure from an indexed data structure. Note that the half-edges are allocated here for the first time.
			If you have a mesh where there are already edges, they will be removed and the data lost, so do not use this function
			to just "update" the topology of half edges.
			**/
			static void ComputeHalfEdgeFromIndexed(MeshType & m){
				assert(HasFVAdjacency(m));
				assert(MeshType::EdgeType::HasHENextAdjacency());
				assert(MeshType::EdgeType::HasHEOppAdjacency());			

                                typename MeshType::template PerFaceAttributeHandle<BitVector> flagVisited =
                                        vcg::tri::Allocator<MeshType>::template AddPerFaceAttribute<BitVector>(m,"");
				std::vector<FacePtrInt > borderEdges;

				// allocate all new half edges
				FaceIterator fi;
				int n_edges = 0;

				// count how many half edge to allocate
				for(fi = m.face.begin(); fi != m.face.end(); ++fi) if(! (*fi).IsD()) 
				{n_edges+=(*fi).VN();	
				 for(int i = 0; i < (*fi).VN(); ++i)
					 if(vcg::face::IsBorder<FaceType>((*fi),(i)))
						++n_edges;
				}

				// allocate the half edges
				typename MeshType::EdgeIterator ei = vcg::tri::Allocator<MeshType>::AddEdges(m,n_edges);
				
				std::vector<VertexPairEdgePtr> all;
				int firstEdge = 0;
				for(fi = m.face.begin(); fi != m.face.end(); ++fi)if(!(*fi).IsD()){
					assert((*fi).VN()>2);
					if(flagVisited[*fi].empty()) {flagVisited[*fi].resize((*fi).VN());}

					for(int i  = 0; i < (*fi).VN(); ++i,++ei)
					{
						(*ei).HEVp()	= (*fi).V(i);
						(*ei).HENp()	= &m.edge[firstEdge + (i +1) % (*fi).VN()];
						if(MeshType::EdgeType::HasEFAdjacency()) 
							(*ei).EFp() = &(*fi);
						if( MeshType::FaceType::HasFHEAdjacency())
							(*fi).FHEp() = &(*ei);
						if(MeshType::EdgeType::HasHEPrevAdjacency())
							(*ei).HEPp()	= &m.edge[firstEdge + (i +(*fi).VN()-1) % (*fi).VN()];
						if(HasVEAdjacency(m))
							(*ei).HEVp()->VEp() = &(*ei);
						all.push_back(VertexPairEdgePtr((*fi).V(i), (*fi).V((*fi).Next(i)),&(*ei)));// it will be used to link the hedges

						if(  vcg::face::IsBorder<FaceType>((*fi),(i)))
							borderEdges.push_back(FacePtrInt(&(*fi),i));
					}
					firstEdge += (*fi).VN();
				}

				// add all the border edges
				int borderLength; 
                                typename std::vector<FacePtrInt >::iterator ebi;
				for( ebi = borderEdges.begin(); ebi != borderEdges.end(); ++ebi)
					if( !flagVisited[(*ebi).f][(*ebi).i])// not already inserted
					{
						 
						borderLength = 0;
						vcg::face::Pos<FaceType> bp((*ebi).f,(*ebi).i);
						FaceType * start = (*ebi).f;
						do{
							all.push_back( VertexPairEdgePtr ( bp.f->V( bp.f->Next(bp.z) ),bp.f->V( bp.z ),&(*ei)));
							(*ei).HEVp()	=  bp.f->V(bp.f->Next(bp.z)) ;
							flagVisited[bp.f][bp.z] = true;
							++ei;
							bp.NextB();
							++borderLength;
							}while (bp.f != start);
				 
						// run over the border edges to link the adjacencies
						for(int be = 0; be < borderLength; ++be){
							if(MeshType::EdgeType::HasEFAdjacency()) 
								m.edge[firstEdge + be].EFp() = NULL;
							if(MeshType::EdgeType::HasHEPrevAdjacency())
								m.edge[firstEdge + be].HEPp()	= &m.edge[firstEdge + (be +borderLength-1) % borderLength];
							m.edge[firstEdge + be].HENp()	= &m.edge[firstEdge + (be +1) % borderLength];
						}
						firstEdge+=borderLength;
				}
                                vcg::tri::Allocator<MeshType>:: template DeletePerFaceAttribute<BitVector>(m,flagVisited );
				
				std::sort(all.begin(),all.end());
				assert(all.size() == n_edges);
				for(int i = 0 ; i < all.size(); )
					if(all[i]  == all[i+1]) 
					{
						all[i].ep->HEOp()	= all[i+1].ep;
						all[i+1].ep->HEOp() = all[i].ep;
						i+=2;
					} 
					else
					{
						all[i].ep->HEOp() = all[i].ep;
						i+=1;
					}

			}
			/**
			builds an indexed data structure from a half-edge data structure.
			Note: if  the half edge have the pointer to face   
			their relation FV (face-vertex) will be computed and the data possibly stored in the
			face will be preserved.
			**/
			static void ComputeIndexedFromHalfEdge(  MeshType & m ){
				assert(HasFVAdjacency(m));
				assert(MeshType::EdgeType::HasHENextAdjacency());
				assert(MeshType::EdgeType::HasHEVAdjacency());
				assert(MeshType::EdgeType::HasHEOppAdjacency());			
				assert(MeshType::FaceType::HasFHEAdjacency());
				bool createFace,hasHEF,hasFHE;

				typename MeshType::template PerEdgeAttributeHandle<bool> hV = Allocator<MeshType>::template AddPerEdgeAttribute<bool>(m,"");

				typename MeshType::EdgeIterator ei;
				typename MeshType::FacePointer fp;
				typename MeshType::FaceIterator fi;
				typename MeshType::EdgePointer ep,epF;
				int vi = 0;

				hasHEF = (MeshType::EdgeType::HasEFAdjacency());
				assert( !hasHEF || (hasHEF && m.fn>0));

				// if the edgetype has the pointer to face   
				// it is assumed the the edget2face pointer (HEFp) are correct
				// and the faces are allocated
					for ( ei = m.edge.begin(); ei != m.edge.end(); ++ei)
					if(!(*ei).IsD())								// it has not been deleted
					if(!hasHEF || ( hasHEF &&  (*ei).EFp()!=NULL))	// if it has a pointer to the face it is 
																	// not null (i.e. it is not a border edge)
					if(!hV[(*ei)] )									// it has not be visited yet
					{
						if(!hasHEF)// if it has 
							fp = &(* Allocator<MeshType>::AddFaces(m,1));
						else
							fp = (*ei).EFp();

						ep = epF = &(*ei);
						std::vector<VertexPointer> vpts;
						do{vpts.push_back((*ep).HEVp()); ep=ep->HENp();}while(ep!=epF);
						int idbg  =fp->VN();
						if(fp->VN() != vpts.size()){
							fp->Dealloc();
							fp ->Alloc(vpts.size());
						}
						int idbg1  =fp->VN();
						for(int i  = 0; i < vpts.size();++i) fp ->V(i) = vpts[i];// set the pointer from face to vertex

  					hV[(*ei)] = true;
				}
				Allocator<MeshType>::DeletePerEdgeAttribute(m,hV);
			}

	/**
	Checks pointers FHEp() are valid
	**/
	static bool CheckConsistency_FHEp(MeshType &  m){
		assert(MeshType::FaceType::HasFHEAdjacency());
		FaceIterator fi;
		for(fi = m.face.begin(); fi != m.face.end(); ++fi)
			if(!(*fi).IsD()){
				if((*fi).FHEp() <  &(*m.edge.begin())) return false;
				if((*fi).FHEp() >  &(m.edge.back())) return false;
			}
		return true;
	}

	/**
	Checks that half edges and face relation are consistent
	**/
	static bool CheckConsistency(MeshType & m){
		assert(MeshType::EdgeType::HasHENextAdjacency());
		assert(MeshType::EdgeType::HasHEOppAdjacency());			
		assert(MeshType::EdgeType::HasHEVAdjacency());
		assert(MeshType::FaceType::HasFHEAdjacency());

		bool hasHEF = ( MeshType::EdgeType::HasEFAdjacency());
		bool hasHEP = ( MeshType::EdgeType::HasHEPrevAdjacency());

		FaceIterator fi;
		EdgePointer ep,ep1;
		int cnt = 0;
		if(( MeshType::EdgeType::HasEFAdjacency())){
			int iDb = 0;
			for(fi = m.face.begin(); fi != m.face.end(); ++fi,++iDb)
				if(!(*fi).IsD())
				{
					ep = ep1 = (*fi).FHEp();
					do{
						if(ep->IsD()) 
							return false; // the edge should not be connected, it has been deleted
	 					if(ep->EFp() != &(*fi)) 
							return false;// edge is not pointing to the rigth face
						ep = ep->HENp();
						if(cnt++ > m.en)
							return false; // edges are ill connected (HENp())
					}while(ep!=ep1);
				}
		}

		EdgePointer epPrev;
		EdgeIterator ei;
		bool extEdge ;
		for( ei  = m.edge.begin(); ei != m.edge.end(); ++ei)
			if(!(*ei).IsD())
		{
			cnt = 0;
			epPrev = ep = ep1 = &(*ei);
			do{
				extEdge = (ep->EFp()==NULL);
				if(hasHEP){
					if( ep->HENp()->HEPp() != ep) 
						return false; // next and prev relation are not mutual
					if( ep->HEPp() == ep) 
						return false; // the previous of an edge cannot be the edge itself
				}
				if( ep->HEOp()  == ep) 
					return false; // opposite relation is not mutual
				if( ep->HEOp()->HEOp() != ep) 
					return false; // opposite relation is not mutual
				if(ep->HENp() == ep)
					return false; //  the next of an edge cannot be the edge itself
				ep = ep->HENp();
				if( ep->HEVp() != epPrev->HEOp()->HEVp()) 
					return false; // the opposite edge points to a vertex different that the vertex of the next edge
				epPrev = ep;
				if(cnt++ > m.en) 
					return false; // edges are ill connected (HENp())
			}while(ep!=ep1);
		}

	return true;
	}

	/** Set the relations HEFp(), FHEp() from a loop of edges to a face
	*/
	private:
	static void SetRelationsLoopFace(EdgeType * e0, FaceType * f){
		assert(EdgeType::HasHENextAdjacency());
		assert(FaceType::HasFHEAdjacency());

		EdgeType *e = e0;
		assert(e!=NULL);
		do{ e->EFp() = f; e = e->HENp(); } while(e != e0);
		f->FHEp() = e0;
	}

	/**
	Merge the two faces. This will probably become a class template or a functor
	*/
	static void MergeFaces(FaceType *, FaceType *){};

	/**
	Find previous hedge in the loop
	*/
	static EdgeType *  PreviousEdge(EdgeType * e0){
		EdgeType *  ep = e0;
		do{
			if(ep->HENp() == e0) return ep;
				ep = ep->HENp();
			}while(ep!=e0);
	}

	public:
	/** Adds an edge between the sources of e0 and e1 and set all the topology relations.
	If the edges store the pointers to the faces then a new face is created.
    <--- e1 ---- X <------e1_HEPp---
                 ^ 	
                 ||
             ei0 || ei1     
                 ||
                  v
	 ----e0_HEPp-> X ----- e0 ------>
	*/
	static void AddEdge(MeshType &m, EdgeType * e0, EdgeType * e1){
		EdgeType *iii =e0->HENp();
		assert(e1!=e0->HENp());
		assert(e0!=e1->HENp());
		EdgePointer tmp;
		bool hasP =  MeshType::EdgeType::HasHEPrevAdjacency();
		assert(e0->HEOp() != e1); // the hedge already exists
		assert(e0!=e1->HENp());

		std::vector<typename MeshType::EdgePointer* > toUpdate;
		toUpdate.push_back(&e0);
		toUpdate.push_back(&e1);
		EdgeIterator ei0  = vcg::tri::Allocator<MeshType>::AddEdges(m,2,toUpdate);

		EdgeIterator ei1 = ei0; ++ei1;
		(*ei0).HENp() = e1;(*ei0).HEVp() = e0->HEVp();
		(*ei1).HENp() = e0;(*ei1).HEVp() = e1->HEVp();

		EdgePointer e0_HEPp = 0,e1_HEPp = 0,ep =0;
		if(hasP){
			e0_HEPp = e0->HEPp();
			e1_HEPp = e1->HEPp();
		}else{// does not have pointer to previous, it must be computed
			ep = e0;
			do{
				if(ep->HENp() == e0) e0_HEPp = ep;
				if(ep->HENp() == e1) e1_HEPp = ep;
				ep = ep->HENp();
			}while(ep!=e0);
		}
		if(hasP){
			(*ei0).HEPp() = e0->HEPp();
			(*ei1).HEPp() = e1->HEPp();
			e0->HEPp() = &(*ei1);
			e1->HEPp() = &(*ei0);
		}
		e0_HEPp -> HENp() = &(*ei0);
		e1_HEPp -> HENp() = &(*ei1);

		(*ei0).HEOp() = &(*ei1);
		(*ei1).HEOp() = &(*ei0);


		if( EdgeType::HasEFAdjacency() && FaceType::HasFHEAdjacency()){
			FaceIterator fi0  = vcg::tri::Allocator<MeshType>::AddFaces(m,1);
			m.face.back().ImportLocal(*e0->EFp());

			SetRelationsLoopFace(&(*ei0),e1->EFp());		// one loop to the old face
			SetRelationsLoopFace(&(*ei1),&m.face.back());	// the other  to the new face
		}
	}

	/** Detach the topology relations of a given edge
    <--- e->HENPp -X --- <---------eO_HEPp---
                   ^ 	
                   ||
               e   || e->HEOp()     
                   ||
                    v
	 ----e_HEPp--> X ----- e->HEOp->HENPp() ------>
	 
	*/
	static void RemoveEdge(MeshType &m, EdgeType * e){
		assert(MeshType::EdgeType::HasHENextAdjacency());
		assert(MeshType::EdgeType::HasHEOppAdjacency());			
		assert(MeshType::FaceType::HasFHEAdjacency());

		bool hasP =  MeshType::EdgeType::HasHEPrevAdjacency();
		EdgePointer e_HEPp,eO_HEPp;

		if(hasP){
			e_HEPp = e->HEPp();
			eO_HEPp = e->HEOp()->HEPp();
		}else{
			e_HEPp = PreviousEdge(e);
			eO_HEPp = PreviousEdge(e->HEOp());
		}

		assert(e_HEPp->HENp() == e);
		assert(eO_HEPp->HENp() == e->HEOp());
		e_HEPp->HENp() = e->HEOp()->HENp();
		eO_HEPp->HENp() = e-> HENp();

		if(hasP) {
			e->HEOp()->HENp()->HEPp() = e_HEPp;
			e->HENp()->HEPp() = eO_HEPp;

			e->HEPp() = NULL;
			e-> HEOp()->HEPp() = NULL;
		}


		// take care of the faces
		if(MeshType::EdgeType::HasEFAdjacency()){
			MergeFaces(e_HEPp->EFp(),eO_HEPp->EFp());
			vcg::tri::Allocator<MeshType>::DeleteFace(m,*eO_HEPp->EFp());
			SetRelationsLoopFace(e_HEPp,e_HEPp->EFp());

		}
		 vcg::tri::Allocator<MeshType>::DeleteEdge(m,*e->HEOp());
		 vcg::tri::Allocator<MeshType>::DeleteEdge(m,*e);

	}

	};// end class EdgeSupport 
} // end namespace vcg
}
#endif // __VCGLIB_EDGE_SUPPORT
