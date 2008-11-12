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
		struct EdgeSupport{
			typedef typename MeshType::VertexType VertexType;
			typedef typename MeshType::VertexPointer VertexPointer;
			typedef typename MeshType::EdgePointer EdgePointer;
			typedef typename  MeshType::FaceIterator FaceIterator;

			struct VertexPairEdgePtr{
				VertexPairEdgePtr(VertexPointer _v0,VertexPointer _v1,EdgePointer _ep):v0(_v0),v1(_v1),ep(_ep){if(v0>v1) std::swap(v0,v1);}
				const bool operator <(const VertexPairEdgePtr &o){return (v0 == o.v0)? (v1<o.v1):(v0<o.v0);}
				const bool operator ==(const VertexPairEdgePtr &o){return (v0 == o.v0)&& (v1==o.v1);}

				VertexPointer v0,v1;
				EdgePointer ep;
			};
			/**
			build a half-edge data structure from an indexed data structure. Note that the half-edges are allocated here for the first time.
			If you have a mesh where there are already edges, they will be removed and the data lost, so do not use this function
			to just "update" the topology of half edges.
			**/
			static void ComputeHalfEdgeFromIndexed(MeshType & m){
				assert(HasFVAdjacency(m));
				assert(MeshType::EdgeType::HasHENextAdjacency());
				assert(MeshType::EdgeType::HasHEOppAdjacency());			

				// allocate all new half edges
				FaceIterator fi;
				int n_edges = 0;
				for(fi = m.face.begin(); fi != m.face.end(); ++fi) if(! (*fi).IsD()) n_edges+=(*fi).VN();		
				MeshType::EdgeIterator ei = vcg::tri::Allocator<MeshType>::AddEdges(m,n_edges);
				
				std::vector<VertexPairEdgePtr> all;
				int firstEdge = 0;
				for(fi = m.face.begin(); fi != m.face.end(); ++fi)if(!(*fi).IsD()){
					assert((*fi).VN()>2);
					for(int i  = 0; i < (*fi).VN(); ++i,++ei)
					{
						(*ei).HEVp()	= (*fi).V(i);
						(*ei).HENp()	= &m.edge[firstEdge + (i +1) % (*fi).VN()];
						if(MeshType::EdgeType::HasEFAdjacency()) 
							(*ei).EFp() = &(*fi);
						if(MeshType::EdgeType::HasHEPrevAdjacency())
							(*ei).HEPp()	= &m.edge[firstEdge + (i +(*fi).VN()-1) % (*fi).VN()];
						all.push_back(VertexPairEdgePtr((*fi).V(i), (*fi).V((*fi).Next(i)),&(*ei)));// it will be used to link the hedges
					}
					firstEdge += (*fi).VN();
				}
				

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
				assert(MeshType::EdgeType::HasHEOppAdjacency());			
				bool createFace,hasHEF;

				typename MeshType::PerEdgeAttributeHandle<bool> hV = Allocator<MeshType>::AddPerEdgeAttribute<bool>(m,"");

				typename MeshType::EdgeIterator ei;
				typename MeshType::FacePointer fp;
				typename MeshType::FaceIterator fi;
				typename MeshType::EdgePointer ep,epF;
				int vi = 0;

				hasHEF = (typename MeshType::EdgeType::HasEFAdjacency());
				assert( !hasHEF || (hasHEF && m.fn>0));

				// if the edgetype has the pointer to face   
				// it is assumed the the edget2face pointer (HEFp) are correct
				// and the faces are allocated
					for ( ei = m.edge.begin(); ei != m.edge.end(); ++ei)
 					if(!hV[(*ei)] )// has not be visited yet
					{
						if(!hasHEF)// if it has 
							fp = &(* Allocator<MeshType>::AddFaces(m,1));
						else
							fp = (*ei).EFp();

						ep = epF = &(*ei);
						ep = ep->HENp();
						std::vector<VertexPointer> vpts;
						while(ep!=epF){vpts.push_back((*ep).HEVp()); ep=ep->HENp();}
						fp ->Alloc(vpts.size());
						for(int i  = 0; i < vpts.size();++i) fp ->V(i) = vpts[i];// set the pointer from face to vertex

 					hV[(*ei)] = true;
				}
				Allocator<MeshType>::DeletePerEdgeAttribute(m,hV);
			}
	};
} // end namespace vcg
}
#endif // __VCGLIB_EDGE_SUPPORT
