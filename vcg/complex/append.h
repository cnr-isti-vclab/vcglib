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

#ifndef __VCGLIB_APPEND
#define __VCGLIB_APPEND

#include <vcg/complex/trimesh/allocate.h>
#include <vcg/complex/trimesh/update/flag.h>

namespace vcg {
namespace tri {

template<class MeshLeft, class MeshRight> 
class Append
{
public:
 typedef typename MeshLeft::ScalarType			ScalarLeft;
 typedef typename MeshLeft::CoordType				CoordLeft;
 typedef typename MeshLeft::VertexType			VertexLeft;
 typedef typename MeshLeft::EdgeType				EdgeLeft;
 typedef typename MeshLeft::FaceType				FaceLeft;
 typedef typename MeshLeft::HEdgeType				HEdgeLeft;
 typedef typename MeshLeft::VertexPointer		VertexPointerLeft;
 typedef typename MeshLeft::VertexIterator	VertexIteratorLeft;
 typedef typename MeshLeft::EdgeIterator		EdgeIteratorLeft;
 typedef typename MeshLeft::HEdgeIterator		HEdgeIteratorLeft;
 typedef typename MeshLeft::FaceIterator		FaceIteratorLeft;


 typedef typename MeshRight::ScalarType			ScalarRight;
 typedef typename MeshRight::CoordType			CoordRight;
 typedef typename MeshRight::VertexType			VertexRight;
 typedef typename MeshRight::EdgeType				EdgeRight;
 typedef typename MeshRight::HEdgeType			HEdgeRight;
 typedef typename MeshRight::FaceType				FaceRight;
 typedef typename MeshRight::VertexPointer  VertexPointerRight;
 typedef typename MeshRight::VertexIterator VertexIteratorRight;
 typedef typename MeshRight::EdgeIterator   EdgeIteratorRight;
 typedef typename MeshRight::HEdgeIterator  HEdgeIteratorRight;
 typedef typename MeshRight::FaceIterator   FaceIteratorRight;
 typedef typename MeshRight::FacePointer    FacePointerRight;

 struct Remap{
		std::vector<int> vert,face,edge, hedge;
 };

 static void ImportVertexAdj(MeshLeft &ml, MeshRight &mr, VertexLeft &vl,   VertexRight &vr, Remap &remap, bool sel ){
		 // Vertex to Edge Adj
		 if(vcg::tri::HasVEAdjacency(ml) && vcg::tri::HasVEAdjacency(mr)){
				 size_t i = Index(mr,vr.cVEp());
         vl.VEp() = (i>ml.edge.size())? 0 : &ml.edge[remap.edge[i]];
				 vl.VEi() = vr.VEi();
		 }

		 if(!sel){
			 // Vertex to Face Adj
			 if(vcg::tri::HasPerVertexVFAdjacency(ml) && vcg::tri::HasPerVertexVFAdjacency(mr) &&
					vcg::tri::HasPerFaceVFAdjacency(ml) && vcg::tri::HasPerFaceVFAdjacency(mr)
					){
					size_t i = Index(mr,vr.cVFp());
			vl.VFp() = (i>ml.face.size())? 0 :&ml.face[remap.face[i]];
					vl.VFi() = vr.VFi();
			 }

			 // Vertex to HEdge Adj
			 if(vcg::tri::HasVHAdjacency(ml) && vcg::tri::HasVHAdjacency(mr)){
					vl.VHp() = &ml.hedge[remap.hedge[Index(mr,vr.cVHp())]];
					vl.VHi() = vr.VHi();
			 }
		 }
 }

 static void ImportEdgeAdj(MeshLeft &ml, MeshRight &mr, EdgeLeft &el, const EdgeRight &er, Remap &remap, bool sel ){

		 // Edge to Vertex  Adj
		 if(vcg::tri::HasEVAdjacency(ml) && vcg::tri::HasEVAdjacency(mr)){
				 el.EVp(0) = &ml.vert[remap.vert[Index(mr,er.cEVp(0))]];
				 el.EVp(1) = &ml.vert[remap.vert[Index(mr,er.cEVp(1))]];
		 }

		 if(!sel){
			 // Edge to Edge  Adj
			 if(vcg::tri::HasEEAdjacency(ml) && vcg::tri::HasEEAdjacency(mr))
			 for(unsigned int vi = 0; vi < 2; ++vi)
					 {
				 size_t i = Index(mr,er.cEEp(vi));
				 el.EEp(i) = (i>ml.edge.size())? 0 : &ml.edge[remap.edge[i]];
							 el.EEi(i) = er.cEEi(i);
					 }

			 // Edge to Face  Adj
			 if(vcg::tri::HasEFAdjacency(ml) && vcg::tri::HasEFAdjacency(mr)){
					 size_t i = Index(mr,er.cEFp());
			 el.EFp() = (i>ml.face.size())? 0 :&ml.face[remap.face[i]];
					 el.EFi() = er.cEFi();
			 }

			 // Edge to HEdge   Adj
			 if(vcg::tri::HasEHAdjacency(ml) && vcg::tri::HasEHAdjacency(mr))
					 el.EHp() = &ml.hedge[remap.hedge[Index(mr,er.cEHp())]];
		 }
 }


 static void ImportFaceAdj(MeshLeft &ml, MeshRight &mr, FaceLeft &fl, const FaceRight &fr, Remap &remap, bool sel ){
		 // Face to Vertex  Adj
		 if(vcg::tri::HasFVAdjacency(ml) && vcg::tri::HasFVAdjacency(mr)){
				 assert(fl.VN() == fr.VN());
				 for( int i = 0; i < fl.VN(); ++i ) 
						 fl.V(i) = &ml.vert[remap.vert[Index(mr,fr.V(i))]];
		 }

		 if(!sel){
			 // Face to Edge  Adj
			 if(vcg::tri::HasFEAdjacency(ml) && vcg::tri::HasFEAdjacency(mr)){
					 assert(fl.VN() == fr.VN());
					 for( int vi = 0; vi < fl.VN(); ++vi ){
				 size_t i = Index(mr,fr.cFEp(vi));
				 fl.FEp(i) = (i>ml.edge.size())? 0 : &ml.edge[remap.edge[i]];
					 }
			 }

			 // Face to Face  Adj
			 if(vcg::tri::HasFFAdjacency(ml) && vcg::tri::HasFFAdjacency(mr)){
					 assert(fl.VN() == fr.VN());
					 for( int vi = 0; vi < fl.VN(); ++vi ){
				 size_t i = Index(mr,fr.cFFp(vi));
							 fl.FFp(vi) = (i>ml.face.size()) ? 0 :&ml.face[remap.face[i]];
							 fl.FFi(vi) = fr.cFFi(vi);
					 }
			 }

			 // Face to HEedge  Adj
			 if(vcg::tri::HasFHAdjacency(ml) && vcg::tri::HasFHAdjacency(mr))
							 fl.FHp() = &ml.hedge[remap.hedge[Index(mr,fr.cFHp())]];
		 }
 }

 static void ImportHEdgeAdj(MeshLeft &ml, MeshRight &mr, HEdgeLeft &hl, const HEdgeRight &hr, Remap &remap, bool /*sel*/ ){
		 // HEdge to Vertex  Adj
		 if(vcg::tri::HasHVAdjacency(ml) && vcg::tri::HasHVAdjacency(mr))
					hl.HVp() = &ml.vert[remap.vert[Index(mr,hr.cHVp())]];

			// HEdge to Edge  Adj
		 if(vcg::tri::HasHEAdjacency(ml) && vcg::tri::HasHEAdjacency(mr)){
					size_t i = Index(mr,hr.cHEp()) ;
          hl.HEp() = (i>ml.edge.size())? 0 : &ml.edge[remap.edge[i]];
			}

		 // HEdge to Face  Adj
		 if(vcg::tri::HasHFAdjacency(ml) && vcg::tri::HasHFAdjacency(mr)){
				 size_t i = Index(mr,hr.cHFp());
         hl.HFp() = (i>ml.face.size())? 0 :&ml.face[remap.face[i]];
			}


		 // HEdge to Opposite HEdge  Adj
		 if(vcg::tri::HasHOppAdjacency(ml) && vcg::tri::HasHOppAdjacency(mr))
					hl.HOp() = &ml.hedge[remap.hedge[Index(mr,hr.cHOp())]];

		 // HEdge to Next  HEdge  Adj
		 if(vcg::tri::HasHNextAdjacency(ml) && vcg::tri::HasHNextAdjacency(mr))
					hl.HNp() = &ml.hedge[remap.hedge[Index(mr,hr.cHNp())]];

		 // HEdge to Next  HEdge  Adj
		 if(vcg::tri::HasHPrevAdjacency(ml) && vcg::tri::HasHPrevAdjacency(mr))
					hl.HPp() = &ml.hedge[remap.hedge[Index(mr,hr.cHPp())]];
 }

// Append Right Mesh to the Left Mesh
// Append::Mesh(ml, mr) is equivalent to ml += mr. 
// Note MeshRigth could be costant...

static void Mesh(MeshLeft& ml, MeshRight& mr, const bool selected = false){

		// phase 1. allocate on ml vert,edge,face, hedge to accomodat those of mr
		// and build the remapping for all

		Remap remap;

		// vertex
		remap.vert.resize(mr.vert.size(),-1);
		VertexIteratorRight vi;
		for(vi=mr.vert.begin();vi!=mr.vert.end();++vi)
				if(!(*vi).IsD() && (!selected || (*vi).IsS())){
						int ind=Index(mr,*vi);
						assert(remap.vert[ind]==-1);
						VertexIteratorLeft vp;
						vp=Allocator<MeshLeft>::AddVertices(ml,1);
						(*vp).ImportData(*(vi));
						remap.vert[ind]=Index(ml,*vp);
					}

		// edge
		remap.edge.resize(mr.edge.size(),-1);
		EdgeIteratorRight ei;
		for(ei=mr.edge.begin(); ei!=mr.edge.end();++ei)
				if(!(*ei).IsD() && (!selected || (*ei).IsS())){
						int ind=Index(mr,*ei);
						assert(remap.edge[ind]==-1);
						EdgeIteratorLeft ep;
						ep=Allocator<MeshLeft>::AddEdges(ml,1);
						(*ep).ImportData(*(ei));
						remap.edge[ind]=Index(ml,*ep);
				}

		// face
		vcg::tri::Allocator<MeshRight>::CompactFaceVector(mr);
		remap.face.resize(mr.face.size(),-1);
		FaceIteratorRight fi;
		for(fi=mr.face.begin();fi!=mr.face.end();++fi)
				if(!(*fi).IsD() && (!selected || (*fi).IsS())){
						int ind=Index(mr,*fi);
						assert(remap.face[ind]==-1);
						FaceIteratorLeft fp;
						fp=Allocator<MeshLeft>::AddFaces(ml,1);
						(*fp).ImportData(*(fi));
						remap.face[ind]=Index(ml,*fp);
				}

		// hedge
		remap.hedge.resize(mr.hedge.size(),-1);
		HEdgeIteratorRight hi;
		for(hi=mr.hedge.begin();hi!=mr.hedge.end();++hi)
				if(!(*hi).IsD() && (!selected || (*hi).IsS())){
						int ind=Index(mr,*hi);
						assert(remap.hedge[ind]==-1); 
						HEdgeIteratorLeft hp;
						hp=Allocator<MeshLeft>::AddHEdges(ml,1);
						(*hp).ImportData(*(hi));
						remap.hedge[ind]=Index(ml,*hp);
				}

		// phase 2.
		// copy data from ml to its corresponding elements in ml and adjacencies

		// vertex
		for(vi=mr.vert.begin();vi!=mr.vert.end();++vi)
				if( !(*vi).IsD() && (!selected || (*vi).IsS())){
				 ml.vert[remap.vert[Index(mr,*vi)]].ImportData(*vi);
				 ImportVertexAdj(ml,mr,ml.vert[remap.vert[Index(mr,*vi)]],*vi,remap,selected);
		}

		// edge
		for(ei=mr.edge.begin();ei!=mr.edge.end();++ei)
				if(!(*ei).IsD() && (!selected || (*ei).IsS())){
				 ml.edge[remap.edge[Index(mr,*ei)]].ImportData(*ei);
				 ImportEdgeAdj(ml,mr,ml.edge[remap.edge[Index(mr,*ei)]],*ei,remap,selected);
		}

		// face
		bool wedgetexcoord = vcg::tri::HasPerWedgeTexCoord(mr); 
		for(fi=mr.face.begin();fi!=mr.face.end();++fi)
				if(!(*fi).IsD() && (!selected || (*fi).IsS())){
				 if(wedgetexcoord)
					 for(int i = 0; i < (*fi).VN(); ++i)
						 (*fi).WT(i).n() += ml.textures.size();
				 ml.face[remap.face[Index(mr,*fi)]].ImportData(*fi);
				 ImportFaceAdj(ml,mr,ml.face[remap.face[Index(mr,*fi)]],*fi,remap,selected);

		}

		// hedge
		for(hi=mr.hedge.begin();hi!=mr.hedge.end();++hi)
				if(!(*hi).IsD() && (!selected || (*hi).IsS())){
				 ml.hedge[remap.hedge[Index(mr,*hi)]].ImportData(*hi);
				 ImportHEdgeAdj(ml,mr,ml.hedge[remap.hedge[Index(mr,*hi)]],*hi,remap,selected);
		}

		// phase 3.
		// take care of other per mesh data: textures,   attributes

		// At the end concatenate the vector with texture names.
		ml.textures.insert(ml.textures.end(),mr.textures.begin(),mr.textures.end());

		// Attributes. Copy only those attributes that are present in both meshes
		// Two attributes in different meshes are considered the same if they have the same
		// name and the same type. This may be deceiving because they could in fact have
		// different semantic, but this is up to the developer.
		// If the left mesh has attributes that are not in the right mesh, their values for the elements
		// of the right mesh will be uninitialized

		// to be done.
		// note: we need to assign attribute values without knowing their type



}


static void Selected(MeshLeft& ml, MeshRight& mr)
{
  Mesh(ml,mr,true);
}

}; // end of class Append





} // End Namespace TriMesh
} // End Namespace vcg


#endif


