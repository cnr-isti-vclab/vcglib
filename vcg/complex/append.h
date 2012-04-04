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

#include <vcg/complex/allocate.h>
#include <vcg/complex/algorithms/update/flag.h>
#include <vcg/complex/algorithms/update/selection.h>
#include <set>

namespace vcg {
namespace tri {

template<class MeshLeft, class ConstMeshRight>
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


 typedef typename ConstMeshRight::ScalarType			ScalarRight;
 typedef typename ConstMeshRight::CoordType			CoordRight;
 typedef typename ConstMeshRight::VertexType			VertexRight;
 typedef typename ConstMeshRight::EdgeType				EdgeRight;
 typedef typename ConstMeshRight::HEdgeType			HEdgeRight;
 typedef typename ConstMeshRight::FaceType				FaceRight;
 typedef typename ConstMeshRight::VertexPointer  VertexPointerRight;
 typedef typename ConstMeshRight::VertexIterator VertexIteratorRight;
 typedef typename ConstMeshRight::EdgeIterator   EdgeIteratorRight;
 typedef typename ConstMeshRight::HEdgeIterator  HEdgeIteratorRight;
 typedef typename ConstMeshRight::FaceIterator   FaceIteratorRight;
 typedef typename ConstMeshRight::FacePointer    FacePointerRight;

 struct Remap{
		std::vector<int> vert,face,edge, hedge;
 };

 static void ImportVertexAdj(MeshLeft &ml, ConstMeshRight &mr, VertexLeft &vl,   VertexRight &vr, Remap &remap ){
   // Vertex to Edge Adj
   if(vcg::tri::HasVEAdjacency(ml) && vcg::tri::HasVEAdjacency(mr) && vr.cVEp() != 0){
     size_t i = Index(mr,vr.cVEp());
     vl.VEp() = (i>ml.edge.size())? 0 : &ml.edge[remap.edge[i]];
     vl.VEi() = vr.VEi();
   }

   // Vertex to Face Adj
   if(vcg::tri::HasPerVertexVFAdjacency(ml) && vcg::tri::HasPerVertexVFAdjacency(mr) && vr.cVFp() != 0 ){
     size_t i = Index(mr,vr.cVFp());
     vl.VFp() = (i>ml.face.size())? 0 :&ml.face[remap.face[i]];
     vl.VFi() = vr.VFi();
   }

   // Vertex to HEdge Adj
   if(vcg::tri::HasVHAdjacency(ml) && vcg::tri::HasVHAdjacency(mr) && vr.cVHp() != 0){
     vl.VHp() = &ml.hedge[remap.hedge[Index(mr,vr.cVHp())]];
     vl.VHi() = vr.VHi();
   }
 }

 static void ImportEdgeAdj(MeshLeft &ml, ConstMeshRight &mr, EdgeLeft &el, const EdgeRight &er, Remap &remap)
 {
     // Edge to Edge  Adj
     if(vcg::tri::HasEEAdjacency(ml) && vcg::tri::HasEEAdjacency(mr))
       for(unsigned int vi = 0; vi < 2; ++vi)
       {
         size_t idx = Index(mr,er.cEEp(vi));
         el.EEp(vi) = (idx>ml.edge.size())? 0 : &ml.edge[remap.edge[idx]];
         el.EEi(vi) = er.cEEi(vi);
       }

     // Edge to Face  Adj
     if(vcg::tri::HasEFAdjacency(ml) && vcg::tri::HasEFAdjacency(mr)){
       size_t idx = Index(mr,er.cEFp());
       el.EFp() = (idx>ml.face.size())? 0 :&ml.face[remap.face[idx]];
       el.EFi() = er.cEFi();
     }

     // Edge to HEdge   Adj
     if(vcg::tri::HasEHAdjacency(ml) && vcg::tri::HasEHAdjacency(mr))
       el.EHp() = &ml.hedge[remap.hedge[Index(mr,er.cEHp())]];
 }


 static void ImportFaceAdj(MeshLeft &ml, ConstMeshRight &mr, FaceLeft &fl, const FaceRight &fr, Remap &remap )
 {
     // Face to Edge  Adj
     if(vcg::tri::HasFEAdjacency(ml) && vcg::tri::HasFEAdjacency(mr)){
       assert(fl.VN() == fr.VN());
       for( int vi = 0; vi < fl.VN(); ++vi ){
         size_t idx = Index(mr,fr.cFEp(vi));
         fl.FEp(vi) = (idx>ml.edge.size())? 0 : &ml.edge[remap.edge[idx]];
       }
     }

     // Face to Face  Adj
     if(vcg::tri::HasFFAdjacency(ml) && vcg::tri::HasFFAdjacency(mr)){
       assert(fl.VN() == fr.VN());
       for( int vi = 0; vi < fl.VN(); ++vi ){
         size_t idx = Index(mr,fr.cFFp(vi));
         fl.FFp(vi) = (idx>ml.face.size()) ? 0 :&ml.face[remap.face[idx]];
         fl.FFi(vi) = fr.cFFi(vi);
       }
     }

     // Face to HEedge  Adj
     if(vcg::tri::HasFHAdjacency(ml) && vcg::tri::HasFHAdjacency(mr))
       fl.FHp() = &ml.hedge[remap.hedge[Index(mr,fr.cFHp())]];
 }

 static void ImportHEdgeAdj(MeshLeft &ml, ConstMeshRight &mr, HEdgeLeft &hl, const HEdgeRight &hr, Remap &remap, bool /*sel*/ ){
   // HEdge to Vertex  Adj
   if(vcg::tri::HasHVAdjacency(ml) && vcg::tri::HasHVAdjacency(mr))
     hl.HVp() = &ml.vert[remap.vert[Index(mr,hr.cHVp())]];

   // HEdge to Edge  Adj
   if(vcg::tri::HasHEAdjacency(ml) && vcg::tri::HasHEAdjacency(mr)){
     size_t idx = Index(mr,hr.cHEp()) ;
     hl.HEp() = (idx>ml.edge.size())? 0 : &ml.edge[remap.edge[idx]];
   }

   // HEdge to Face  Adj
   if(vcg::tri::HasHFAdjacency(ml) && vcg::tri::HasHFAdjacency(mr)){
     size_t idx = Index(mr,hr.cHFp());
     hl.HFp() = (idx>ml.face.size())? 0 :&ml.face[remap.face[idx]];
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

static void Mesh(MeshLeft& ml, ConstMeshRight& mr, const bool selected = false, const bool adjFlag = false)
{
  // Note that if the the selection of the vertexes is not consistent with the face selection
  // the append could build faces referencing non existent vertices
  // so it is mandatory that the selection of the vertices reflects the loose selection
  // from edges and faces (e.g. if a face is selected all its vertices must be selected).
  if(selected)
  {
    assert(adjFlag == false); // It is rather meaningless to partially copy adj relations.
    tri::UpdateSelection<ConstMeshRight>::VertexFromEdgeLoose(mr,true);
    tri::UpdateSelection<ConstMeshRight>::VertexFromFaceLoose(mr,true);
  }

  // phase 1. allocate on ml vert,edge,face, hedge to accomodat those of mr
  // and build the remapping for all

  Remap remap;

  // vertex
  remap.vert.resize(mr.vert.size(),-1);
  VertexIteratorRight vi;
  VertexIteratorLeft vp;
  int svn = UpdateSelection<ConstMeshRight>::VertexCount(mr);
  if(selected) vp=Allocator<MeshLeft>::AddVertices(ml,svn);
          else vp=Allocator<MeshLeft>::AddVertices(ml,mr.vn);

  for(vi=mr.vert.begin();vi!=mr.vert.end();++vi)
    if(!(*vi).IsD() && (!selected || (*vi).IsS())){
      int ind=Index(mr,*vi);
      remap.vert[ind]=Index(ml,*vp);
      ++vp;
    }

  // edge
  remap.edge.resize(mr.edge.size(),-1);
  EdgeIteratorRight ei;
  EdgeIteratorLeft ep;
  int sen = UpdateSelection<ConstMeshRight>::EdgeCount(mr);
  if(selected) ep=Allocator<MeshLeft>::AddEdges(ml,sen);
          else ep=Allocator<MeshLeft>::AddEdges(ml,mr.en);

  for(ei=mr.edge.begin(); ei!=mr.edge.end();++ei)
    if(!(*ei).IsD() && (!selected || (*ei).IsS())){
      int ind=Index(mr,*ei);
      remap.edge[ind]=Index(ml,*ep);
      ++ep;
    }

  // face
  remap.face.resize(mr.face.size(),-1);
  FaceIteratorRight fi;
  FaceIteratorLeft fp;
  int sfn = UpdateSelection<ConstMeshRight>::FaceCount(mr);
  if(selected) fp=Allocator<MeshLeft>::AddFaces(ml,sfn);
          else fp=Allocator<MeshLeft>::AddFaces(ml,mr.fn);

  for(fi=mr.face.begin();fi!=mr.face.end();++fi)
    if(!(*fi).IsD() && (!selected || (*fi).IsS())){
      int ind=Index(mr,*fi);
      remap.face[ind]=Index(ml,*fp);
      ++fp;
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
      if(adjFlag) ImportVertexAdj(ml,mr,ml.vert[remap.vert[Index(mr,*vi)]],*vi,remap);
    }

  // edge
  for(ei=mr.edge.begin();ei!=mr.edge.end();++ei)
    if(!(*ei).IsD() && (!selected || (*ei).IsS())){
      ml.edge[remap.edge[Index(mr,*ei)]].ImportData(*ei);
      // Edge to Vertex  Adj
      EdgeLeft &el = ml.edge[remap.edge[Index(mr,*ei)]];
      if(vcg::tri::HasEVAdjacency(ml) && vcg::tri::HasEVAdjacency(mr)){
        el.V(0) = &ml.vert[remap.vert[Index(mr,ei->cV(0))]];
        el.V(1) = &ml.vert[remap.vert[Index(mr,ei->cV(1))]];
      }
      if(adjFlag) ImportEdgeAdj(ml,mr,el,*ei,remap);
    }

  // face
  const int textureOffset =  ml.textures.size();
  bool WTFlag = vcg::tri::HasPerWedgeTexCoord(mr) && (textureOffset>0);
  for(fi=mr.face.begin();fi!=mr.face.end();++fi)
    if(!(*fi).IsD() && (!selected || (*fi).IsS()))
    {
      FaceLeft &fl = ml.face[remap.face[Index(mr,*fi)]];
      if(WTFlag)
        for(int i = 0; i < 3; ++i)
          fl.WT(i).n() +=textureOffset;
      if(vcg::tri::HasFVAdjacency(ml) && vcg::tri::HasFVAdjacency(mr)){
        fl.V(0) = &ml.vert[remap.vert[Index(mr,fi->cV(0))]];
        fl.V(1) = &ml.vert[remap.vert[Index(mr,fi->cV(1))]];
        fl.V(2) = &ml.vert[remap.vert[Index(mr,fi->cV(2))]];
      }
      ml.face[remap.face[Index(mr,*fi)]].ImportData(*fi);
      if(adjFlag)  ImportFaceAdj(ml,mr,ml.face[remap.face[Index(mr,*fi)]],*fi,remap);

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

		unsigned int id_r;
		typename std::set< PointerToAttribute  >::iterator al, ar;

		// per vertex attributes
		for(al = ml.vert_attr.begin(); al != ml.vert_attr.end(); ++al)
			if(!(*al)._name.empty()){
				ar =    mr.vert_attr.find(*al);
				if(ar!= mr.vert_attr.end()){
					id_r = 0;
					for(vi=mr.vert.begin();vi!=mr.vert.end();++vi,++id_r)
						if( !(*vi).IsD() && (!selected || (*vi).IsS()))
							memcpy((*al)._handle->At(remap.vert[Index(mr,*vi)]),(*ar)._handle->At(id_r),
								(*al)._handle->SizeOf());
				}
			}

		// per edge attributes
		for(al = ml.edge_attr.begin(); al != ml.edge_attr.end(); ++al)
			if(!(*al)._name.empty()){
				ar =    mr.edge_attr.find(*al);
				if(ar!= mr.edge_attr.end()){
					id_r = 0;
					for(ei=mr.edge.begin();ei!=mr.edge.end();++ei,++id_r)
						if( !(*ei).IsD() && (!selected || (*ei).IsS()))
							memcpy((*al)._handle->At(remap.vert[Index(mr,*ei)]),(*ar)._handle->At(id_r),
								(*al)._handle->SizeOf());
				}
			}

		// per face attributes
		for(al = ml.face_attr.begin(); al != ml.face_attr.end(); ++al)
			if(!(*al)._name.empty()){
				ar =    mr.face_attr.find(*al);
				if(ar!= mr.face_attr.end()){
					id_r = 0;
					for(fi=mr.face.begin();fi!=mr.face.end();++fi,++id_r)
						if( !(*fi).IsD() && (!selected || (*fi).IsS()))
							memcpy((*al)._handle->At(remap.vert[Index(mr,*fi)]),(*ar)._handle->At(id_r),
								(*al)._handle->SizeOf());
				}
			}

                // per mesh attributes
                // if both ml and mr have an attribute with the same name, no action is done
                // if mr has an attribute that is NOT present in ml, the attribute is added to ml
                //for(ar = mr.mesh_attr.begin(); ar != mr.mesh_attr.end(); ++ar)
                //        if(!(*ar)._name.empty()){
                //                al =    ml.mesh_attr.find(*ar);
                //                if(al== ml.mesh_attr.end())
                //                        //...
                //        }
}

static void MeshCopy(MeshLeft& ml, ConstMeshRight& mr, bool selected=false)
{
  ml.Clear();
  Mesh(ml,mr,selected);
  ml.bbox=mr.bbox;
}

static void Selected(MeshLeft& ml, ConstMeshRight& mr)
{
  Mesh(ml,mr,true);
}

}; // end of class Append





} // End Namespace tri
} // End Namespace vcg


#endif


