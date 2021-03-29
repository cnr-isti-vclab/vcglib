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
#ifndef __VCGLIB_APPEND
#define __VCGLIB_APPEND

#include <vcg/complex/allocate.h>
#include <vcg/complex/algorithms/update/selection.h>

namespace vcg {
namespace tri {
/** \ingroup trimesh */
/*! \brief  Class to safely duplicate and append (portion of) meshes.

Adding elements to a mesh, like faces and vertices can involve the reallocation of the vectors of the involved elements.
This class provide the only safe methods to add elements of a mesh to another one.
\sa \ref allocation
*/
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
 typedef typename MeshLeft::TetraType       TetraLeft;
 typedef typename MeshLeft::VertexPointer		VertexPointerLeft;
 typedef typename MeshLeft::VertexIterator	VertexIteratorLeft;
 typedef typename MeshLeft::EdgeIterator		EdgeIteratorLeft;
 typedef typename MeshLeft::HEdgeIterator		HEdgeIteratorLeft;
 typedef typename MeshLeft::FaceIterator		FaceIteratorLeft;
 typedef typename MeshLeft::TetraIterator   TetraIteratorLeft;


 typedef typename ConstMeshRight::ScalarType     ScalarRight;
 typedef typename ConstMeshRight::CoordType      CoordRight;
 typedef typename ConstMeshRight::VertexType     VertexRight;
 typedef typename ConstMeshRight::EdgeType       EdgeRight;
 typedef typename ConstMeshRight::HEdgeType      HEdgeRight;
 typedef typename ConstMeshRight::FaceType       FaceRight;
 typedef typename ConstMeshRight::TetraType      TetraRight;
 typedef typename ConstMeshRight::TetraPointer   TetraPointerRight;
 typedef typename ConstMeshRight::TetraIterator  TetraIteratorRight;
 typedef typename ConstMeshRight::VertexPointer  VertexPointerRight;
 typedef typename ConstMeshRight::VertexIterator VertexIteratorRight;
 typedef typename ConstMeshRight::EdgeIterator   EdgeIteratorRight;
 typedef typename ConstMeshRight::HEdgeIterator  HEdgeIteratorRight;
 typedef typename ConstMeshRight::FaceIterator   FaceIteratorRight;
 typedef typename ConstMeshRight::FacePointer    FacePointerRight;

 struct Remap{
   static size_t InvalidIndex() { return std::numeric_limits<size_t>::max(); }
        std::vector<size_t> vert, face, edge, hedge, tetra;
 };

 static void ImportVertexAdj(MeshLeft &ml, const ConstMeshRight &mr, VertexLeft &vl,   const VertexRight &vr, Remap &remap ){
   // Vertex to Edge Adj
   if(HasVEAdjacency(ml) && HasVEAdjacency(mr) && vr.cVEp() != 0){
     size_t i = Index(mr,vr.cVEp());
     vl.VEp() = (i>ml.edge.size())? 0 : &ml.edge[remap.edge[i]];
     vl.VEi() = vr.VEi();
   }

   // Vertex to Face Adj
   if(HasPerVertexVFAdjacency(ml) && HasPerVertexVFAdjacency(mr) && vr.cVFp() != 0 ){
     size_t i = Index(mr,vr.cVFp());
     vl.VFp() = (i>ml.face.size())? 0 :&ml.face[remap.face[i]];
     vl.VFi() = vr.VFi();
   }

   // Vertex to HEdge Adj
   if(HasVHAdjacency(ml) && HasVHAdjacency(mr) && vr.cVHp() != 0){
     vl.VHp() = &ml.hedge[remap.hedge[Index(mr,vr.cVHp())]];
     vl.VHi() = vr.VHi();
   }

   // Vertex to Tetra Adj
   if(HasVTAdjacency(ml) && HasVTAdjacency(mr) && vr.cVTp() != 0){
     size_t i = Index(mr, vr.cVTp());
     vl.VTp() = (i > ml.edge.size()) ? 0 : &ml.tetra[remap.tetra[i]];
     vl.VTi() = vr.VTi();
   }
 }

 static void ImportEdgeAdj(MeshLeft &ml, const ConstMeshRight &mr, EdgeLeft &el, const EdgeRight &er, Remap &remap)
 {
     // Edge to Edge  Adj
     if(HasEEAdjacency(ml) && HasEEAdjacency(mr))
       for(unsigned int vi = 0; vi < 2; ++vi)
       {
         size_t idx = Index(mr,er.cEEp(vi));
         el.EEp(vi) = (idx>ml.edge.size())? 0 : &ml.edge[remap.edge[idx]];
         el.EEi(vi) = er.cEEi(vi);
       }

     // Edge to Face  Adj
     if(HasEFAdjacency(ml) && HasEFAdjacency(mr)){
       size_t idx = Index(mr,er.cEFp());
       el.EFp() = (idx>ml.face.size())? 0 :&ml.face[remap.face[idx]];
       el.EFi() = er.cEFi();
     }

     // Edge to HEdge   Adj
     if(HasEHAdjacency(ml) && HasEHAdjacency(mr))
       el.EHp() = &ml.hedge[remap.hedge[Index(mr,er.cEHp())]];
 }


 static void ImportFaceAdj(MeshLeft &ml, const ConstMeshRight &mr, FaceLeft &fl, const FaceRight &fr, Remap &remap )
 {
   // Face to Edge  Adj
   if(HasFEAdjacency(ml) && HasFEAdjacency(mr)){
     assert(fl.VN() == fr.VN());
     for( int vi = 0; vi < fl.VN(); ++vi ){
       size_t idx = remap.edge[Index(mr,fr.cFEp(vi))];
       if(idx!=Remap::InvalidIndex())
         fl.FEp(vi) = &ml.edge[idx];
     }
   }

   // Face to Face  Adj
   if(HasFFAdjacency(ml) && HasFFAdjacency(mr)){
     assert(fl.VN() == fr.VN());
     for( int vi = 0; vi < fl.VN(); ++vi ){
       size_t idx = remap.face[Index(mr,fr.cFFp(vi))];
       if(idx!=Remap::InvalidIndex()){
		 assert(idx >= 0 && idx < ml.face.size());
         fl.FFp(vi) = &ml.face[idx];
         fl.FFi(vi) = fr.cFFi(vi);
       }
     }
   }

	// Vertex to Face Adj
	if(HasPerFaceVFAdjacency(ml) && HasPerFaceVFAdjacency(mr))
	{
		assert(fl.VN() == fr.VN());
		for (int vi = 0; vi < fl.VN(); ++vi)
		{
			const auto * fp     = fr.cVFp(vi);
			const auto   vfindex = fr.cVFi(vi);
			size_t fidx = (fp == nullptr) ? Remap::InvalidIndex() : remap.face[Index(mr,fp)];

			if (fidx == Remap::InvalidIndex()) // end of VF chain (or not initialized)
			{
				fl.VFClear(vi);
				assert(fl.cVFi(vi) == -1);
			}
			else
			{
				assert(fidx >= 0 && fidx < ml.face.size());
				fl.VFp(vi) = &ml.face[fidx];
				fl.VFi(vi) = vfindex;
			}
		}
	}

   // Face to HEedge  Adj
   if(HasFHAdjacency(ml) && HasFHAdjacency(mr))
     fl.FHp() = &ml.hedge[remap.hedge[Index(mr,fr.cFHp())]];
 }

 static void ImportHEdgeAdj(MeshLeft &ml, const ConstMeshRight &mr, HEdgeLeft &hl, const HEdgeRight &hr, Remap &remap, bool /*sel*/ ){
   // HEdge to Vertex  Adj
   if(HasHVAdjacency(ml) && HasHVAdjacency(mr))
     hl.HVp() = &ml.vert[remap.vert[Index(mr,hr.cHVp())]];

   // HEdge to Edge  Adj
   if(HasHEAdjacency(ml) && HasHEAdjacency(mr)){
     size_t idx = Index(mr,hr.cHEp()) ;
     hl.HEp() = (idx>ml.edge.size())? 0 : &ml.edge[remap.edge[idx]];
   }

   // HEdge to Face  Adj
   if(HasHFAdjacency(ml) && HasHFAdjacency(mr)){
     size_t idx = Index(mr,hr.cHFp());
     hl.HFp() = (idx>ml.face.size())? 0 :&ml.face[remap.face[idx]];
   }


   // HEdge to Opposite HEdge  Adj
   if(HasHOppAdjacency(ml) && HasHOppAdjacency(mr))
     hl.HOp() = &ml.hedge[remap.hedge[Index(mr,hr.cHOp())]];

   // HEdge to Next  HEdge  Adj
   if(HasHNextAdjacency(ml) && HasHNextAdjacency(mr))
     hl.HNp() = &ml.hedge[remap.hedge[Index(mr,hr.cHNp())]];

   // HEdge to Next  HEdge  Adj
   if(HasHPrevAdjacency(ml) && HasHPrevAdjacency(mr))
     hl.HPp() = &ml.hedge[remap.hedge[Index(mr,hr.cHPp())]];
 }

  static void ImportTetraAdj(MeshLeft &ml, const ConstMeshRight &mr, TetraLeft &tl, const TetraRight &tr, Remap &remap )
 {
   // Tetra to Tetra  Adj
   if(HasTTAdjacency(ml) && HasTTAdjacency(mr)){
     for( int vi = 0; vi < 4; ++vi ){
       size_t idx = remap.tetra[Index(mr,tr.cTTp(vi))];
       if(idx != Remap::InvalidIndex()){
         tl.TTp(vi) = &ml.tetra[idx];
         tl.TTi(vi) = tr.cTTi(vi);
       }
     }
   }
}

// Append Right Mesh to the Left Mesh
// Append::Mesh(ml, mr) is equivalent to ml += mr.
// Note MeshRigth could be costant...
 /*! \brief %Append the second mesh to the first one.

   The first mesh is not destroyed and no attempt of avoid duplication of already present elements is done.
   If requested only the selected elements are appended to the first one.
   The second mesh is not changed at all (it could be constant) with the exception of the selection (see below note).

   \note  If the the selection of the vertexes is not consistent with the face selection
   the append could build faces referencing non existent vertices
   so it is mandatory that the selection of the vertices reflects the loose selection
   from edges and faces (e.g. if a face is selected then all its vertices must be selected).

   \note Attributes. This function will copy only those attributes that are present in both meshes.
   Two attributes in different meshes are considered the same iff they have the same
   name and the same type. This may be deceiving because they could in fact have
   different semantic, but this is up to the developer.
   If the left mesh has attributes that are not in the right mesh, their values for the elements
   of the right mesh will be uninitialized

 */

static void Mesh(MeshLeft& ml, ConstMeshRight& mr, const bool selected = false, const bool adjFlag = false)
{
  // Note that if the the selection of the vertexes is not consistent with the face selection
  // the append could build faces referencing non existent vertices
  // so it is mandatory that the selection of the vertices reflects the loose selection
  // from edges and faces (e.g. if a face is selected all its vertices must be selected).
  // note the use of the parameter for preserving existing vertex selection.
  if(selected)
  {
    assert(adjFlag == false || ml.IsEmpty()); // It is rather meaningless to partially copy adj relations.
    tri::UpdateSelection<ConstMeshRight>::VertexFromEdgeLoose(mr,true);
    tri::UpdateSelection<ConstMeshRight>::VertexFromFaceLoose(mr,true);
  }

  // phase 1. allocate on ml vert,edge,face, hedge to accomodat those of mr
  // and build the remapping for all

  Remap remap;

  // vertex
  remap.vert.resize(mr.vert.size(), Remap::InvalidIndex());
  VertexIteratorLeft vp;
  size_t svn = UpdateSelection<ConstMeshRight>::VertexCount(mr);
  if(selected)
      vp=Allocator<MeshLeft>::AddVertices(ml,int(svn));
  else
      vp=Allocator<MeshLeft>::AddVertices(ml,mr.vn);

  for(VertexIteratorRight vi=mr.vert.begin(); vi!=mr.vert.end(); ++vi)
  {
    if(!(*vi).IsD() && (!selected || (*vi).IsS()))
    {
      size_t ind=Index(mr,*vi);
      remap.vert[ind]=int(Index(ml,*vp));
      ++vp;
    }
  }
  // edge
  remap.edge.resize(mr.edge.size(), Remap::InvalidIndex());
  EdgeIteratorLeft ep;
  size_t sen = UpdateSelection<ConstMeshRight>::EdgeCount(mr);
  if(selected) ep=Allocator<MeshLeft>::AddEdges(ml,sen);
          else ep=Allocator<MeshLeft>::AddEdges(ml,mr.en);

  for(EdgeIteratorRight ei=mr.edge.begin(); ei!=mr.edge.end(); ++ei)
    if(!(*ei).IsD() && (!selected || (*ei).IsS())){
      size_t ind=Index(mr,*ei);
      remap.edge[ind]=int(Index(ml,*ep));
      ++ep;
    }

  // face
  remap.face.resize(mr.face.size(), Remap::InvalidIndex());
  FaceIteratorLeft fp;
  size_t sfn = UpdateSelection<ConstMeshRight>::FaceCount(mr);
  if(selected) fp=Allocator<MeshLeft>::AddFaces(ml,sfn);
          else fp=Allocator<MeshLeft>::AddFaces(ml,mr.fn);

  for(FaceIteratorRight fi=mr.face.begin(); fi!=mr.face.end(); ++fi)
    if(!(*fi).IsD() && (!selected || (*fi).IsS())){
      size_t ind=Index(mr,*fi);
      remap.face[ind]=int(Index(ml,*fp));
      ++fp;
    }

  // hedge
  remap.hedge.resize(mr.hedge.size(),Remap::InvalidIndex());
  for(HEdgeIteratorRight hi=mr.hedge.begin(); hi!=mr.hedge.end(); ++hi)
    if(!(*hi).IsD() && (!selected || (*hi).IsS())){
      size_t ind=Index(mr,*hi);
      assert(remap.hedge[ind]==Remap::InvalidIndex());
      HEdgeIteratorLeft hp = Allocator<MeshLeft>::AddHEdges(ml,1);
      (*hp).ImportData(*(hi));
      remap.hedge[ind]=Index(ml,*hp);
    }
  
  remap.tetra.resize(mr.tetra.size(), Remap::InvalidIndex());
  for (TetraIteratorRight ti = mr.tetra.begin(); ti != mr.tetra.end(); ++ti)
    if (!(*ti).IsD() && (!selected || (*ti).IsS())) {
      size_t idx = Index(mr, *ti);
      assert (remap.tetra[idx] == Remap::InvalidIndex());
      TetraIteratorLeft tp = Allocator<MeshLeft>::AddTetras(ml, 1);
      (*tp).ImportData(*ti);
      remap.tetra[idx] = Index(ml, *tp);
    }

  // phase 2.
  // copy data from mr to its corresponding elements in ml and adjacencies

  // vertex
  for(VertexIteratorRight vi=mr.vert.begin();vi!=mr.vert.end();++vi)
    if( !(*vi).IsD() && (!selected || (*vi).IsS())){
      ml.vert[remap.vert[Index(mr,*vi)]].ImportData(*vi);
      if(adjFlag) ImportVertexAdj(ml,mr,ml.vert[remap.vert[Index(mr,*vi)]],*vi,remap);
    }

  // edge
  for(EdgeIteratorRight ei=mr.edge.begin();ei!=mr.edge.end();++ei)
    if(!(*ei).IsD() && (!selected || (*ei).IsS())){
      ml.edge[remap.edge[Index(mr,*ei)]].ImportData(*ei);
      // Edge to Vertex  Adj
      EdgeLeft &el = ml.edge[remap.edge[Index(mr,*ei)]];
      if(HasEVAdjacency(ml) && HasEVAdjacency(mr)){
        el.V(0) = &ml.vert[remap.vert[Index(mr,ei->cV(0))]];
        el.V(1) = &ml.vert[remap.vert[Index(mr,ei->cV(1))]];
      }
      if(adjFlag) ImportEdgeAdj(ml,mr,el,*ei,remap);
    }

  // face
  const size_t textureOffset =  ml.textures.size();
  bool WTFlag = HasPerWedgeTexCoord(mr) && (textureOffset>0);
  for(FaceIteratorRight fi=mr.face.begin();fi!=mr.face.end();++fi)
    if(!(*fi).IsD() && (!selected || (*fi).IsS()))
    {
      FaceLeft &fl = ml.face[remap.face[Index(mr,*fi)]];
      fl.Alloc(fi->VN());
      if(HasFVAdjacency(ml) && HasFVAdjacency(mr)){
        for(int i = 0; i < fl.VN(); ++i)
          fl.V(i) = &ml.vert[remap.vert[Index(mr,fi->cV(i))]];
      }
	  fl.ImportData(*fi);
      if(WTFlag)
        for(int i = 0; i < fl.VN(); ++i)
          fl.WT(i).n() += short(textureOffset);
      if(adjFlag)  ImportFaceAdj(ml,mr,ml.face[remap.face[Index(mr,*fi)]],*fi,remap);

    }

  // hedge
  for(HEdgeIteratorRight hi=mr.hedge.begin();hi!=mr.hedge.end();++hi)
    if(!(*hi).IsD() && (!selected || (*hi).IsS())){
      ml.hedge[remap.hedge[Index(mr,*hi)]].ImportData(*hi);
      ImportHEdgeAdj(ml,mr,ml.hedge[remap.hedge[Index(mr,*hi)]],*hi,remap,selected);
    }
  
  //tetra
  for(TetraIteratorRight ti = mr.tetra.begin(); ti != mr.tetra.end(); ++ti)
    if(!(*ti).IsD() && (!selected || (*ti).IsS()))
    {
      TetraLeft &tl = ml.tetra[remap.tetra[Index(mr,*ti)]];
      
      if(HasFVAdjacency(ml) && HasFVAdjacency(mr)){
        for(int i = 0; i < 4; ++i)
          tl.V(i) = &ml.vert[remap.vert[Index(mr,ti->cV(i))]];
      }
	    tl.ImportData(*ti);
      if(adjFlag)  ImportTetraAdj(ml, mr, ml.tetra[remap.tetra[Index(mr,*ti)]], *ti, remap);

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
		typename std::set< PointerToAttribute >::iterator al, ar;

		// per vertex attributes
		for(al = ml.vert_attr.begin(); al != ml.vert_attr.end(); ++al)
			if(!(*al)._name.empty()){
				ar =    mr.vert_attr.find(*al);
				if(ar!= mr.vert_attr.end()){
					id_r = 0;
					for(VertexIteratorRight vi=mr.vert.begin();vi!=mr.vert.end();++vi,++id_r)
						if( !(*vi).IsD() && (!selected || (*vi).IsS()))
							(*al)._handle->CopyValue(remap.vert[Index(mr,*vi)], id_r, (*ar)._handle);
				}
			}

		// per edge attributes
		for(al = ml.edge_attr.begin(); al != ml.edge_attr.end(); ++al)
			if(!(*al)._name.empty()){
				ar =    mr.edge_attr.find(*al);
				if(ar!= mr.edge_attr.end()){
					id_r = 0;
					for(EdgeIteratorRight ei=mr.edge.begin();ei!=mr.edge.end();++ei,++id_r)
						if( !(*ei).IsD() && (!selected || (*ei).IsS()))
							(*al)._handle->CopyValue(remap.edge[Index(mr,*ei)], id_r, (*ar)._handle);
				}
			}

		// per face attributes
		for(al = ml.face_attr.begin(); al != ml.face_attr.end(); ++al)
			if(!(*al)._name.empty()){
				ar =    mr.face_attr.find(*al);
				if(ar!= mr.face_attr.end()){
					id_r = 0;
					for(FaceIteratorRight fi=mr.face.begin();fi!=mr.face.end();++fi,++id_r)
						if( !(*fi).IsD() && (!selected || (*fi).IsS()))
							(*al)._handle->CopyValue(remap.face[Index(mr,*fi)], id_r, (*ar)._handle);
				}
			}

		// per tetra attributes
		for(al = ml.tetra_attr.begin(); al != ml.tetra_attr.end(); ++al)
			if(!(*al)._name.empty()){
				ar =    mr.tetra_attr.find(*al);
				if(ar!= mr.tetra_attr.end()){
					id_r = 0;
					for(TetraIteratorRight ti = mr.tetra.begin(); ti != mr.tetra.end(); ++ti, ++id_r)
						if( !(*ti).IsD() && (!selected || (*ti).IsS()))
							(*al)._handle->CopyValue(remap.tetra[Index(mr, *ti)], id_r, (*ar)._handle);
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

/**
 * @brief MeshAppendConst
 * @param ml
 * @param mr
 *
 * This is function is similar with the Mesh function, but does not
 * never update selections. In some cases, after the append,
 * selection of vertices may be inconsistent with face selection,
 * as explained above.
 * To avoid this, before using this function, call the following functions:
 *
 * \code{.cpp}
 * vcg::tri::UpdateSelection<MyMesh>::VertexFromEdgeLoose(mr,true);
 * vcg::tri::UpdateSelection<MyMesh>::VertexFromFaceLoose(mr,true);
 * \endcode
 *
 * or, use the Mesh function that takes a non-const Right Mesh argument.
 */
static void MeshAppendConst(
        MeshLeft& ml,
        const ConstMeshRight& mr,
        const bool selected = false,
        const bool adjFlag = false)
{
  // phase 1. allocate on ml vert,edge,face, hedge to accomodat those of mr
  // and build the remapping for all

  Remap remap;

  // vertex
  remap.vert.resize(mr.vert.size(), Remap::InvalidIndex());
  VertexIteratorLeft vp;
  size_t svn = UpdateSelection<ConstMeshRight>::VertexCount(mr);
  if(selected)
      vp=Allocator<MeshLeft>::AddVertices(ml,int(svn));
  else
      vp=Allocator<MeshLeft>::AddVertices(ml,mr.vn);

  ForEachVertex(mr, [&](const VertexRight& v)
  {
    if(!selected || v.IsS())
    {
      size_t ind=Index(mr,v);
      remap.vert[ind]=int(Index(ml,*vp));
      ++vp;
    }
  });
  // edge
  remap.edge.resize(mr.edge.size(), Remap::InvalidIndex());
  EdgeIteratorLeft ep;
  size_t sen = UpdateSelection<ConstMeshRight>::EdgeCount(mr);
  if(selected) ep=Allocator<MeshLeft>::AddEdges(ml,sen);
          else ep=Allocator<MeshLeft>::AddEdges(ml,mr.en);

  ForEachEdge(mr, [&](const EdgeRight& e)
  {
    if(!selected || e.IsS()){
      size_t ind=Index(mr,e);
      remap.edge[ind]=int(Index(ml,*ep));
      ++ep;
    }
  });

  // face
  remap.face.resize(mr.face.size(), Remap::InvalidIndex());
  FaceIteratorLeft fp;
  size_t sfn = UpdateSelection<ConstMeshRight>::FaceCount(mr);
  if(selected) fp=Allocator<MeshLeft>::AddFaces(ml,sfn);
          else fp=Allocator<MeshLeft>::AddFaces(ml,mr.fn);

  ForEachFace(mr, [&](const FaceRight& f)
  {
    if(!selected || f.IsS()){
      size_t ind=Index(mr,f);
      remap.face[ind]=int(Index(ml,*fp));
      ++fp;
    }
  });

  // hedge
  remap.hedge.resize(mr.hedge.size(),Remap::InvalidIndex());

  ForEachHEdge(mr, [&](const HEdgeRight& he)
  {
    if(!selected || he.IsS()){
      size_t ind=Index(mr,he);
      assert(remap.hedge[ind]==Remap::InvalidIndex());
      HEdgeIteratorLeft hp = Allocator<MeshLeft>::AddHEdges(ml,1);
      (*hp).ImportData(he);
      remap.hedge[ind]=Index(ml,*hp);
    }
  });

  remap.tetra.resize(mr.tetra.size(), Remap::InvalidIndex());

  ForEachTetra(mr, [&](const TetraRight& t)
  {
    if (!selected || t.IsS()) {
      size_t idx = Index(mr, t);
      assert (remap.tetra[idx] == Remap::InvalidIndex());
      TetraIteratorLeft tp = Allocator<MeshLeft>::AddTetras(ml, 1);
      (*tp).ImportData(t);
      remap.tetra[idx] = Index(ml, *tp);
    }
  });

  // phase 2.
  // copy data from mr to its corresponding elements in ml and adjacencies

  // vertex
  ForEachVertex(mr, [&](const VertexRight& v)
  {
    if(!selected || v.IsS()){
      ml.vert[remap.vert[Index(mr,v)]].ImportData(v);
      if(adjFlag) ImportVertexAdj(ml,mr,ml.vert[remap.vert[Index(mr,v)]],v,remap);
    }
  });

  // edge
  ForEachEdge(mr, [&](const EdgeRight& e)
  {
    if(!selected || e.IsS()){
      ml.edge[remap.edge[Index(mr,e)]].ImportData(e);
      // Edge to Vertex  Adj
      EdgeLeft &el = ml.edge[remap.edge[Index(mr,e)]];
      if(HasEVAdjacency(ml) && HasEVAdjacency(mr)){
        el.V(0) = &ml.vert[remap.vert[Index(mr,e.cV(0))]];
        el.V(1) = &ml.vert[remap.vert[Index(mr,e.cV(1))]];
      }
      if(adjFlag) ImportEdgeAdj(ml,mr,el,e,remap);
    }
  });

  // face
  const size_t textureOffset =  ml.textures.size();
  bool WTFlag = HasPerWedgeTexCoord(mr) && (textureOffset>0);
  ForEachFace(mr, [&](const FaceRight& f)
  {
    if(!selected || f.IsS())
    {
      FaceLeft &fl = ml.face[remap.face[Index(mr,f)]];
      fl.Alloc(f.VN());
      if(HasFVAdjacency(ml) && HasFVAdjacency(mr)){
        for(int i = 0; i < fl.VN(); ++i)
          fl.V(i) = &ml.vert[remap.vert[Index(mr,f.cV(i))]];
      }
      fl.ImportData(f);
      if(WTFlag)
        for(int i = 0; i < fl.VN(); ++i)
          fl.WT(i).n() += short(textureOffset);
      if(adjFlag)  ImportFaceAdj(ml,mr,ml.face[remap.face[Index(mr,f)]],f,remap);

    }
  });

  // hedge
  ForEachHEdge(mr, [&](const HEdgeRight& he)
  {
    if(!selected || he.IsS()){
      ml.hedge[remap.hedge[Index(mr,he)]].ImportData(he);
      ImportHEdgeAdj(ml,mr,ml.hedge[remap.hedge[Index(mr,he)]],he,remap,selected);
    }
  });

  //tetra
  ForEachTetra(mr, [&](const TetraRight& t)
  {
    if(!selected || t.IsS())
    {
      TetraLeft &tl = ml.tetra[remap.tetra[Index(mr,t)]];

      if(HasFVAdjacency(ml) && HasFVAdjacency(mr)){
        for(int i = 0; i < 4; ++i)
          tl.V(i) = &ml.vert[remap.vert[Index(mr,t.cV(i))]];
      }
        tl.ImportData(t);
      if(adjFlag)  ImportTetraAdj(ml, mr, ml.tetra[remap.tetra[Index(mr,t)]], t, remap);

    }
  });

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
	typename std::set< PointerToAttribute >::iterator al, ar;

	// per vertex attributes
	for (al = ml.vert_attr.begin(); al != ml.vert_attr.end(); ++al)
		if(!(*al)._name.empty())
		{
			ar = mr.vert_attr.find(*al);
			if (ar != mr.vert_attr.end())
			{
				id_r = 0;
				for (const auto & v : mr.vert)
				{
					if( !v.IsD() && (!selected || v.IsS()))
						(*al)._handle->CopyValue(remap.vert[Index(mr,v)], id_r, (*ar)._handle);
					++id_r;
				}
			}
		}

	// per edge attributes
	for (al = ml.edge_attr.begin(); al != ml.edge_attr.end(); ++al)
		if (!(*al)._name.empty())
		{
			ar = mr.edge_attr.find(*al);
			if (ar!= mr.edge_attr.end())
			{
				id_r = 0;
				for (const auto & e : mr.edge)
				{
					if( !e.IsD() && (!selected || e.IsS()))
						(*al)._handle->CopyValue(remap.edge[Index(mr,e)], id_r, (*ar)._handle);
					++id_r;
				}
			}
		}

	// per face attributes
	for (al = ml.face_attr.begin(); al != ml.face_attr.end(); ++al)
		if (!(*al)._name.empty())
		{
			ar = mr.face_attr.find(*al);
			if (ar!= mr.face_attr.end())
			{
				id_r = 0;
				for (const auto & f : mr.face)
				{
					if( !f.IsD() && (!selected || f.IsS()))
						(*al)._handle->CopyValue(remap.face[Index(mr,f)], id_r, (*ar)._handle);
					++id_r;
				}
			}
		}

	// per tetra attributes
	for (al = ml.tetra_attr.begin(); al != ml.tetra_attr.end(); ++al)
		if (!(*al)._name.empty())
		{
			ar = mr.tetra_attr.find(*al);
			if (ar!= mr.tetra_attr.end())
			{
				id_r = 0;
				for (const auto & t: mr.tetra)
				{
					if( !t.IsD() && (!selected || t.IsS()))
						(*al)._handle->CopyValue(remap.tetra[Index(mr, t)], id_r, (*ar)._handle);
					++id_r;
				}
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

/*! \brief Copy the second mesh over the first one.
  The first mesh is destroyed. If requested only the selected elements are copied.
*/
static void MeshCopy(MeshLeft& ml, ConstMeshRight& mr, bool selected=false, const bool adjFlag = false)
{
  ml.Clear();
  Mesh(ml,mr,selected,adjFlag);
  ml.bbox.Import(mr.bbox);
}

static void MeshCopyConst(MeshLeft& ml, const ConstMeshRight& mr, bool selected=false, const bool adjFlag = false)
{
  ml.Clear();
  MeshAppendConst(ml,mr,selected,adjFlag);
  ml.bbox.Import(mr.bbox);
}
/*! \brief %Append only the selected elements of second mesh to the first one.

  It is just a wrap of the main Append::Mesh()
  */
static void Selected(MeshLeft& ml, ConstMeshRight& mr)
{
  Mesh(ml,mr,true);
}

}; // end of class Append





} // End Namespace tri
} // End Namespace vcg


#endif


