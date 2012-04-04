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

#if defined(_MSC_VER)
#pragma warning( disable : 4804 )
#endif
#include <assert.h>
#include <string>
#include <vector>
#include <set>
#include <vcg/space/box3.h>
#include <vcg/space/color4.h>
#include <vcg/math/shot.h>

#include <vcg/container/simple_temporary_data.h>
#include <vcg/simplex/vertex/base.h>
#include <vcg/simplex/edge/base.h>
#include <vcg/simplex/face/base.h>
#include <vcg/connectors/hedge.h>
#include <vcg/complex/used_types.h>
#include <vcg/container/derivation_chain.h>
#include <vcg/complex/allocate.h>


#ifndef __VCG_MESH
#define __VCG_MESH

namespace vcg {
namespace tri {
/** \addtogroup trimesh */
/*@{*/
/*@{*/
/** Class Mesh.
    This is class for definition of a mesh.
		@param VertContainerType (Template Parameter) Specifies the type of the vertices container any the vertex type.
		@param FaceContainer (Template Parameter) Specifies the type of the faces container any the face type.
 */


 /* MeshTypeHolder is a class which is used to define the types in the mesh
*/

		template <class TYPESPOOL>
		struct BaseMeshTypeHolder{

				typedef bool ScalarType;
				typedef std::vector< typename TYPESPOOL::VertexType  >	CONTV;
				typedef std::vector< typename TYPESPOOL::EdgeType  >		CONTE;
				typedef std::vector< typename TYPESPOOL::FaceType >		CONTF;
				typedef std::vector< typename TYPESPOOL::HEdgeType  >		CONTH;

				typedef CONTV									VertContainer;
				typedef _Vertex 								VertexType;
				typedef typename TYPESPOOL::VertexPointer		VertexPointer;
				typedef const typename TYPESPOOL::VertexPointer ConstVertexPointer;
				typedef bool									CoordType;
				typedef typename CONTV::iterator				VertexIterator;
				typedef typename CONTV::const_iterator	ConstVertexIterator;

				typedef CONTE										EdgeContainer;
				typedef typename CONTE::value_type					EdgeType;
				typedef typename  TYPESPOOL::EdgePointer	EdgePointer;
				typedef typename CONTE::iterator				EdgeIterator;
				typedef typename CONTE::const_iterator	ConstEdgeIterator;

				typedef CONTF														FaceContainer;
				typedef typename CONTF::value_type					FaceType;
				typedef typename CONTF::const_iterator				ConstFaceIterator;
				typedef typename CONTF::iterator					FaceIterator;
				typedef typename TYPESPOOL::FacePointer	FacePointer;
				typedef const typename TYPESPOOL::FacePointer		ConstFacePointer;

				typedef CONTH														HEdgeContainer;
				typedef typename CONTH::value_type			HEdgeType;
				typedef typename TYPESPOOL::HEdgePointer					HEdgePointer;
				typedef typename CONTH::iterator				HEdgeIterator;
				typedef typename CONTH::const_iterator	ConstHEdgeIterator;


		};



		template <class T, typename CONT, class TRAIT >
						struct MeshTypeHolder: public T {};

		template <class T, typename CONT>
						struct MeshTypeHolder<T, CONT, AllTypes::AVertexType>: public T {
								typedef CONT VertContainer;
								typedef typename VertContainer::value_type VertexType;
								typedef VertexType * VertexPointer;
								typedef const VertexType * ConstVertexPointer;
								typedef typename VertexType::ScalarType ScalarType;
								typedef typename VertexType::CoordType CoordType;
								typedef typename VertContainer::iterator VertexIterator;
								typedef typename VertContainer::const_iterator ConstVertexIterator;
		};


	template <typename T, class CONT>
					struct MeshTypeHolder< T, CONT, AllTypes::AEdgeType>: public T{
								typedef CONT EdgeContainer;
								typedef typename EdgeContainer::value_type EdgeType;
								typedef typename EdgeContainer::value_type *  EdgePointer;
								typedef typename EdgeContainer::iterator EdgeIterator;
								typedef typename EdgeContainer::const_iterator ConstEdgeIterator;
};

	template <typename T, class CONT>
					struct MeshTypeHolder< T, CONT,  AllTypes::AFaceType>:public T {
								typedef CONT FaceContainer;
								typedef typename FaceContainer::value_type FaceType;
								typedef typename FaceContainer::const_iterator ConstFaceIterator;
								typedef typename FaceContainer::iterator FaceIterator;
								typedef FaceType * FacePointer;
								typedef const FaceType * ConstFacePointer;
				};

	template <typename T, class CONT>
					struct MeshTypeHolder< T, CONT, AllTypes::AHEdgeType>: public T{
								typedef CONT HEdgeContainer;
								typedef typename HEdgeContainer::value_type			HEdgeType;
								typedef typename HEdgeContainer::value_type *		HEdgePointer;
								typedef typename HEdgeContainer::iterator				HEdgeIterator;
								typedef typename HEdgeContainer::const_iterator ConstHEdgeIterator;
};


/*struct DummyContainer {};
template <class CONT> struct Deriver: public MeshTypeHolder<CONT, typename CONT::value_type::IAm>{};
template <> struct Deriver<DummyContainer>{}*/;

template <typename T, typename CONT> struct Der: public MeshTypeHolder<T,CONT, typename CONT::value_type::IAm>{};
struct DummyContainer{struct value_type{ typedef int IAm;}; };

template < class Container0 = DummyContainer, class Container1 = DummyContainer, class Container2 = DummyContainer, class Container3 = DummyContainer >
class TriMesh
	: public  MArity4<   BaseMeshTypeHolder<typename Container0::value_type::TypesPool>, Container0, Der ,Container1, Der, Container2, Der, Container3, Der>{
	public:

		typedef typename TriMesh::ScalarType		ScalarType;
		typedef typename TriMesh::VertContainer VertContainer;
		typedef typename TriMesh::EdgeContainer EdgeContainer;
		typedef typename TriMesh::FaceContainer FaceContainer;

		// types for vertex
		typedef typename TriMesh::VertexType						VertexType;
		typedef typename TriMesh::VertexPointer					VertexPointer;
		typedef typename TriMesh::ConstVertexPointer		ConstVertexPointer;
		typedef typename TriMesh::CoordType							CoordType;
		typedef typename TriMesh::VertexIterator				VertexIterator;
		typedef typename TriMesh::ConstVertexIterator		ConstVertexIterator;

		// types for edge
		typedef typename TriMesh::EdgeType							EdgeType;
		typedef typename TriMesh::EdgePointer						EdgePointer;
		typedef typename TriMesh::EdgeIterator					EdgeIterator;
		typedef typename TriMesh::ConstEdgeIterator			ConstEdgeIterator;

		//types for face
		typedef typename TriMesh::FaceType							FaceType;
		typedef typename TriMesh::ConstFaceIterator			ConstFaceIterator;
		typedef typename TriMesh::FaceIterator					FaceIterator;
		typedef typename TriMesh::FacePointer						FacePointer;
		typedef typename TriMesh::ConstFacePointer			ConstFacePointer;

		// types for hedge
		typedef typename TriMesh::HEdgeType							HEdgeType;
		typedef typename TriMesh::HEdgePointer					HEdgePointer;
		typedef typename TriMesh::HEdgeIterator					HEdgeIterator;
		typedef typename TriMesh::HEdgeContainer				HEdgeContainer;
		typedef typename TriMesh::ConstHEdgeIterator		ConstHEdgeIterator;

		typedef vcg::PointerToAttribute PointerToAttribute;

	typedef TriMesh<Container0, Container1,Container2,Container3> MeshType;

	typedef Box3<ScalarType> BoxType;

	/// Set of vertices 
	VertContainer vert;
	/// Actual number of vertices
	int vn;
	/// Set of faces
	FaceContainer face;
	/// Actual number of faces
	int fn;
	/// Set of edges
	EdgeContainer edge;
	/// Actual number of edges
	int en;
	/// Set of hedges
	HEdgeContainer hedge;
	/// Actual number of hedges
	int hn;
	/// Bounding box of the mesh
	Box3<ScalarType> bbox;
	
  /// Nomi di textures
	//
  std::vector<std::string> textures;
	//
  std::vector<std::string> normalmaps;

	int attrn;	// total numer of attribute created

	
	std::set< PointerToAttribute > vert_attr;
 	std::set< PointerToAttribute > edge_attr;
 	std::set< PointerToAttribute > face_attr;
	std::set< PointerToAttribute > mesh_attr;



	template <class ATTR_TYPE, class CONT>
	class AttributeHandle{
	public:
		AttributeHandle(){_handle=(SimpleTempData<CONT,ATTR_TYPE> *)NULL;}
		AttributeHandle( void *ah,const int & n):_handle ( (SimpleTempData<CONT,ATTR_TYPE> *)ah ),n_attr(n){}
                AttributeHandle operator = ( const PointerToAttribute & pva){
			_handle = (SimpleTempData<CONT,ATTR_TYPE> *)pva._handle;
			n_attr = pva.n_attr;
			return (*this);
		}

		//pointer to the SimpleTempData that stores the attribute
		SimpleTempData<CONT,ATTR_TYPE> * _handle;

		// its attribute number
		int n_attr; 

		// access function
		template <class RefType>
		ATTR_TYPE & operator [](const RefType  & i){return (*_handle)[i];}
	};

	template <class ATTR_TYPE>
	class PerVertexAttributeHandle: public AttributeHandle<ATTR_TYPE,VertContainer>{
	public:
		PerVertexAttributeHandle():AttributeHandle<ATTR_TYPE,VertContainer>(){}
                PerVertexAttributeHandle( void *ah,const int & n):AttributeHandle<ATTR_TYPE,VertContainer>(ah,n){}
	};

	template <class ATTR_TYPE>
	class PerFaceAttributeHandle: public AttributeHandle<ATTR_TYPE,FaceContainer>{
	public:
		PerFaceAttributeHandle():AttributeHandle<ATTR_TYPE,FaceContainer>(){}
                PerFaceAttributeHandle( void *ah,const int & n):AttributeHandle<ATTR_TYPE,FaceContainer>(ah,n){}
	};

	template <class ATTR_TYPE>
	class PerEdgeAttributeHandle:  public AttributeHandle<ATTR_TYPE,EdgeContainer>{
	public: 
		PerEdgeAttributeHandle():AttributeHandle<ATTR_TYPE,EdgeContainer>(){}
                PerEdgeAttributeHandle( void *ah,const int & n):AttributeHandle<ATTR_TYPE,EdgeContainer>(ah,n){}
	};

	template <class ATTR_TYPE>
	class PerMeshAttributeHandle{
	public:
		PerMeshAttributeHandle(){_handle=NULL;}
		PerMeshAttributeHandle(void *ah,const int & n):_handle ( (Attribute<ATTR_TYPE> *)ah ),n_attr(n){}
		PerMeshAttributeHandle operator = ( const PerMeshAttributeHandle & pva){ 
			_handle = (Attribute<ATTR_TYPE> *)pva._handle;
			n_attr = pva.n_attr;
			return (*this);
		}

		Attribute<ATTR_TYPE> * _handle;
		int n_attr;
		ATTR_TYPE & operator ()(){ return *((Attribute<ATTR_TYPE> *)_handle)->attribute;}
	};


	// the camera member (that should keep the intrinsics) is no more needed since 2006, when intrisncs moved into the Shot structure
	//Camera<ScalarType> camera; // intrinsic
	Shot<ScalarType> shot;		// intrinsic && extrinsic

private:
	/// The per-mesh color. Not very useful and meaningful...
	Color4b c;
public:

	inline const Color4b & C() const
	{
		return c;
	}

	inline Color4b & C()
	{
		return c;
	}


	/// Default constructor
	TriMesh()
	{
    Clear();
	}

	/// destructor
	~TriMesh()
	{
		typename std::set< PointerToAttribute>::iterator i;
		for( i = vert_attr.begin(); i != vert_attr.end(); ++i) 
			delete ((SimpleTempDataBase*)(*i)._handle);
		for( i = edge_attr.begin(); i != edge_attr.end(); ++i)
			delete ((SimpleTempDataBase*)(*i)._handle);
		for( i = face_attr.begin(); i != face_attr.end(); ++i)
			delete ((SimpleTempDataBase*)(*i)._handle);
		for( i = mesh_attr.begin(); i != mesh_attr.end(); ++i) 
			delete ((SimpleTempDataBase*)(*i)._handle);

		FaceIterator fi;
		for(fi = face.begin(); fi != face.end(); ++fi) (*fi).Dealloc();
	}

	 int Mem(const int & nv, const int & nf) const  {
		typename std::set< PointerToAttribute>::const_iterator i;
		int size = 0;
		size += sizeof(TriMesh)+sizeof(VertexType)*nv+sizeof(FaceType)*nf;

		for( i = vert_attr.begin(); i != vert_attr.end(); ++i) 
			size += ((SimpleTempDataBase*)(*i)._handle)->SizeOf()*nv;
		for( i = edge_attr.begin(); i != edge_attr.end(); ++i)
			size += ((SimpleTempDataBase*)(*i)._handle)->SizeOf()*en;
		for( i = face_attr.begin(); i != face_attr.end(); ++i)
			size +=  ((SimpleTempDataBase*)(*i)._handle)->SizeOf()*nf;
		for( i = mesh_attr.begin(); i != mesh_attr.end(); ++i) 
			size +=  ((SimpleTempDataBase*)(*i)._handle)->SizeOf();

		return size;
	}
	int MemUsed() const  {return Mem(vert.size(),face.size());}
	inline int MemNeeded() const {return Mem(vn,fn);}



/// Function to destroy the mesh
void Clear()
{
	vert.clear();
	face.clear();
  edge.clear();
//	textures.clear();
//	normalmaps.clear();
	vn = 0;
	en = 0;
	fn = 0;
	hn = 0;
  imark = 0;
  attrn = 0;
  C()=Color4b::Gray;
}


int & SimplexNumber(){ return fn;}
int & VertexNumber(){ return vn;}

/// The incremental mark
int imark;

/// Calcolo del volume di una mesh chiusa
ScalarType Volume()
{
  FaceIterator fi;
  int j,k;
  ScalarType V = 0;
  CoordType T,N,B;
 
  for(fi = face.begin(); fi!=face.end(); ++fi)
  {
	  for(j = 0; j < 3; ++j)
	  {
	    /*calcolo tangente, normale e binormale (6 volte)*/
	    k = (j+1)%3;
	    T = (*fi).P(k) - (*fi).P(j);
	    T.Normalize();
	    B = ( (*fi).P( k     ) - (*fi).P(j) ) ^
	        ( (*fi).P((k+1)%3) - (*fi).P(j) ) ;
	    B.Normalize();
			N = T ^ B;
     
	    CoordType pj = (*fi).P(j);
	    CoordType pk = (*fi).P(k);
   

	    V +=  (pk*  T )*(pk*N)*(pk*B);
	    V +=  (pj*(-T))*(pj*N)*(pj*B);
     }
  }
	return V/6.0;
}

private:
	// TriMesh cannot be copied. Use Append (see vcg/complex/append.h)
  TriMesh operator =(const TriMesh &  /*m*/){assert(0);return TriMesh();}
	TriMesh(const TriMesh & ){}

};	// end class Mesh

/// Initialize the imark-system of the faces
template <class MeshType> inline  void InitFaceIMark(MeshType & m)
{
	typename MeshType::FaceIterator f;
	
	for(f=m.face.begin();f!=m.face.end();++f)
		if( !(*f).IsD() && (*f).IsR() && (*f).IsW() )
			(*f).InitIMark();
}

/// Initialize the imark-system of the vertices
template <class MeshType> inline  void InitVertexIMark(MeshType & m)
{
	typename MeshType::VertexIterator vi;

	for(vi=m.vert.begin();vi!=m.vert.end();++vi)
		if( !(*vi).IsD() && (*vi).IsRW() )
			(*vi).InitIMark();
}
/** Access function to the incremental mark. 
*/
template <class MeshType> inline int & IMark(MeshType & m){return m.imark;}

/** Check if the vertex incremental mark matches the one of the mesh. 
	@param v Vertex pointer
*/
template <class MeshType> inline bool IsMarked(MeshType & m, typename MeshType::ConstVertexPointer  v )  { return v->IMark() == m.imark; }
/** Check if the face incremental mark matches the one of the mesh. 
	@param v Face pointer
*/
template <class MeshType> inline bool IsMarked( MeshType & m,typename MeshType::ConstFacePointer f )  { return f->IMark() == m.imark; }
/** Set the vertex incremental mark of the vertex to the one of the mesh.
	@param v Vertex pointer
*/
template <class MeshType> inline void Mark(MeshType & m, typename MeshType::VertexPointer v )  { v->IMark() = m.imark; }
/** Set the face incremental mark of the vertex to the one of the mesh.
	@param v Vertex pointer
*/
template <class MeshType> inline void Mark(MeshType & m, typename MeshType::FacePointer f )  { f->IMark() = m.imark; }
/// Unmark the mesh
template <class MeshType> inline void UnMarkAll(MeshType & m) { ++m.imark; }


template < class CType0, class CType1 , class CType2, class CType3>
bool HasPerVertexVEAdjacency (const TriMesh < CType0, CType1, CType2, CType3> & /*m*/) {return TriMesh < CType0 , CType1, CType2, CType3>::VertContainer::value_type::HasVEAdjacency();}
template < class  CType0, class CType1, class CType2 , class CType3>
bool HasPerEdgeVEAdjacency   (const TriMesh < CType0, CType1, CType2, CType3> & /*m*/) {return TriMesh < CType0 , CType1, CType2, CType3>::EdgeContainer::value_type::HasVEAdjacency();}

template < class VertexType> bool VertexVectorHasPerFaceVFAdjacency     (const std::vector<VertexType> &) {  return VertexType::HasVFAdjacency(); }
template < class FaceType  > bool   FaceVectorHasPerFaceVFAdjacency     (const std::vector<FaceType  > &) {  return FaceType::HasVFAdjacency(); }

template < class TriMeshType> bool   HasPerFaceVFAdjacency     (const TriMeshType &m) { return tri::FaceVectorHasPerFaceVFAdjacency  (m.vert); }
template < class TriMeshType> bool HasPerVertexVFAdjacency     (const TriMeshType &m) { return tri::VertexVectorHasPerFaceVFAdjacency(m.vert); }


template < class VertexType> bool VertexVectorHasPerVertexQuality     (const std::vector<VertexType> &) {  return VertexType::HasQuality     (); }
template < class VertexType> bool VertexVectorHasPerVertexNormal      (const std::vector<VertexType> &) {  return VertexType::HasNormal      (); }
template < class VertexType> bool VertexVectorHasPerVertexColor       (const std::vector<VertexType> &) {  return VertexType::HasColor       (); }
template < class VertexType> bool VertexVectorHasPerVertexMark        (const std::vector<VertexType> &) {  return VertexType::HasMark        (); }
template < class VertexType> bool VertexVectorHasPerVertexFlags       (const std::vector<VertexType> &) {  return VertexType::HasFlags       (); }
template < class VertexType> bool VertexVectorHasPerVertexRadius      (const std::vector<VertexType> &) {  return VertexType::HasRadius      (); }
template < class VertexType> bool VertexVectorHasPerVertexCurvature   (const std::vector<VertexType> &) {  return VertexType::HasCurvature   (); }
template < class VertexType> bool VertexVectorHasPerVertexCurvatureDir(const std::vector<VertexType> &) {  return VertexType::HasCurvatureDir(); }
template < class VertexType> bool VertexVectorHasPerVertexTexCoord    (const std::vector<VertexType> &) {  return VertexType::HasTexCoord    (); }

template < class TriMeshType> bool HasPerVertexQuality     (const TriMeshType &m) { return tri::VertexVectorHasPerVertexQuality     (m.vert); }
template < class TriMeshType> bool HasPerVertexNormal      (const TriMeshType &m) { return tri::VertexVectorHasPerVertexNormal      (m.vert); }
template < class TriMeshType> bool HasPerVertexColor       (const TriMeshType &m) { return tri::VertexVectorHasPerVertexColor       (m.vert); }
template < class TriMeshType> bool HasPerVertexMark        (const TriMeshType &m) { return tri::VertexVectorHasPerVertexMark        (m.vert); }
template < class TriMeshType> bool HasPerVertexFlags       (const TriMeshType &m) { return tri::VertexVectorHasPerVertexFlags       (m.vert); }
template < class TriMeshType> bool HasPerVertexRadius      (const TriMeshType &m) { return tri::VertexVectorHasPerVertexRadius      (m.vert); }
template < class TriMeshType> bool HasPerVertexCurvature   (const TriMeshType &m) { return tri::VertexVectorHasPerVertexCurvature   (m.vert); }
template < class TriMeshType> bool HasPerVertexCurvatureDir(const TriMeshType &m) { return tri::VertexVectorHasPerVertexCurvatureDir(m.vert); }
template < class TriMeshType> bool HasPerVertexTexCoord    (const TriMeshType &m) { return tri::VertexVectorHasPerVertexTexCoord    (m.vert); }


template < class FaceType>    bool FaceVectorHasPerWedgeColor   (const std::vector<FaceType> &) {  return FaceType::HasWedgeColor   (); }
template < class FaceType>    bool FaceVectorHasPerWedgeNormal  (const std::vector<FaceType> &) {  return FaceType::HasWedgeNormal  (); }
template < class FaceType>    bool FaceVectorHasPerWedgeTexCoord(const std::vector<FaceType> &) {  return FaceType::HasWedgeTexCoord(); }

template < class TriMeshType> bool HasPerWedgeColor   (const TriMeshType &m) { return tri::FaceVectorHasPerWedgeColor   (m.face); }
template < class TriMeshType> bool HasPerWedgeNormal  (const TriMeshType &m) { return tri::FaceVectorHasPerWedgeNormal  (m.face); }
template < class TriMeshType> bool HasPerWedgeTexCoord(const TriMeshType &m) { return tri::FaceVectorHasPerWedgeTexCoord(m.face); }

template < class  CType0, class CType1, class CType2 , class CType3>
bool HasPolyInfo (const TriMesh < CType0, CType1, CType2, CType3> & /*m*/) {return TriMesh < CType0 , CType1, CType2, CType3>::FaceContainer::value_type::HasPolyInfo();}

template < class FaceType>    bool FaceVectorHasPerFaceFlags  (const std::vector<FaceType> &) {  return FaceType::HasFlags      (); }
template < class FaceType>    bool FaceVectorHasPerFaceNormal (const std::vector<FaceType> &) {  return FaceType::HasFaceNormal (); }
template < class FaceType>    bool FaceVectorHasPerFaceColor  (const std::vector<FaceType> &) {  return FaceType::HasFaceColor  (); }
template < class FaceType>    bool FaceVectorHasPerFaceMark   (const std::vector<FaceType> &) {  return FaceType::HasFaceMark   (); }
template < class FaceType>    bool FaceVectorHasPerFaceQuality(const std::vector<FaceType> &) {  return FaceType::HasFaceQuality(); }
template < class FaceType>    bool FaceVectorHasFFAdjacency   (const std::vector<FaceType> &) {  return FaceType::HasFFAdjacency(); }
template < class FaceType>    bool FaceVectorHasFEAdjacency   (const std::vector<FaceType> &) {  return FaceType::HasFEAdjacency(); }
template < class FaceType>    bool FaceVectorHasFVAdjacency   (const std::vector<FaceType> &) {  return FaceType::HasFVAdjacency(); }
template < class FaceType>    bool FaceVectorHasPerFaceCurvatureDir   (const std::vector<FaceType> &) {  return FaceType::HasCurvatureDir(); }

template < class TriMeshType> bool HasPerFaceFlags       (const TriMeshType &m) { return tri::FaceVectorHasPerFaceFlags       (m.face); }
template < class TriMeshType> bool HasPerFaceNormal      (const TriMeshType &m) { return tri::FaceVectorHasPerFaceNormal      (m.face); }
template < class TriMeshType> bool HasPerFaceColor       (const TriMeshType &m) { return tri::FaceVectorHasPerFaceColor       (m.face); }
template < class TriMeshType> bool HasPerFaceMark        (const TriMeshType &m) { return tri::FaceVectorHasPerFaceMark        (m.face); }
template < class TriMeshType> bool HasPerFaceQuality     (const TriMeshType &m) { return tri::FaceVectorHasPerFaceQuality     (m.face); }
template < class TriMeshType> bool HasPerFaceCurvatureDir(const TriMeshType &m) { return tri::FaceVectorHasPerFaceCurvatureDir(m.face); }
template < class TriMeshType> bool HasFFAdjacency   (const TriMeshType &m) { return tri::FaceVectorHasFFAdjacency   (m.face); }
template < class TriMeshType> bool HasFEAdjacency   (const TriMeshType &m) { return tri::FaceVectorHasFEAdjacency   (m.face); }
template < class TriMeshType> bool HasFVAdjacency   (const TriMeshType &m) { return tri::FaceVectorHasFVAdjacency   (m.face); }



template < class  CType0, class CType1, class CType2 , class CType3>
bool HasVEAdjacency (const TriMesh < CType0, CType1, CType2, CType3> & /*m*/) {return TriMesh < CType0 , CType1, CType2, CType3>::VertContainer::value_type::HasVEAdjacency();}

template < class  CType0, class CType1, class CType2 , class CType3>
bool HasVHAdjacency (const TriMesh < CType0, CType1, CType2, CType3> & /*m*/) {return TriMesh < CType0 , CType1, CType2, CType3>::VertContainer::value_type::HasVHAdjacency();}

template < class  CType0, class CType1, class CType2 , class CType3>
bool HasEVAdjacency (const TriMesh < CType0, CType1, CType2, CType3> & /*m*/) {return TriMesh < CType0 , CType1, CType2, CType3>::EdgeType::HasEVAdjacency();}

template < class  CType0, class CType1, class CType2 , class CType3>
bool HasEEAdjacency (const TriMesh < CType0, CType1, CType2, CType3> & /*m*/) {return TriMesh < CType0 , CType1, CType2, CType3>::EdgeType::HasEEAdjacency();}

template < class  CType0, class CType1, class CType2 , class CType3>
bool HasEFAdjacency (const TriMesh < CType0, CType1, CType2, CType3> & /*m*/) {return TriMesh < CType0 , CType1, CType2, CType3>::EdgeType::HasEFAdjacency();}

template < class  CType0, class CType1, class CType2 , class CType3>
bool HasEHAdjacency (const TriMesh < CType0, CType1, CType2, CType3> & /*m*/) {return TriMesh < CType0 , CType1, CType2, CType3>::EdgeType::HasEHAdjacency();}

template < class  CType0, class CType1, class CType2 , class CType3>
bool HasFHAdjacency (const TriMesh < CType0, CType1, CType2, CType3> & /*m*/) {return TriMesh < CType0 , CType1, CType2, CType3>::FaceType::HasFHAdjacency();}

template < class  CType0, class CType1, class CType2 , class CType3>
bool HasHVAdjacency (const TriMesh < CType0, CType1, CType2, CType3> & /*m*/) {return TriMesh < CType0 , CType1, CType2, CType3>::HEdgeType::HasHVAdjacency();}

template < class  CType0, class CType1, class CType2 , class CType3>
bool HasHEAdjacency (const TriMesh < CType0, CType1, CType2, CType3> & /*m*/) {return TriMesh < CType0 , CType1, CType2, CType3>::HEdgeType::HasHEAdjacency();}

template < class  CType0, class CType1, class CType2 , class CType3>
bool HasHFAdjacency (const TriMesh < CType0, CType1, CType2, CType3> & /*m*/) {return TriMesh < CType0 , CType1, CType2, CType3>::HEdgeType::HasHFAdjacency();}

template < class  CType0, class CType1, class CType2 , class CType3>
bool HasHNextAdjacency (const TriMesh < CType0, CType1, CType2, CType3> & /*m*/) {return TriMesh< CType0,   CType1,   CType2 ,  CType3>::HEdgeType::HasHNextAdjacency();}

template < class  CType0, class CType1, class CType2 , class CType3>
bool HasHPrevAdjacency (const TriMesh < CType0, CType1, CType2, CType3> & /*m*/) {return TriMesh< CType0,   CType1,   CType2 , CType3>::HEdgeType::HasHPrevAdjacency();}

template < class  CType0, class CType1, class CType2 , class CType3>
bool HasHOppAdjacency (const TriMesh < CType0, CType1, CType2, CType3> & /*m*/)  {return TriMesh< CType0,   CType1,   CType2 ,  CType3>::HEdgeType::HasHOppAdjacency();}

//template < class CType0, class CType1 , class CType2, class CType3>
//bool HasVFAdjacency (const TriMesh < CType0 , CType1,   CType2, CType3> &  m ) {
//		// gcc 4.4: if the expressions assigned to a1 and a2 are replaced in the assert we get a compilation error
//		// for the macro assert
//		bool a1 =  TriMesh < CType0 , CType1,   CType2, CType3>::FaceContainer::value_type::HasVFAdjacency();
//		bool a2 =  TriMesh < CType0 , CType1,   CType2, CType3>::VertContainer::value_type::HasVFAdjacency();
//		// a1 and a2 are still evaluated but not referenced, this causes a warning
//		(void)a1;
//		(void)a2;
//		assert(a1==a2);
//
//		return   vcg::tri::HasPerVertexVFAdjacency<   CType0,   CType1 ,   CType2,   CType3>(m) &&
//						 vcg::tri::HasPerFaceVFAdjacency<   CType0,   CType1 ,   CType2,   CType3>(m) ;
//}

template <class MESH_TYPE>
bool HasPerVertexAttribute(const MESH_TYPE &m,   std::string   name){
		typename std::set< typename MESH_TYPE::PointerToAttribute>::const_iterator ai;
		typename MESH_TYPE::PointerToAttribute h; 
		h._name = name;
		ai = m.vert_attr.find(h);
		return (ai!= m.vert_attr.end() ) ;
}
template <class MESH_TYPE>
bool HasPerFaceAttribute(const MESH_TYPE &m,   std::string   name){
		typename std::set< typename MESH_TYPE::PointerToAttribute>::const_iterator ai;
		typename MESH_TYPE::PointerToAttribute h; 
		h._name = name;
		ai = m.face_attr.find(h);
		return (ai!= m.face_attr.end() ) ;
}

template <class MESH_TYPE>
bool HasPerMeshAttribute(const MESH_TYPE &m,   std::string   name){
		typename std::set< typename MESH_TYPE::PointerToAttribute>::const_iterator ai;
		typename MESH_TYPE::PointerToAttribute h; 
		h._name = name;
		ai = m.mesh_attr.find(h);
		return (ai!= m.mesh_attr.end() ) ;
}

/*@}*/
/*@}*/
}	 // end namespace
}	 // end namespace


#endif

