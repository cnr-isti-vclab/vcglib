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

#ifndef __VCG_COMPLEX_BASE
#define __VCG_COMPLEX_BASE

#include <typeindex>
#include <set>

#include <vcg/container/simple_temporary_data.h>

#include "exception.h"
#include "used_types.h"

namespace vcg {

class PointerToAttribute
{
public:
	SimpleTempDataBase * _handle;       // pointer to the SimpleTempData that stores the attribute
	std::string _name;					// name of the attribute
	int _sizeof;						// size of the attribute type (used only with VMI loading)
	int _padding;						// padding 	(used only with VMI loading)

	int n_attr;							// unique ID of the attribute
	std::type_index _type;
	void Resize(size_t sz){((SimpleTempDataBase *)_handle)->Resize(sz);}
	void Reorder(std::vector<size_t> & newVertIndex){((SimpleTempDataBase *)_handle)->Reorder(newVertIndex);}
	bool operator<(const  PointerToAttribute    b) const {	return(_name.empty()&&b._name.empty())?(_handle < b._handle):( _name < b._name);}

	PointerToAttribute(): _type(typeid(void)) { };
};


namespace tri {
/** \addtogroup trimesh */
/*@{*/


/* MeshTypeHolder is a class which is used to define the types in the mesh
*/

template <class TYPESPOOL>
struct BaseMeshTypeHolder{

	typedef bool ScalarType;
	typedef std::vector< typename TYPESPOOL::VertexType > CONTV;
	typedef std::vector< typename TYPESPOOL::EdgeType >   CONTE;
	typedef std::vector< typename TYPESPOOL::FaceType >   CONTF;
	typedef std::vector< typename TYPESPOOL::HEdgeType >  CONTH;
	typedef std::vector<typename TYPESPOOL::TetraType>    CONTT;


	typedef CONTV									VertContainer;
	typedef typename CONTV::value_type 								VertexType;
	typedef typename TYPESPOOL::VertexPointer		VertexPointer;
	typedef const typename TYPESPOOL::VertexPointer ConstVertexPointer;
	typedef bool									CoordType;
	typedef typename CONTV::iterator				VertexIterator;
	typedef typename CONTV::const_iterator	ConstVertexIterator;

	typedef CONTE										EdgeContainer;
	typedef typename CONTE::value_type					EdgeType;
	typedef typename  TYPESPOOL::EdgePointer	EdgePointer;
	typedef const typename  TYPESPOOL::EdgePointer	ConstEdgePointer;
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

	typedef CONTT TetraContainer;
	typedef typename CONTT::value_type TetraType;
	typedef typename TYPESPOOL::TetraPointer TetraPointer;
	typedef const typename TYPESPOOL::TetraPointer ConstTetraPointer;
	typedef typename CONTT::iterator TetraIterator;
	typedef typename CONTT::const_iterator ConstTetraIterator;


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
	typedef EdgeType *  EdgePointer;
	typedef const EdgeType * ConstEdgePointer;
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

template <typename T, class CONT>
struct MeshTypeHolder<T, CONT, AllTypes::ATetraType> : public T
{
	typedef CONT TetraContainer;
	typedef typename TetraContainer::value_type TetraType;
	typedef TetraType *TetraPointer;
	typedef const TetraType *ConstTetraPointer;
	typedef typename TetraContainer::iterator TetraIterator;
	typedef typename TetraContainer::const_iterator ConstTetraIterator;
};

template <typename T, typename CONT> struct Der: public MeshTypeHolder<T,CONT, typename CONT::value_type::IAm>{};
struct DummyContainer{struct value_type{ typedef int IAm;}; };
/** \brief The official \b mesh class

As explained in \ref basic_concepts, this class is templated over a list of container of simplexes (like vertex, face, edges)
 */

template < class Container0 = DummyContainer, class Container1 = DummyContainer, class Container2 = DummyContainer, class Container3 = DummyContainer, class Container4 = DummyContainer >
class TriMesh
        : public  MArity5<   BaseMeshTypeHolder<typename Container0::value_type::TypesPool>, Container0, Der ,Container1, Der, Container2, Der, Container3, Der, Container4, Der >{
public:

	typedef typename TriMesh::ScalarType		ScalarType;
	typedef typename TriMesh::VertContainer VertContainer;
	typedef typename TriMesh::EdgeContainer EdgeContainer;
	typedef typename TriMesh::FaceContainer FaceContainer;
	typedef typename TriMesh::TetraContainer TetraContainer;

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
	typedef typename TriMesh::ConstEdgePointer					ConstEdgePointer;
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

	// types for tetra
	typedef typename TriMesh::TetraType						TetraType;
	typedef typename TriMesh::TetraPointer					TetraPointer;
	typedef typename TriMesh::TetraIterator					TetraIterator;
	typedef typename TriMesh::ConstTetraIterator		    ConstTetraIterator;

	typedef vcg::PointerToAttribute PointerToAttribute;

	typedef TriMesh<Container0, Container1, Container2, Container3, Container4> MeshType;

	typedef Box3<ScalarType> BoxType;

	/// Container of vertices, usually a vector.
	VertContainer vert;
	/// Current number of vertices; this member is for internal use only. You should always use the VN() member
	int vn;
	/// Current number of vertices
	inline int VN() const { return vn; }

	/// Container of edges, usually a vector.
	EdgeContainer edge;
	/// Current number of edges; this member is for internal use only. You should always use the EN() member
	int en;
	/// Current number of edges
	inline int EN() const { return en; }

	/// Container of faces, usually a vector.
	FaceContainer face;
	/// Current number of faces; this member is for internal use only. You should always use the FN() member
	int fn;
	/// Current number of faces
	inline int FN() const { return fn; }

	/// Container of half edges, usually a vector.
	HEdgeContainer hedge;
	/// Current number of halfedges; this member is for internal use only. You should always use the HN() member
	int hn;
	/// Current number of halfedges;
	inline int HN() const { return hn; }

	/// Container of tetras, usually a vector.
	TetraContainer tetra;
	/// Current number of tetras; this member is for internal use only. You should always use the TN() member
	int tn;
	/// Current number of tetras;
	inline int TN() const { return tn; }

	/// Bounding box of the mesh
	Box3<typename TriMesh::VertexType::CoordType::ScalarType> bbox;

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
	std::set< PointerToAttribute > tetra_attr;


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
		void resize(size_t /*size*/) { };
	};

	template <class ATTR_TYPE, class CONT>
	class ConstAttributeHandle{
	public:
		ConstAttributeHandle(){_handle=(SimpleTempData<CONT,ATTR_TYPE> *)nullptr;}
		ConstAttributeHandle( const void *ah,const int & n):_handle ( (const SimpleTempData<CONT,ATTR_TYPE> *)ah ),n_attr(n){}


		//pointer to the SimpleTempData that stores the attribute
		const SimpleTempData<CONT,ATTR_TYPE> * _handle;

		// its attribute number
		int n_attr;

		// access function
		template <class RefType>
		const ATTR_TYPE & operator [](const RefType  & i) const {return (*_handle)[i];}
		void resize(size_t /*size*/) { };
	};

	template <class ATTR_TYPE>
	class PerVertexAttributeHandle: public AttributeHandle<ATTR_TYPE,VertContainer>{
	public:
		PerVertexAttributeHandle():AttributeHandle<ATTR_TYPE,VertContainer>(){}
		PerVertexAttributeHandle( void *ah,const int & n):AttributeHandle<ATTR_TYPE,VertContainer>(ah,n){}
	};

	template <class ATTR_TYPE>
	class ConstPerVertexAttributeHandle: public ConstAttributeHandle<ATTR_TYPE,VertContainer>{
	public:
		ConstPerVertexAttributeHandle():ConstAttributeHandle<ATTR_TYPE,VertContainer>(){}
		ConstPerVertexAttributeHandle( const void *ah,const int & n):ConstAttributeHandle<ATTR_TYPE,VertContainer>(ah,n){}
	};

	template <class ATTR_TYPE>
	class PerFaceAttributeHandle: public AttributeHandle<ATTR_TYPE,FaceContainer>{
	public:
		PerFaceAttributeHandle():AttributeHandle<ATTR_TYPE,FaceContainer>(){}
		PerFaceAttributeHandle( void *ah,const int & n):AttributeHandle<ATTR_TYPE,FaceContainer>(ah,n){}
	};

	template <class ATTR_TYPE>
	class ConstPerFaceAttributeHandle: public ConstAttributeHandle<ATTR_TYPE,FaceContainer>{
	public:
		ConstPerFaceAttributeHandle():ConstAttributeHandle<ATTR_TYPE,FaceContainer>(){}
		ConstPerFaceAttributeHandle( void *ah,const int & n):ConstAttributeHandle<ATTR_TYPE,FaceContainer>(ah,n){}
	};

	template <class ATTR_TYPE>
	class PerEdgeAttributeHandle:  public AttributeHandle<ATTR_TYPE,EdgeContainer>{
	public:
		PerEdgeAttributeHandle():AttributeHandle<ATTR_TYPE,EdgeContainer>(){}
		PerEdgeAttributeHandle( void *ah,const int & n):AttributeHandle<ATTR_TYPE,EdgeContainer>(ah,n){}
	};

	template <class ATTR_TYPE>
	class ConstPerEdgeAttributeHandle:  public ConstAttributeHandle<ATTR_TYPE,EdgeContainer>{
	public:
		ConstPerEdgeAttributeHandle():ConstAttributeHandle<ATTR_TYPE,EdgeContainer>(){}
		ConstPerEdgeAttributeHandle( void *ah,const int & n):ConstAttributeHandle<ATTR_TYPE,EdgeContainer>(ah,n){}
	};

	template <class ATTR_TYPE>
	class PerTetraAttributeHandle : public AttributeHandle<ATTR_TYPE, TetraContainer>
	{
	public:
		PerTetraAttributeHandle() : AttributeHandle<ATTR_TYPE, TetraContainer>() {}
		PerTetraAttributeHandle(void *ah, const int &n) : AttributeHandle<ATTR_TYPE, TetraContainer>(ah, n) {}
	};

	template <class ATTR_TYPE>
	class ConstPerTetraAttributeHandle : public ConstAttributeHandle<ATTR_TYPE, TetraContainer>
	{
	public:
		ConstPerTetraAttributeHandle() : ConstAttributeHandle<ATTR_TYPE, TetraContainer>() {}
		ConstPerTetraAttributeHandle(void *ah, const int &n) : ConstAttributeHandle<ATTR_TYPE, TetraContainer>(ah, n) {}
	};

	template <class ATTR_TYPE>
	class PerMeshAttributeHandle{
	public:
		PerMeshAttributeHandle(){_handle=NULL;}
		PerMeshAttributeHandle(void *ah,const int & n):_handle ( (Attribute<ATTR_TYPE> *)ah ),n_attr(n){}
		//PerMeshAttributeHandle operator = ( const PerMeshAttributeHandle & pva){
		//	_handle = (Attribute<ATTR_TYPE> *)pva._handle;
		//	n_attr = pva.n_attr;
		//	return (*this);
		//}

		Attribute<ATTR_TYPE> * _handle;
		int n_attr;
		ATTR_TYPE & operator ()(){ return *((ATTR_TYPE*) (_handle->DataBegin()));}
	};

	template <class ATTR_TYPE>
	class ConstPerMeshAttributeHandle{
	public:
		ConstPerMeshAttributeHandle(){_handle=nullptr;}
		ConstPerMeshAttributeHandle(const void *ah,const int & n):_handle ( (const Attribute<ATTR_TYPE> *)ah ),n_attr(n){}

		const Attribute<ATTR_TYPE> * _handle;
		int n_attr;
		const ATTR_TYPE & operator ()(){ return *((const ATTR_TYPE*)(_handle->DataBegin()));}
	};

	// Some common Handle typedefs to simplify use
	typedef typename MeshType::template PerVertexAttributeHandle<ScalarType> PerVertexScalarHandle;
	typedef typename MeshType::template PerVertexAttributeHandle<int>        PerVertexIntHandle;
	typedef typename MeshType::template PerVertexAttributeHandle<bool>       PerVertexBoolHandle;
	typedef typename MeshType::template PerVertexAttributeHandle<CoordType>  PerVertexCoordHandle;

	typedef typename MeshType::template PerFaceAttributeHandle<ScalarType> PerFaceScalarHandle;
	typedef typename MeshType::template PerFaceAttributeHandle<int>        PerFaceIntHandle;
	typedef typename MeshType::template PerFaceAttributeHandle<bool>       PerFaceBoolHandle;
	typedef typename MeshType::template PerFaceAttributeHandle<CoordType>  PerFaceCoordHandle;

	typedef typename MeshType::template PerTetraAttributeHandle<ScalarType> PerTetraScalarHandle;
	typedef typename MeshType::template PerTetraAttributeHandle<int> PerTetraIntHandle;
	typedef typename MeshType::template PerTetraAttributeHandle<bool> PerTetraBoolHandle;
	typedef typename MeshType::template PerTetraAttributeHandle<CoordType> PerTetraCoordHandle;


	// the camera member (that should keep the intrinsics) is no more needed since 2006, when intrisncs moved into the Shot structure
	//Camera<ScalarType> camera; // intrinsic
	Shot<ScalarType> shot;		// intrinsic && extrinsic

private:
	/// The per-mesh color. Not very useful and meaningful...
	Color4b c;
public:

	inline const Color4b &C() const	{ return c; }
	inline       Color4b &C()       { return c;  }
	inline       Color4b cC() const { return c;  }

	/// Default constructor
	TriMesh()
	{
		Clear();
	}

	/// destructor
	virtual ~TriMesh()
	{
//		ClearAttributes();
		Clear();
	}

	int Mem(const int & nv, const int & nf, const int & nt) const  {
		typename std::set< PointerToAttribute>::const_iterator i;
		int size = 0;
		size += sizeof(TriMesh)+sizeof(VertexType)*nv+sizeof(FaceType)*nf;

		for( i = vert_attr.begin(); i != vert_attr.end(); ++i)
			size += ((SimpleTempDataBase*)(*i)._handle)->SizeOf()*nv;
		for( i = edge_attr.begin(); i != edge_attr.end(); ++i)
			size += ((SimpleTempDataBase*)(*i)._handle)->SizeOf()*en;
		for( i = face_attr.begin(); i != face_attr.end(); ++i)
			size +=  ((SimpleTempDataBase*)(*i)._handle)->SizeOf()*nf;
		for (i = tetra_attr.begin(); i != tetra_attr.end(); ++i)
			size += ((SimpleTempDataBase *)(*i)._handle)->SizeOf() * nt;
		for( i = mesh_attr.begin(); i != mesh_attr.end(); ++i)
			size +=  ((SimpleTempDataBase*)(*i)._handle)->SizeOf();

		return size;
	}
	int MemUsed() const  {return Mem(vert.size(),face.size(), tetra.size());}
	inline int MemNeeded() const {return Mem(vn,fn);}



	/// Function to destroy the mesh
	void Clear()
	{
		for(FaceIterator fi = face.begin(); fi != face.end(); ++fi)
			(*fi).Dealloc();
		vert.clear();
		face.clear();
		edge.clear();
		tetra.clear();
		textures.clear();
		normalmaps.clear();
		vn = 0;
		en = 0;
		fn = 0;
		hn = 0;
		tn = 0;
		attrn = 0;
		imark = 0;
		C()=Color4b::Gray;

		typename std::set< PointerToAttribute>::iterator i;
		for (i = vert_attr.begin(); i != vert_attr.end(); ++i)
			((SimpleTempDataBase*)(*i)._handle)->Resize(0);

		for (i = edge_attr.begin(); i != edge_attr.end(); ++i)
			((SimpleTempDataBase*)(*i)._handle)->Resize(0);

		for (i = face_attr.begin(); i != face_attr.end(); ++i)
			((SimpleTempDataBase*)(*i)._handle)->Resize(0);

		for (i = tetra_attr.begin(); i != tetra_attr.end(); ++i)
			((SimpleTempDataBase*)(*i)._handle)->Resize(0);

	}


	void ClearAttributes()
	{
		// Clear attributes
		typename std::set< PointerToAttribute>::iterator i;
		for (i = vert_attr.begin(); i != vert_attr.end(); ++i)
			delete ((SimpleTempDataBase*)(*i)._handle);
		vert_attr.clear();

		for (i = edge_attr.begin(); i != edge_attr.end(); ++i)
			delete ((SimpleTempDataBase*)(*i)._handle);
		edge_attr.clear();

		for (i = face_attr.begin(); i != face_attr.end(); ++i)
			delete ((SimpleTempDataBase*)(*i)._handle);
		face_attr.clear();

		for (i = tetra_attr.begin(); i != tetra_attr.end(); ++i)
			delete ((SimpleTempDataBase *)(*i)._handle);
		tetra_attr.clear();

		for (i = mesh_attr.begin(); i != mesh_attr.end(); ++i)
			delete ((SimpleTempDataBase*)(*i)._handle);
		mesh_attr.clear();
		attrn = 0;
	}

	bool IsEmpty() const
	{
		return vert.empty() && edge.empty() && face.empty() && tetra.empty();
	}

	int & SimplexNumber(){ return fn;}
	int & VertexNumber(){ return vn;}

	/// The incremental mark
	int imark;

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
/// Initialize the imark-system of the edges
template <class MeshType> inline  void InitEdgeIMark(MeshType & m)
{
	typename MeshType::EdgeIterator ei;

	for (ei = m.edge.begin(); ei != m.edge.end(); ++ei)
		if( !(*ei).IsD() && (*ei).IsRW() )
			(*ei).InitIMark();
}

///initialize the imark-sysyem of the tetras
template <class MeshType>
inline void InitTetraIMark(MeshType &m)
{
	typename MeshType::TetraIterator ti;

	for (ti = m.tetra.begin(); ti != m.tetra.end(); ++ti)
		if (!(*ti).IsD() && (*ti).IsRW())
			(*ti).InitIMark();
}

/** \brief Access function to the incremental mark.
  You should not use this member directly. In most of the case just use IsMarked() and Mark()
*/
template <class MeshType> inline int & IMark(MeshType & m){return m.imark;}

/** \brief Check if the vertex incremental mark matches the one of the mesh.
	@param m the mesh containing the element
	@param v Vertex pointer */
template <class MeshType> inline bool IsMarked(const MeshType & m, typename MeshType::ConstVertexPointer  v )  { return v->cIMark() == m.imark; }

/** \brief Check if the edge incremental mark matches the one of the mesh.
	@param m the mesh containing the element
	@param e edge pointer */
template <class MeshType> inline bool IsMarked(const MeshType & m, typename MeshType::ConstEdgePointer e )  { return e->cIMark() == m.imark; }

/** \brief Check if the face incremental mark matches the one of the mesh.
	@param m the mesh containing the element
	@param f Face pointer */
template <class MeshType> inline bool IsMarked(const MeshType & m,typename MeshType::ConstFacePointer f )  { return f->cIMark() == m.imark; }

/** \brief Check if the tetra incremental mark matches the one of the mesh.
	@param m the mesh containing the element
	@param t tetra pointer */
template <class MeshType>
inline bool IsMarked(const MeshType &m, typename MeshType::ConstTetraPointer t) { return t->cIMark() == m.imark; }

/** \brief Set the vertex incremental mark of the vertex to the one of the mesh.
	@param m the mesh containing the element
	@param v Vertex pointer */
template <class MeshType> inline void Mark(MeshType & m, typename MeshType::VertexPointer v )  { v->IMark() = m.imark; }

/** \brief Set the edge incremental mark of the edge to the one of the mesh.
	@param m the mesh containing the element
	@param e edge pointer */
template <class MeshType> inline void Mark(MeshType & m, typename MeshType::EdgePointer e )  { e->IMark() = m.imark; }

/** \brief Set the face incremental mark of the vertex to the one of the mesh.
	@param m the mesh containing the element
	@param f Vertex pointer */
template <class MeshType> inline void Mark(MeshType & m, typename MeshType::FacePointer f )  { f->IMark() = m.imark; }

/** \brief Set the tetra incremental mark to the one of the mesh.
	@param m the mesh containing the element
	@param t tetra pointer */
template <class MeshType>
inline void Mark(MeshType &m, typename MeshType::TetraPointer t) { t->IMark() = m.imark; }


/** \brief Unmark, in constant time, all the elements (face and vertices) of a mesh.
	@param m the mesh containing the element

	In practice this function just increment the internal counter that stores the value for which an element is considered marked;
	therefore all the mesh elements become immediately un-mmarked.
	*/
template <class MeshType> inline void UnMarkAll(MeshType & m)
{
	++m.imark;
}


//template < class CType0, class CType1 , class CType2, class CType3>
//bool HasPerVertexVEAdjacency (const TriMesh < CType0, CType1, CType2, CType3> & /*m*/) {return TriMesh < CType0 , CType1, CType2, CType3>::VertContainer::value_type::HasVEAdjacency();}
//template < class  CType0, class CType1, class CType2 , class CType3>
//bool HasPerEdgeVEAdjacency   (const TriMesh < CType0, CType1, CType2, CType3> & /*m*/) {return TriMesh < CType0 , CType1, CType2, CType3>::EdgeContainer::value_type::HasVEAdjacency();}

template < class VertexType> bool VertexVectorHasVFAdjacency     (const std::vector<VertexType> &) {  return VertexType::HasVFAdjacency(); }
template < class VertexType> bool VertexVectorHasVEAdjacency     (const std::vector<VertexType> &) {  return VertexType::HasVEAdjacency(); }
template < class VertexType> bool VertexVectorHasVTAdjacency     (const std::vector<VertexType> &) { return VertexType::HasVTAdjacency(); }
template < class EdgeType  > bool   EdgeVectorHasVEAdjacency     (const std::vector<EdgeType  > &) {  return EdgeType::HasVEAdjacency(); }
template < class EdgeType  > bool   EdgeVectorHasEEAdjacency     (const std::vector<EdgeType> &)   {  return EdgeType::HasEEAdjacency(); }
template < class EdgeType  > bool   EdgeVectorHasEFAdjacency     (const std::vector<EdgeType> &)   {  return EdgeType::HasEFAdjacency(); }
template < class FaceType  > bool   FaceVectorHasVFAdjacency     (const std::vector<FaceType  > &) {  return FaceType::HasVFAdjacency(); }

template < class TriMeshType> bool HasPerVertexVFAdjacency     (const TriMeshType &m) { return tri::VertexVectorHasVFAdjacency(m.vert); }
template < class TriMeshType> bool HasPerVertexVEAdjacency     (const TriMeshType &m) { return tri::VertexVectorHasVEAdjacency(m.vert); }
template < class TriMeshType> bool HasPerVertexVTAdjacency     (const TriMeshType &m) { return tri::VertexVectorHasVTAdjacency(m.vert); }
template < class TriMeshType> bool   HasPerEdgeVEAdjacency     (const TriMeshType &m) { return tri::EdgeVectorHasVEAdjacency  (m.edge); }
template < class TriMeshType> bool   HasPerEdgeEFAdjacency     (const TriMeshType &m) { return tri::EdgeVectorHasEFAdjacency  (m.edge); }
template < class TriMeshType> bool   HasPerFaceVFAdjacency     (const TriMeshType &m) { return tri::FaceVectorHasVFAdjacency  (m.face); }


template < class VertexType> bool VertexVectorHasPerVertexQuality     (const std::vector<VertexType> &) {  return VertexType::HasQuality     (); }
template < class VertexType> bool VertexVectorHasPerVertexNormal      (const std::vector<VertexType> &) {  return VertexType::HasNormal      (); }
template < class VertexType> bool VertexVectorHasPerVertexColor       (const std::vector<VertexType> &) {  return VertexType::HasColor       (); }
template < class VertexType> bool VertexVectorHasPerVertexMark        (const std::vector<VertexType> &) {  return VertexType::HasMark        (); }
template < class VertexType> bool VertexVectorHasPerVertexFlags       (const std::vector<VertexType> &) {  return VertexType::HasFlags       (); }
template < class VertexType> bool VertexVectorHasPerVertexRadius      (const std::vector<VertexType> &) {  return VertexType::HasRadius      (); }
template < class VertexType> bool VertexVectorHasPerVertexCurvatureDir(const std::vector<VertexType> &) {  return VertexType::HasCurvatureDir(); }
template < class VertexType> bool VertexVectorHasPerVertexTexCoord    (const std::vector<VertexType> &) {  return VertexType::HasTexCoord    (); }

template < class TriMeshType> bool HasPerVertexQuality     (const TriMeshType &m) { return tri::VertexVectorHasPerVertexQuality     (m.vert); }
template < class TriMeshType> bool HasPerVertexNormal      (const TriMeshType &m) { return tri::VertexVectorHasPerVertexNormal      (m.vert); }
template < class TriMeshType> bool HasPerVertexColor       (const TriMeshType &m) { return tri::VertexVectorHasPerVertexColor       (m.vert); }
template < class TriMeshType> bool HasPerVertexMark        (const TriMeshType &m) { return tri::VertexVectorHasPerVertexMark        (m.vert); }
template < class TriMeshType> bool HasPerVertexFlags       (const TriMeshType &m) { return tri::VertexVectorHasPerVertexFlags       (m.vert); }
template < class TriMeshType> bool HasPerVertexRadius      (const TriMeshType &m) { return tri::VertexVectorHasPerVertexRadius      (m.vert); }
template < class TriMeshType> bool HasPerVertexCurvatureDir(const TriMeshType &m) { return tri::VertexVectorHasPerVertexCurvatureDir(m.vert); }
template < class TriMeshType> bool HasPerVertexTexCoord    (const TriMeshType &m) { return tri::VertexVectorHasPerVertexTexCoord    (m.vert); }

template < class EdgeType> bool EdgeVectorHasPerEdgeQuality     (const std::vector<EdgeType> &) {  return EdgeType::HasQuality     (); }
template < class EdgeType> bool EdgeVectorHasPerEdgeNormal      (const std::vector<EdgeType> &) {  return EdgeType::HasNormal      (); }
template < class EdgeType> bool EdgeVectorHasPerEdgeColor       (const std::vector<EdgeType> &) {  return EdgeType::HasColor       (); }
template < class EdgeType> bool EdgeVectorHasPerEdgeMark        (const std::vector<EdgeType> &) {  return EdgeType::HasMark        (); }
template < class EdgeType> bool EdgeVectorHasPerEdgeFlags       (const std::vector<EdgeType> &) {  return EdgeType::HasFlags       (); }

template < class TriMeshType> bool HasPerEdgeQuality     (const TriMeshType &m) { return tri::EdgeVectorHasPerEdgeQuality     (m.edge); }
template < class TriMeshType> bool HasPerEdgeNormal      (const TriMeshType &m) { return tri::EdgeVectorHasPerEdgeNormal      (m.edge); }
template < class TriMeshType> bool HasPerEdgeColor       (const TriMeshType &m) { return tri::EdgeVectorHasPerEdgeColor       (m.edge); }
template < class TriMeshType> bool HasPerEdgeMark        (const TriMeshType &m) { return tri::EdgeVectorHasPerEdgeMark        (m.edge); }
template < class TriMeshType> bool HasPerEdgeFlags       (const TriMeshType &m) { return tri::EdgeVectorHasPerEdgeFlags       (m.edge); }


template < class FaceType>    bool FaceVectorHasPerWedgeColor   (const std::vector<FaceType> &) {  return FaceType::HasWedgeColor   (); }
template < class FaceType>    bool FaceVectorHasPerWedgeNormal  (const std::vector<FaceType> &) {  return FaceType::HasWedgeNormal  (); }
template < class FaceType>    bool FaceVectorHasPerWedgeTexCoord(const std::vector<FaceType> &) {  return FaceType::HasWedgeTexCoord(); }

template < class TriMeshType> bool HasPerWedgeColor   (const TriMeshType &m) { return tri::FaceVectorHasPerWedgeColor   (m.face); }
template < class TriMeshType> bool HasPerWedgeNormal  (const TriMeshType &m) { return tri::FaceVectorHasPerWedgeNormal  (m.face); }
template < class TriMeshType> bool HasPerWedgeTexCoord(const TriMeshType &m) { return tri::FaceVectorHasPerWedgeTexCoord(m.face); }

template < class  CType0, class CType1, class CType2 , class CType3>
bool HasPolyInfo (const TriMesh < CType0, CType1, CType2, CType3> & /*m*/) {return TriMesh < CType0 , CType1, CType2, CType3>::FaceContainer::value_type::HasPolyInfo();}

template < class FaceType>    bool FaceVectorHasPerFaceFlags  (const std::vector<FaceType> &) {  return FaceType::HasFlags  (); }
template < class FaceType>    bool FaceVectorHasPerFaceNormal (const std::vector<FaceType> &) {  return FaceType::HasNormal (); }
template < class FaceType>    bool FaceVectorHasPerFaceColor  (const std::vector<FaceType> &) {  return FaceType::HasColor  (); }
template < class FaceType>    bool FaceVectorHasPerFaceMark   (const std::vector<FaceType> &) {  return FaceType::HasMark   (); }
template < class FaceType>    bool FaceVectorHasPerFaceQuality(const std::vector<FaceType> &) {  return FaceType::HasQuality(); }
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


template < class TetraType> bool TetraVectorHasPerTetraFlags  (const std::vector<TetraType> &) { return TetraType::HasFlags  (); }
template < class TetraType> bool TetraVectorHasPerTetraColor  (const std::vector<TetraType> &) { return TetraType::HasColor  (); }
template < class TetraType> bool TetraVectorHasPerTetraMark   (const std::vector<TetraType> &) { return TetraType::HasMark   (); }
template < class TetraType> bool TetraVectorHasPerTetraQuality(const std::vector<TetraType> &) { return TetraType::HasQuality(); }
template < class TetraType> bool TetraVectorHasVTAdjacency   (const std::vector<TetraType> &) { return TetraType::HasVTAdjacency(); }
template < class TetraType> bool TetraVectorHasTTAdjacency   (const std::vector<TetraType> &) { return TetraType::HasTTAdjacency(); }

template < class TriMeshType> bool HasPerTetraFlags       (const TriMeshType &m) { return tri::TetraVectorHasPerTetraFlags       (m.tetra); }
template < class TriMeshType> bool HasPerTetraColor       (const TriMeshType &m) { return tri::TetraVectorHasPerTetraColor       (m.tetra); }
template < class TriMeshType> bool HasPerTetraMark        (const TriMeshType &m) { return tri::TetraVectorHasPerTetraMark        (m.tetra); }
template < class TriMeshType> bool HasPerTetraQuality     (const TriMeshType &m) { return tri::TetraVectorHasPerTetraQuality     (m.tetra); }

template < class TriMeshType> bool HasFFAdjacency   (const TriMeshType &m) { return tri::FaceVectorHasFFAdjacency   (m.face); }
template < class TriMeshType> bool HasEEAdjacency   (const TriMeshType &m) { return tri::EdgeVectorHasEEAdjacency   (m.edge); }
template < class TriMeshType> bool HasFEAdjacency   (const TriMeshType &m) { return tri::FaceVectorHasFEAdjacency   (m.face); }
template < class TriMeshType> bool HasEFAdjacency   (const TriMeshType &m) { return tri::FaceVectorHasFEAdjacency(m.face) && tri::EdgeVectorHasEFAdjacency(m.edge); }
template < class TriMeshType> bool HasFVAdjacency   (const TriMeshType &m) { return tri::FaceVectorHasFVAdjacency   (m.face); }
template < class TriMeshType> bool HasVFAdjacency   (const TriMeshType &m) { return tri::FaceVectorHasVFAdjacency   (m.face) && tri::VertexVectorHasVFAdjacency(m.vert);  }
template < class TriMeshType> bool HasVEAdjacency   (const TriMeshType &m) { return tri::EdgeVectorHasVEAdjacency   (m.edge) && tri::VertexVectorHasVEAdjacency(m.vert);  }
template < class TriMeshType> bool HasTVAdjacency   (const TriMeshType &m) { return tri::TetraVectorHasVTAdjacency(m.tetra); }
template < class TriMeshType> bool HasVTAdjacency   (const TriMeshType &m) { return tri::VertexVectorHasVTAdjacency(m.vert) && tri::TetraVectorHasVTAdjacency(m.tetra); }
template < class TriMeshType> bool HasTTAdjacency   (const TriMeshType &m) { return tri::TetraVectorHasTTAdjacency(m.tetra); }


//template < class  CType0, class CType1, class CType2 , class CType3>
//bool HasVEAdjacency (const TriMesh < CType0, CType1, CType2, CType3> & /*m*/) {return TriMesh < CType0 , CType1, CType2, CType3>::VertContainer::value_type::HasVEAdjacency();}

template < class  CType0, class CType1, class CType2 , class CType3>
bool HasVHAdjacency (const TriMesh < CType0, CType1, CType2, CType3> & /*m*/) {return TriMesh < CType0 , CType1, CType2, CType3>::VertContainer::value_type::HasVHAdjacency();}

template < class  CType0, class CType1, class CType2 , class CType3>
bool HasEVAdjacency (const TriMesh < CType0, CType1, CType2, CType3> & /*m*/) {return TriMesh < CType0 , CType1, CType2, CType3>::EdgeType::HasEVAdjacency();}

//template < class  CType0, class CType1, class CType2 , class CType3>
//bool HasEEAdjacency (const TriMesh < CType0, CType1, CType2, CType3> & /*m*/) {return TriMesh < CType0 , CType1, CType2, CType3>::EdgeType::HasEEAdjacency();}

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

template <class MeshType>
bool HasPerVertexAttribute(const MeshType &m,   std::string   name){
	typename std::set< typename MeshType::PointerToAttribute>::const_iterator ai;
	typename MeshType::PointerToAttribute h;
	h._name = name;
	ai = m.vert_attr.find(h);
	return (ai!= m.vert_attr.end() ) ;
}
template <class MeshType>
bool HasPerFaceAttribute(const MeshType &m,   std::string   name){
	typename std::set< typename MeshType::PointerToAttribute>::const_iterator ai;
	typename MeshType::PointerToAttribute h;
	h._name = name;
	ai = m.face_attr.find(h);
	return (ai!= m.face_attr.end() ) ;
}

template <class MeshType>
bool HasPerTetraAttribute(const MeshType &m, std::string name)
{
	typename std::set<typename MeshType::PointerToAttribute>::const_iterator ai;
	typename MeshType::PointerToAttribute h;
	h._name = name;
	ai = m.tetra_attr.find(h);
	return (ai != m.tetra_attr.end());
}

template <class MeshType>
bool HasPerMeshAttribute(const MeshType &m,   std::string   name){
	typename std::set< typename MeshType::PointerToAttribute>::const_iterator ai;
	typename MeshType::PointerToAttribute h;
	h._name = name;
	ai = m.mesh_attr.find(h);
	return (ai!= m.mesh_attr.end() ) ;
}

template <class MeshType> void RequireVertexCompactness (const MeshType &m) {
	if(m.vert.size()!=size_t(m.vn)) throw vcg::MissingCompactnessException("Vertex Vector Contains deleted elements");
}
template <class MeshType> void RequireFaceCompactness   (const MeshType &m) {
	if(m.face.size()!=size_t(m.fn)) throw vcg::MissingCompactnessException("Face Vector Contains deleted elements");
}
template <class MeshType> void RequireEdgeCompactness   (const MeshType &m) {
	if(m.edge.size()!=size_t(m.en)) throw vcg::MissingCompactnessException("Edge Vector Contains deleted elements");
}
template <class MeshType> void RequireTetraCompactness(const MeshType &m) {
	if (m.tetra.size() != size_t(m.tn)) throw vcg::MissingCompactnessException("Tetra Vector Contains deleted elements");
}

template <class MeshType>
void RequireCompactness(const MeshType &m)
{
	RequireVertexCompactness<MeshType>(m);
	RequireFaceCompactness<MeshType>(m);
	RequireEdgeCompactness<MeshType>(m);
	RequireTetraCompactness<MeshType>(m);
}

//todo require tetramesh
template <class MeshType> void RequireTriangularMesh (const MeshType &m ) { if( tri::HasPolyInfo( m ) ) throw vcg::MissingTriangularRequirementException("");}
template <class MeshType> void RequirePolygonalMesh (const MeshType &m )  { if(!tri::HasPolyInfo( m ) ) throw vcg::MissingPolygonalRequirementException("");}

template <class MeshType> void RequireVFAdjacency    (const MeshType &m) { if(!tri::HasVFAdjacency   (m)) throw vcg::MissingComponentException("VFAdjacency"); }
template <class MeshType> void RequireVEAdjacency    (const MeshType &m) { if(!tri::HasVEAdjacency   (m)) throw vcg::MissingComponentException("VEAdjacency"); }
template <class MeshType> void RequireFFAdjacency    (const MeshType &m) { if(!tri::HasFFAdjacency   (m)) throw vcg::MissingComponentException("FFAdjacency"); }
template <class MeshType> void RequireEEAdjacency    (const MeshType &m) { if(!tri::HasEEAdjacency   (m)) throw vcg::MissingComponentException("EEAdjacency"); }
template <class MeshType> void RequireFEAdjacency    (const MeshType &m) { if(!tri::HasFEAdjacency   (m)) throw vcg::MissingComponentException("FEAdjacency"); }
template <class MeshType> void RequireFHAdjacency    (const MeshType &m) { if(!tri::HasFHAdjacency   (m)) throw vcg::MissingComponentException("FHAdjacency"); }
template <class MeshType> void RequireVTAdjacency    (const MeshType &m) { if(!tri::HasVTAdjacency   (m)) throw vcg::MissingComponentException("VTAdjacency"); }
template <class MeshType> void RequireTTAdjacency    (const MeshType &m) { if(!tri::HasTTAdjacency   (m)) throw vcg::MissingComponentException("TTAdjacency"); }

template <class MeshType> void RequirePerVertexQuality     (const MeshType &m) { if(!tri::HasPerVertexQuality     (m)) throw vcg::MissingComponentException("PerVertexQuality     "); }
template <class MeshType> void RequirePerVertexNormal      (const MeshType &m) { if(!tri::HasPerVertexNormal      (m)) throw vcg::MissingComponentException("PerVertexNormal      "); }
template <class MeshType> void RequirePerVertexColor       (const MeshType &m) { if(!tri::HasPerVertexColor       (m)) throw vcg::MissingComponentException("PerVertexColor       "); }
template <class MeshType> void RequirePerVertexMark        (const MeshType &m) { if(!tri::HasPerVertexMark        (m)) throw vcg::MissingComponentException("PerVertexMark        "); }
template <class MeshType> void RequirePerVertexFlags       (const MeshType &m) { if(!tri::HasPerVertexFlags       (m)) throw vcg::MissingComponentException("PerVertexFlags       "); }
template <class MeshType> void RequirePerVertexRadius      (const MeshType &m) { if(!tri::HasPerVertexRadius      (m)) throw vcg::MissingComponentException("PerVertexRadius      "); }
template <class MeshType> void RequirePerVertexCurvatureDir(const MeshType &m) { if(!tri::HasPerVertexCurvatureDir(m)) throw vcg::MissingComponentException("PerVertexCurvatureDir"); }
template <class MeshType> void RequirePerVertexTexCoord    (const MeshType &m) { if(!tri::HasPerVertexTexCoord    (m)) throw vcg::MissingComponentException("PerVertexTexCoord    "); }

template <class MeshType> void RequirePerEdgeQuality (const MeshType &m) { if(!tri::HasPerEdgeQuality (m)) throw vcg::MissingComponentException("PerEdgeQuality "); }
template <class MeshType> void RequirePerEdgeNormal  (const MeshType &m) { if(!tri::HasPerEdgeNormal  (m)) throw vcg::MissingComponentException("PerEdgeNormal  "); }
template <class MeshType> void RequirePerEdgeColor   (const MeshType &m) { if(!tri::HasPerEdgeColor   (m)) throw vcg::MissingComponentException("PerEdgeColor   "); }
template <class MeshType> void RequirePerEdgeMark    (const MeshType &m) { if(!tri::HasPerEdgeMark    (m)) throw vcg::MissingComponentException("PerEdgeMark    "); }
template <class MeshType> void RequirePerEdgeFlags   (const MeshType &m) { if(!tri::HasPerEdgeFlags   (m)) throw vcg::MissingComponentException("PerEdgeFlags   "); }

template <class MeshType> void RequirePerFaceFlags       (const MeshType &m) { if(!tri::HasPerFaceFlags       (m)) throw vcg::MissingComponentException("PerFaceFlags       "); }
template <class MeshType> void RequirePerFaceNormal      (const MeshType &m) { if(!tri::HasPerFaceNormal      (m)) throw vcg::MissingComponentException("PerFaceNormal      "); }
template <class MeshType> void RequirePerFaceColor       (const MeshType &m) { if(!tri::HasPerFaceColor       (m)) throw vcg::MissingComponentException("PerFaceColor       "); }
template <class MeshType> void RequirePerFaceMark        (const MeshType &m) { if(!tri::HasPerFaceMark        (m)) throw vcg::MissingComponentException("PerFaceMark        "); }
template <class MeshType> void RequirePerFaceQuality     (const MeshType &m) { if(!tri::HasPerFaceQuality     (m)) throw vcg::MissingComponentException("PerFaceQuality     "); }
template <class MeshType> void RequirePerFaceCurvatureDir(const MeshType &m) { if(!tri::HasPerFaceCurvatureDir(m)) throw vcg::MissingComponentException("PerFaceCurvatureDir"); }

template <class MeshType> void RequirePerFaceWedgeColor   (const MeshType &m) { if(!tri::HasPerWedgeColor   (m)) throw vcg::MissingComponentException("PerFaceWedgeColor   "); }
template <class MeshType> void RequirePerFaceWedgeNormal  (const MeshType &m) { if(!tri::HasPerWedgeNormal  (m)) throw vcg::MissingComponentException("PerFaceWedgeNormal  "); }
template <class MeshType> void RequirePerFaceWedgeTexCoord(const MeshType &m) { if(!tri::HasPerWedgeTexCoord(m)) throw vcg::MissingComponentException("PerFaceWedgeTexCoord"); }

template <class MeshType> void RequirePerTetraFlags       (const MeshType &m) { if(!tri::HasPerTetraFlags       (m)) throw vcg::MissingComponentException("PerTetraFlags       "); }
template <class MeshType> void RequirePerTetraColor       (const MeshType &m) { if(!tri::HasPerTetraColor       (m)) throw vcg::MissingComponentException("PerTetraColor       "); }
template <class MeshType> void RequirePerTetraMark        (const MeshType &m) { if(!tri::HasPerTetraMark        (m)) throw vcg::MissingComponentException("PerTetraMark        "); }
template <class MeshType> void RequirePerTetraQuality     (const MeshType &m) { if(!tri::HasPerTetraQuality     (m)) throw vcg::MissingComponentException("PerTetraQuality     "); }

template <class MeshType> void RequirePerVertexAttribute(const MeshType &m, const char *name) { if(!HasPerVertexAttribute(m,name)) throw vcg::MissingComponentException("PerVertex attribute"); }
template <class MeshType> void RequirePerEdgeAttribute(const MeshType &m, const char *name)   { if(!HasPerEdgeAttribute(m,name))   throw vcg::MissingComponentException("PerEdge attribute"); }
template <class MeshType> void RequirePerFaceAttribute(const MeshType &m, const char *name)   { if(!HasPerFaceAttribute(m,name))   throw vcg::MissingComponentException("PerFace attribute"); }
template <class MeshType> void RequirePerTetraAttribute(const MeshType &m, const char *name)  { if(!HasPerTetraAttribute(m,name))  throw vcg::MissingComponentException("PerTetra attribute"); }
template <class MeshType> void RequirePerMeshAttribute(const MeshType &m, const char *name)   { if(!HasPerMeshAttribute(m,name))   throw vcg::MissingComponentException("PerMesh attribute"); }

/*@}*/
/*@}*/
}	 // end namespace
}	 // end namespace




#endif // BASE_H
