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
#include <vcg/complex/used_types.h>
#include <vcg/container/derivation_chain.h>

/*
People should subclass his vertex class from these one...
*/

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
		@param FaceContainerType (Template Parameter) Specifies the type of the faces container any the face type.
 */


 /* MeshTypeHolder is a class which is used to define the types in the mesh
*/

		template <class TYPESPOOL>
		struct BaseMeshTypeHolder{

				typedef bool ScalarType;
				typedef std::vector< Vertex<TYPESPOOL> >	CONTV;
				typedef std::vector< Edge<TYPESPOOL> >		CONTE;
				typedef std::vector< Face<TYPESPOOL> >		CONTF;
//				typedef std::vector< HEdge<TYPESPOOL> >		CONTHE;

				typedef CONTV														VertContainer;
				typedef Vertex<TYPESPOOL>								VertexType;
				typedef VertexType *										VertexPointer;
				typedef const VertexType *							ConstVertexPointer;
				typedef bool														CoordType;
				typedef typename CONTV::iterator				VertexIterator;
				typedef typename CONTV::const_iterator	ConstVertexIterator;

				typedef CONTE														EdgeContainer;
				typedef typename CONTE::value_type			EdgeType;
				typedef typename CONTE::value_type*			EdgePointer;
				typedef typename CONTE::iterator				EdgeIterator;
				typedef typename CONTE::const_iterator	ConstEdgeIterator;

				typedef CONTF														FaceContainer;
				typedef typename CONTF::value_type			FaceType;
				typedef typename CONTF::const_iterator	ConstFaceIterator;
				typedef typename CONTF::iterator				FaceIterator;
				typedef typename CONTF::value_type*			FacePointer;
				typedef const typename CONTF::value_type*ConstFacePointer;

				//typedef CONTHE													HEdgeContainer;
				//typedef typename CONTHE::value_type			HEdgeType;
				//typedef typename CONTHE::value_type*		HEdgePointer;
				//typedef typename CONTHE::iterator				HEdgeIterator;
				//typedef typename CONTHE::const_iterator	ConstHEdgeIterator;


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


/*struct DummyContainer {};
template <class CONT> struct Deriver: public MeshTypeHolder<CONT, typename CONT::value_type::IAm>{};
template <> struct Deriver<DummyContainer>{}*/;

template <typename T, typename CONT> struct Der: public MeshTypeHolder<T,CONT, typename CONT::value_type::IAm>{};
struct DummyContainer{struct value_type{ typedef int IAm;}; };

template < class Container0 = DummyContainer, class Container1 = DummyContainer, class Container2 = DummyContainer/*, class Container3 = DummyContainer*/ >
class TriMesh
		: public  MArity3<   BaseMeshTypeHolder<typename Container0::value_type::TypesPool>, Container0, Der ,Container1, Der, Container2, Der/*, Container3, Der*/>{
	public:

		typedef typename TriMesh::ScalarType ScalarType;
		typedef typename TriMesh::VertContainer VertContainer;
		typedef typename TriMesh::EdgeContainer EdgeContainer;
		typedef typename TriMesh::FaceContainer FaceContainer;

		// types for vertex
		typedef typename TriMesh::VertexType VertexType;
		typedef typename TriMesh::VertexPointer VertexPointer;
		typedef typename TriMesh::ConstVertexPointer ConstVertexPointer;
		typedef typename TriMesh::CoordType CoordType;
		typedef typename TriMesh::VertexIterator VertexIterator;
		typedef typename TriMesh::ConstVertexIterator ConstVertexIterator;

		// types for edge
		typedef typename TriMesh::EdgeType EdgeType;
		typedef typename TriMesh::EdgePointer EdgePointer;
		typedef typename TriMesh::EdgeIterator EdgeIterator;
		typedef typename TriMesh::ConstEdgeIterator ConstEdgeIterator;

		//types for face
		typedef typename TriMesh::FaceType FaceType;
		typedef typename TriMesh::ConstFaceIterator ConstFaceIterator;
		typedef typename TriMesh::FaceIterator FaceIterator;
		typedef typename TriMesh::FacePointer FacePointer;
		typedef typename TriMesh::ConstFacePointer ConstFacePointer;


	typedef TriMesh<Container0, Container1,Container2/*,Container3*/> MeshType;

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
	/// Actual number of faces
	int en;
	/// Bounding box of the mesh
	Box3<ScalarType> bbox;
	
  /// Nomi di textures
	//
  std::vector<std::string> textures;
	//
  std::vector<std::string> normalmaps;

	int attrn;	// total numer of attribute created

	class PointerToAttribute{
	public:
		void * _handle;			// pointer to the SimpleTempData that stores the attribute
		std::string _name;		// name of the attribute
 		int _sizeof;			// size of the attribute type (used only with VMI loading)	
		int _padding;			// padding 	(used only with VMI loading)

		int n_attr;				// unique ID of the attribute
		void Resize(const int & sz){((SimpleTempDataBase<VertContainer>*)_handle)->Resize(sz);}
		void Reorder(std::vector<size_t> & newVertIndex){((SimpleTempDataBase<VertContainer>*)_handle)->Reorder(newVertIndex);}
                bool operator<(const  PointerToAttribute    b) const {	return(_name.empty()&&b._name.empty())?(_handle < b._handle):( _name < b._name);}
	};
	
	std::set< PointerToAttribute > vert_attr;
 	std::set< PointerToAttribute > edge_attr;
 	std::set< PointerToAttribute > face_attr;
	std::set< PointerToAttribute > mesh_attr;


	template <class ATTR_TYPE, class CONT>
	class AttributeHandle{
	public:
		AttributeHandle(){_handle=(SimpleTempData<CONT,ATTR_TYPE> *)NULL;}
		AttributeHandle( void *ah,const int & n):_handle ( (SimpleTempData<CONT,ATTR_TYPE> *)ah ),n_attr(n){}
		AttributeHandle operator = ( const CONT & pva){ 
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
    PerVertexAttributeHandle( void *ah,const int & n):AttributeHandle<ATTR_TYPE,VertContainer>(ah,n){};
	};

	template <class ATTR_TYPE>
	class PerFaceAttributeHandle: public AttributeHandle<ATTR_TYPE,FaceContainer>{
	public:
		PerFaceAttributeHandle():AttributeHandle<ATTR_TYPE,FaceContainer>(){}
		PerFaceAttributeHandle( void *ah,const int & n):AttributeHandle<ATTR_TYPE,FaceContainer>(ah,n){};
	};

	template <class ATTR_TYPE>
	class PerEdgeAttributeHandle:  public AttributeHandle<ATTR_TYPE,EdgeContainer>{
	public: 
		PerEdgeAttributeHandle():AttributeHandle<ATTR_TYPE,EdgeContainer>(){}
		PerEdgeAttributeHandle( void *ah,const int & n):AttributeHandle<ATTR_TYPE,EdgeContainer>(ah,n){};
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
		fn = vn = 0;
		imark = 0;
		attrn = 0;
		C()=Color4b::Gray;
	}

	/// destructor
	~TriMesh()
	{
		typename std::set< PointerToAttribute>::iterator i;
		for( i = vert_attr.begin(); i != vert_attr.end(); ++i) 
			delete ((SimpleTempDataBase<VertContainer>*)(*i)._handle);
		for( i = face_attr.begin(); i != face_attr.end(); ++i) 
			delete ((SimpleTempDataBase<FaceContainer>*)(*i)._handle);
		for( i = mesh_attr.begin(); i != mesh_attr.end(); ++i) 
			delete ((AttributeBase*)(*i)._handle);

		FaceIterator fi;
		for(fi = face.begin(); fi != face.end(); ++fi) (*fi).Dealloc();
	}

	 int Mem(const int & nv, const int & nf) const  {
		typename std::set< PointerToAttribute>::const_iterator i;
		int size = 0;
		size += sizeof(TriMesh)+sizeof(VertexType)*nv+sizeof(FaceType)*nf;

		for( i = vert_attr.begin(); i != vert_attr.end(); ++i) 
			size += ((SimpleTempDataBase<VertContainer>*)(*i)._handle)->SizeOf()*nv;
		for( i = face_attr.begin(); i != face_attr.end(); ++i) 
			size +=  ((SimpleTempDataBase<FaceContainer>*)(*i)._handle)->SizeOf()*nf;
		for( i = mesh_attr.begin(); i != mesh_attr.end(); ++i) 
			size +=  ((AttributeBase*)(*i)._handle)->SizeOf();

		return size;
	}
	int MemUsed() const  {return Mem(vert.size(),face.size());}
	inline int MemNeeded() const {return Mem(vn,fn);}



/// Function to destroy the mesh
void Clear()
{
	vert.clear();
	face.clear();
//	textures.clear();
//	normalmaps.clear();
	vn = 0;
	fn = 0;
}

/// Reflection functions that speak about vertex and face properties.
static bool HasPerVertexNormal()  { return VertexType::HasNormal() ; }
static bool HasPerVertexColor()   { return VertexType::HasColor()  ; }
static bool HasPerVertexMark()    { return VertexType::HasMark()   ; }
static bool HasPerVertexQuality() { return VertexType::HasQuality(); }
static bool HasPerVertexTexCoord(){ return VertexType::HasTexCoord(); }
static bool HasPerVertexRadius()  { return VertexType::HasRadius(); }
static bool HasPerVertexFlags()   { return VertexType::HasFlags(); }

static bool HasPolyInfo()					{ return FaceType::HasPolyInfo() ; }
static bool HasPerFaceColor()     { return FaceType::HasFaceColor() ; }
static bool HasPerFaceNormal()    { return FaceType::HasFaceNormal()  ; }
static bool HasPerFaceMark()      { return FaceType::HasFaceMark()   ; }
static bool HasPerFaceQuality()   { return FaceType::HasFaceQuality(); }
static bool HasPerFaceFlags()			{ return FaceType::HasFlags(); }

static bool HasPerWedgeColor()     { return FaceType::HasWedgeColor()  ; }
static bool HasPerWedgeNormal()    { return FaceType::HasWedgeNormal()  ; }
static bool HasPerWedgeMark()      { return FaceType::HasWedgeMark()   ; }
static bool HasPerWedgeQuality()   { return FaceType::HasWedgeQuality(); }
static bool HasPerWedgeTexCoord()  { return FaceType::HasWedgeTexCoord(); }

static bool HasFFTopology()       { return FaceType::HasFFAdjacency();  }
static bool HasVFTopology()       { return ((FaceType::HasVFAdjacency())&&(VertexType::HasVFAdjacency())); }
static bool HasTopology()         { return HasFFTopology() || HasVFTopology(); }

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
	// TriMesh cannot be copied. Use Append (see vcg/complex/trimesh/append.h)
	TriMesh operator =(const TriMesh &  m){assert(0);return TriMesh();}
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



template < class VertContainerType, class FaceContainerType , class EdgeContainerType>
bool HasPerVertexQuality (const TriMesh < VertContainerType , FaceContainerType, EdgeContainerType> & /*m*/) {return VertContainerType::value_type::HasQuality();}

template < class VertContainerType, class FaceContainerType , class EdgeContainerType>
bool HasPerVertexMark (const TriMesh < VertContainerType , FaceContainerType, EdgeContainerType> & /*m*/) {return VertContainerType::value_type::HasMark();}

template < class VertContainerType, class FaceContainerType, class EdgeContainerType >
bool HasPerVertexCurvature (const TriMesh < VertContainerType , FaceContainerType, EdgeContainerType> & /*m*/) {return VertContainerType::value_type::HasCurvature();}

template < class VertContainerType, class FaceContainerType , class EdgeContainerType>
bool HasPerVertexCurvatureDir (const TriMesh < VertContainerType , FaceContainerType, EdgeContainerType> & /*m*/) {return VertContainerType::value_type::HasCurvatureDir();}

template < class VertContainerType, class FaceContainerType , class EdgeContainerType>
bool HasPerVertexColor (const TriMesh < VertContainerType , FaceContainerType, EdgeContainerType> & /*m*/) {return VertContainerType::value_type::HasColor();}

template < class VertContainerType, class FaceContainerType, class EdgeContainerType >
bool HasPerVertexTexCoord (const TriMesh < VertContainerType , FaceContainerType, EdgeContainerType> & /*m*/) {return VertContainerType::value_type::HasTexCoord();}

template < class VertContainerType, class FaceContainerType, class EdgeContainerType >
bool HasPerVertexFlags (const TriMesh < VertContainerType , FaceContainerType, EdgeContainerType> & /*m*/) {return VertContainerType::value_type::HasFlags();}

template < class VertContainerType, class FaceContainerType, class EdgeContainerType >
bool HasPerVertexNormal (const TriMesh < VertContainerType , FaceContainerType, EdgeContainerType> & /*m*/) {return VertContainerType::value_type::HasNormal();}

template < class VertContainerType, class FaceContainerType, class EdgeContainerType >
bool HasPerVertexRadius (const TriMesh < VertContainerType , FaceContainerType, EdgeContainerType> & /*m*/) {return VertContainerType::value_type::HasRadius();}

template < class VertContainerType, class FaceContainerType, class EdgeContainerType >
bool HasPerWedgeTexCoord (const TriMesh < VertContainerType , FaceContainerType, EdgeContainerType> & /*m*/) {return FaceContainerType::value_type::HasWedgeTexCoord();}

template < class VertContainerType, class FaceContainerType, class EdgeContainerType >
bool HasPerWedgeNormal (const TriMesh < VertContainerType , FaceContainerType, EdgeContainerType> & /*m*/) {return FaceContainerType::value_type::HasWedgeNormal();}

template < class VertContainerType, class FaceContainerType, class EdgeContainerType >
bool HasPerWedgeColor (const TriMesh < VertContainerType , FaceContainerType, EdgeContainerType> & /*m*/) {return FaceContainerType::value_type::HasWedgeColor();}

template < class VertContainerType, class FaceContainerType, class EdgeContainerType >
bool HasPerFaceFlags (const TriMesh < VertContainerType , FaceContainerType, EdgeContainerType> & /*m*/) {return FaceContainerType::value_type::HasFlags();}

template < class VertContainerType, class FaceContainerType, class EdgeContainerType >
bool HasPerFaceNormal (const TriMesh < VertContainerType , FaceContainerType, EdgeContainerType> & /*m*/) {return FaceContainerType::value_type::HasFaceNormal();}

template < class VertContainerType, class FaceContainerType, class EdgeContainerType >
bool HasPerFaceColor (const TriMesh < VertContainerType , FaceContainerType, EdgeContainerType> & /*m*/) {return FaceContainerType::value_type::HasFaceColor();}

template < class VertContainerType, class FaceContainerType, class EdgeContainerType >
bool HasPerFaceMark (const TriMesh < VertContainerType , FaceContainerType, EdgeContainerType> & /*m*/) {return FaceContainerType::value_type::HasMark();}

template < class VertContainerType, class FaceContainerType, class EdgeContainerType >
bool HasPerFaceQuality (const TriMesh < VertContainerType , FaceContainerType, EdgeContainerType> & /*m*/) {return FaceContainerType::value_type::HasFaceQuality();}

template < class VertContainerType, class FaceContainerType, class EdgeContainerType >
bool HasFFAdjacency (const TriMesh < VertContainerType , FaceContainerType, EdgeContainerType> & /*m*/) {return FaceContainerType::value_type::HasFFAdjacency();}

template < class VertContainerType, class FaceContainerType, class EdgeContainerType >
bool HasFVAdjacency (const TriMesh < VertContainerType , FaceContainerType, EdgeContainerType> & /*m*/) {return FaceContainerType::value_type::HasFVAdjacency();}

template < class VertContainerType, class FaceContainerType, class EdgeContainerType >
bool HasVEAdjacency (const TriMesh < VertContainerType , FaceContainerType, EdgeContainerType> & /*m*/) {return VertContainerType::value_type::HasVEAdjacency();}

template < class VertContainerType, class FaceContainerType, class EdgeContainerType >
bool HasEVAdjacency (const TriMesh < VertContainerType , FaceContainerType, EdgeContainerType> & /*m*/) {return TriMesh < VertContainerType , FaceContainerType, EdgeContainerType>::EdgeType::HasEVAdjacency();}

template < class VertContainerType, class FaceContainerType, class EdgeContainerType >
bool HasEFAdjacency (const TriMesh < VertContainerType , FaceContainerType, EdgeContainerType> & /*m*/) {return TriMesh < VertContainerType , FaceContainerType, EdgeContainerType>::EdgeType::HasEFAdjacency();}

template < class VertContainerType, class FaceContainerType, class EdgeContainerType >
bool HasFHEAdjacency (const TriMesh < VertContainerType , FaceContainerType, EdgeContainerType> & /*m*/) {return TriMesh < VertContainerType , FaceContainerType, EdgeContainerType>::FaceType::HasFHEAdjacency();}

template < class VertContainerType, class FaceContainerType, class EdgeContainerType >
bool HasHEVAdjacency (const TriMesh < VertContainerType , FaceContainerType, EdgeContainerType> & /*m*/) {return TriMesh < VertContainerType , FaceContainerType, EdgeContainerType>::EdgeType::HasHEVAdjacency();}

template < class VertContainerType, class FaceContainerType, class EdgeContainerType >
bool HasHENextAdjacency (const TriMesh < VertContainerType , FaceContainerType, EdgeContainerType> & /*m*/) {return EdgeContainerType::value_type::HasHENextAdjacency();}

template < class VertContainerType, class FaceContainerType, class EdgeContainerType >
bool HasHEPrevAdjacency (const TriMesh < VertContainerType , FaceContainerType, EdgeContainerType> & /*m*/) {return EdgeContainerType::value_type::HasHEPrevAdjacency();}

template < class VertContainerType, class FaceContainerType, class EdgeContainerType >
bool HasHEOppAdjacency (const TriMesh < VertContainerType , FaceContainerType, EdgeContainerType> & /*m*/) {return EdgeContainerType::value_type::HasHEOppAdjacency();}

template < class VertContainerType, class FaceContainerType , class EdgeContainerType>
bool HasVFAdjacency (const TriMesh < VertContainerType , FaceContainerType,   EdgeContainerType> & /*m*/) {
  assert(FaceContainerType::value_type::HasVFAdjacency() == VertContainerType::value_type::HasVFAdjacency());
  return FaceContainerType::value_type::HasVFAdjacency();
}

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

