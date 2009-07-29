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
Revision 1.34  2008/05/16 10:07:36  ganovelli
added Trimesh destructor to purge unremoved PerVertex[PerFace]Attribute

Revision 1.33  2008/05/15 16:32:27  ganovelli
PerVertexAttribute and PerFaceAttribute added to Trimesh

Revision 1.32  2008/04/15 10:34:07  cignoni
added  HasPerVertexTexCoord ( mesh )

Revision 1.31  2008/02/21 17:27:06  cignoni
Added HasPerVertexColor static function

Revision 1.30  2008/01/28 14:46:03  cignoni
added hasPerWedgeColor and HasPerWedgeNormal

Revision 1.29  2008/01/28 08:42:07  cignoni
added HasPerFaceNormal and HasPerVertexNormal

Revision 1.28  2007/03/12 15:38:03  tarini
Texture coord name change!  "TCoord" and "Texture" are BAD. "TexCoord" is GOOD.

Revision 1.27  2007/02/22 09:18:41  cignoni
Added guards on msvc pragmas

Revision 1.26  2007/02/14 15:31:41  ganovelli
Added HasPerVertexFlag

Revision 1.25  2006/11/28 22:35:29  cignoni
Added Consistency check in the HasVFAdj static function

Revision 1.24  2006/11/07 11:29:23  cignoni
Corrected some errors in the reflections Has*** functions

Revision 1.23  2006/10/27 11:08:18  ganovelli
added override to HasFFAdjacency , HasVFAdjacency for the optional attributes (see also complex/trimesh/allocate.h)

Revision 1.22  2006/07/10 14:26:22  cignoni
Minor. Added a disambiguating this at the constructor of trimesh

Revision 1.21  2006/05/25 04:40:57  cignoni
Updated HasPerFaceColor/Quality to the new style with mesh param.

Revision 1.20  2006/05/03 21:35:31  cignoni
Added new style HasPerFaceColor(m) and HasPerFaceMark(m)

Revision 1.19  2005/12/02 00:05:34  cignoni
Added HasFlags and a couple of missing include files

Revision 1.18  2005/11/26 00:16:03  cignoni
added  HasPerWedgeTexture  taking mesh as input. (needed for optional components)

Revision 1.17  2005/11/16 22:35:47  cignoni
Added missing includes (color and assert)
Added texture name members

Revision 1.16  2005/11/15 12:09:17  rita_borgo
Changed Volume Routine, before was returning negative values

Revision 1.15  2005/10/03 16:00:08  rita_borgo
Minor changes

Revision 1.14  2005/03/18 16:37:46  fiorin
Minor changes

Revision 1.13  2005/01/28 17:56:57  pietroni
changed HasVFTopology function... control if both vertex and face define the vf topology

Revision 1.12  2004/10/28 00:54:34  cignoni
Better Doxygen documentation

Revision 1.11  2004/10/07 14:25:38  ganovelli
added camera and shot

Revision 1.10  2004/09/08 15:15:05  ganovelli
changes for gcc

Revision 1.9  2004/07/15 12:03:50  ganovelli
access to imark added

Revision 1.8  2004/07/15 11:39:24  ganovelli
IsDeleted to IsD

Revision 1.7  2004/07/09 10:18:19  ganovelli
added access functions to vn and fn

Revision 1.6  2004/05/04 02:29:54  ganovelli
removed  Const from ConstFacePointer and ConstVertexPointer in the arguement function Mark, which are meant to be changed

Revision 1.5  2004/03/18 16:00:10  cignoni
minor changes

Revision 1.4  2004/03/10 00:57:44  cignoni
minor changes

Revision 1.3  2004/03/07 21:54:56  cignoni
some more reflection functions

Revision 1.2  2004/03/04 00:08:15  cignoni
First working version!

Revision 1.1  2004/02/19 13:11:06  cignoni
Initial commit


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
#include <vcg/simplex/edge/base.h>
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

	// awful series of tricks to make the edge container as default parameter
	

	// TriMeshEdgeHolder is used to declare everything is needed for the edge type
	// If a real container of edge is passed then it simply  makes the handy  declarations
template < class VertContainerType, class FaceContainerType, class EdgeContainerType>
class TriMeshEdgeHolder{
	public:
		typedef EdgeContainerType EdgeContainer;
		typedef typename EdgeContainer::value_type EdgeType;
		typedef EdgeType*  EdgePointer;
		typedef const EdgeType * ConstEdgePointer;
	};

	// a dummy class is used  to provide the interface of a stl container
	class DummyClass:public std::vector<int>{};
	// If the DummyClass is passed it provides the interfaces to compile
template < class VertContainerType, class FaceContainerType >
class TriMeshEdgeHolder<VertContainerType,FaceContainerType,DummyClass>{
	public:
	class EdgeType: public EdgeSimp2<	typename VertContainerType::value_type,EdgeType,typename FaceContainerType::value_type >{};
	typedef typename FaceContainerType::value_type::EdgeType EdgeTypeExternal;
	struct EdgePointer {};	
	struct ConstEdgePointer {};	
	typedef std::vector< EdgeType > EdgeContainerType;
	typedef typename std::vector< EdgeType > EdgeContainer;
};



template < class VertContainerType, class FaceContainerType, class EdgeConts = DummyClass >
class TriMesh: public TriMeshEdgeHolder<VertContainerType,FaceContainerType,EdgeConts>{
	public:
	typedef typename TriMeshEdgeHolder<VertContainerType,FaceContainerType,EdgeConts>::EdgeContainer EdgeContainer;


	typedef TriMesh<VertContainerType, FaceContainerType> MeshType;
	typedef FaceContainerType FaceContainer;
	typedef VertContainerType VertContainer;
	typedef typename VertContainer::value_type VertexType;
	typedef typename FaceContainer::value_type FaceType;
	typedef typename VertexType::ScalarType ScalarType;
	typedef typename VertexType::CoordType CoordType;
	typedef typename VertContainer::iterator VertexIterator;
	typedef typename EdgeContainer::iterator EdgeIterator;
	typedef typename FaceContainer::iterator FaceIterator;
	typedef typename VertContainer::const_iterator ConstVertexIterator;
	typedef typename EdgeContainer::const_iterator ConstEdgeIterator;
	typedef typename FaceContainer::const_iterator ConstFaceIterator;
	typedef VertexType * VertexPointer;
	typedef const VertexType * ConstVertexPointer;
	typedef FaceType * FacePointer;
	typedef const FaceType * ConstFacePointer;
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

/// Initialize the imark-system of the faces
void InitFaceIMark()
{
	FaceIterator f;
	
	for(f=face.begin();f!=face.end();++f)
		if( !(*f).IsD() && (*f).IsR() && (*f).IsW() )
			(*f).InitIMark();
}

/// Initialize the imark-system of the vertices
void InitVertexIMark()
{
	VertexIterator vi;

	for(vi=vert.begin();vi!=vert.end();++vi)
		if( !(*vi).IsD() && (*vi).IsRW() )
			(*vi).InitIMark();
}

/// The incremental mark
int imark;
/** Access function to the incremental mark. 
*/
inline int & IMark(){return imark;}
/** Check if the vertex incremental mark matches the one of the mesh. 
	@param v Vertex pointer
*/
inline bool IsMarked( ConstVertexPointer  v ) const { return v->IMark() == imark; }
/** Check if the face incremental mark matches the one of the mesh. 
	@param v Face pointer
*/
inline bool IsMarked( ConstFacePointer f ) const { return f->IMark() == imark; }
/** Set the vertex incremental mark of the vertex to the one of the mesh.
	@param v Vertex pointer
*/
inline void Mark( VertexPointer v ) const { v->IMark() = imark; }
/** Set the face incremental mark of the vertex to the one of the mesh.
	@param v Vertex pointer
*/
inline void Mark( FacePointer f ) const { f->IMark() = imark; }
/// Unmark the mesh
inline void UnMarkAll() { ++imark; }



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

