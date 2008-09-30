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
Revision 1.7  2007/03/12 15:38:02  tarini
Texture coord name change!  "TCoord" and "Texture" are BAD. "TexCoord" is GOOD.

Revision 1.6  2005/09/14 14:09:16  spinelli
ConstVertexPointer --> VertexPointer
ConstEdgePointer --> EdgePointer

Revision 1.5  2005/05/17 21:19:37  ganovelli
some std::and typename  missing  (CRS4)

Revision 1.4  2004/10/28 00:47:42  cignoni
Better Doxygen documentation

Revision 1.3  2004/09/20 09:30:03  cignoni
Better Doxygen docs

Revision 1.2  2004/05/10 14:41:45  ganovelli
name of adhacency function updated

Revision 1.1  2004/04/26 19:10:04  ganovelli
created


****************************************************************************/

#pragma warning( disable : 4804 )

/*
People should subclass his vertex class from these one...
*/

#ifndef __VCGLIB_EDGEMESH
#define __VCGLIB_EDGEMESH

namespace vcg {
namespace edg {
/** \addtogroup edgemesh */
/*@{*/

/** \class EdgeMesh.
    This is class for definition of a mesh.
		@param VertContainerType (Template Parameter) Specifies the type of the vertices container any the vertex type.
		@param EdgeContainerType (Template Parameter) Specifies the type of the faces container any the face type.
 */
template < class VertContainerType, class EdgeContainerType >
class EdgeMesh{
	public:
	typedef EdgeContainerType EdgeContainer;
	typedef VertContainerType VertContainer;
	typedef typename VertContainer::value_type VertexType;
	typedef typename EdgeContainerType::value_type EdgeType;
	typedef typename VertexType::ScalarType ScalarType;
	typedef typename VertexType::CoordType CoordType;
	typedef typename VertContainer::iterator VertexIterator;
	typedef typename EdgeContainerType::iterator EdgeIterator;
	typedef typename VertContainer::const_iterator ConstVertexIterator;
	typedef typename EdgeContainerType::const_iterator ConstEdgeIterator;
	typedef VertexType * VertexPointer;
	typedef const VertexType * ConstVertexPointer;
	typedef EdgeType * EdgePointer;
	typedef const EdgeType * ConstEdgePointer;
	typedef Box3<ScalarType> BoxType;

	/// Set of vertices 
	VertContainer vert;
	/// Real number of vertices
	int vn;
	/// Set of faces
	EdgeContainer edges;
	/// Real number of faces
	int en;
	/// Bounding box of the mesh
	Box3<ScalarType> bbox;
	
  /// Nomi di textures
	//std::vector<string> textures;
	//std::vector<string> normalmaps;

		/// La camera
	//Camera<ScalarType> camera;

		/// Il colore della mesh
private:
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
	EdgeMesh()
	{
		en = vn = 0;
		imark = 0;
	}

	inline int MemUsed() const
	{
		return sizeof(EdgeMesh)+sizeof(VertexType)*vert.size()+sizeof(EdgeType)*edges.size();
	}

	inline int MemNeeded() const
	{
		return sizeof(EdgeMesh)+sizeof(VertexType)*vn+sizeof(EdgeType)*en;
	}



/// Function to destroy the mesh
void Clear()
{
	vert.clear();
	edges.clear();
//	textures.clear();
//	normalmaps.clear();
	vn = 0;
	en = 0;
}

/// Reflection functions that speak about vertex and face properties.
static bool HasPerVertexNormal()  { return VertexType::HasNormal() ; }
static bool HasPerVertexColor()   { return VertexType::HasColor()  ; }
static bool HasPerVertexMark()    { return VertexType::HasMark()   ; }
static bool HasPerVertexQuality() { return VertexType::HasQuality(); }
static bool HasPerVertexTexCoord(){ return VertexType::HasTexCoord(); }

static bool HasPerEdgeColor()     { return EdgeType::HasEdgeColor() ; }
static bool HasPerEdgeNormal()    { return EdgeType::HasEdgeNormal()  ; }
static bool HasPerEdgeMark()      { return EdgeType::HasEdgeMark()   ; }
static bool HasPerEdgeQuality()   { return EdgeType::HasEdgeQuality(); }

static bool HasEETopology()       { return EdgeType::HasEEAdjacency();  }
static bool HasVETopology()       { return EdgeType::HasVEAdjacency(); }
static bool HasTopology()         { return HasEETopology() || HasVETopology(); }


/// Initialize the imark-system of the faces
void InitEdgeIMark()
{
	EdgeIterator f;
	
	for(f=edges.begin();f!=edges.end();++f)
		if( !(*f).IsDeleted() && (*f).IsR() && (*f).IsW() )
			(*f).InitIMark();
}

/// Initialize the imark-system of the vertices
void InitVertexIMark()
{
	VertexIterator vi;

	for(vi=vert.begin();vi!=vert.end();++vi)
		if( !(*vi).IsDeleted() && (*vi).IsRW() )
			(*vi).InitIMark();
}

/// The incremental mark
int imark;

/** Check if the vertex incremental mark matches the one of the mesh. 
*/
inline bool IsMarked( ConstVertexPointer  v ) const { return v->IMark() == imark; }
/** Check if the face incremental mark matches the one of the mesh. 
*/
inline bool IsMarked( ConstEdgePointer f ) const { return f->IMark() == imark; }
/** Set the vertex incremental mark of the vertex to the one of the mesh.
*/
inline void Mark( VertexPointer v ) const { v->IMark() = imark; }
/** Set the face incremental mark of the vertex to the one of the mesh.
*/
inline void Mark( EdgePointer f ) const { f->IMark() = imark; }
/// Unmark the mesh
inline void UnMarkAll() { ++imark; }

};	// end class EdgeMesh

/*@}*/
}	 // end namespace
}	 // end namespace


#endif

