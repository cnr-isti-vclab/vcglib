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
Revision 1.13  2007/03/12 15:38:03  tarini
Texture coord name change!  "TCoord" and "Texture" are BAD. "TexCoord" is GOOD.

Revision 1.12  2005/12/12 11:10:35  ganovelli
modifications to compile with gcc

Revision 1.11  2005/01/12 11:06:54  ganovelli
added InitVertexIMark

Revision 1.10  2004/10/11 17:44:07  ganovelli
added include of color4

Revision 1.9  2004/07/15 00:16:37  cignoni
Better doxigen documentation

Revision 1.8  2004/07/09 10:18:49  ganovelli
added access functions to vn and tn

Revision 1.7  2004/06/25 11:27:21  pietroni
added function to access temporary mark for decimation

Revision 1.6  2004/06/01 17:17:29  ganovelli
pragma once removed ,
load ts removed to be put in io_tetramesh
CLear() added

Revision 1.5  2004/05/13 12:16:12  pietroni
first version... add vertex to mesh

Revision 1.4  2004/05/13 07:41:47  turini
Chenged #include <space\\box3.h> in #include <vcg\\space\\box3.h>

Revision 1.3  2004/05/06 10:57:49  pietroni
changed names to topology functions

Revision 1.2  2004/04/28 11:37:14  pietroni
*** empty log message ***

Revision 1.1  2004/04/20 12:41:39  pietroni
*** empty log message ***

Revision 1.1  2004/04/15 08:54:20  pietroni
*** empty log message ***


***************************************************************************/


#ifndef __VCG_TETRAMESH
#define __VCG_TETRAMESH
#include <vcg/space/box3.h>
#include <vcg/space/color4.h>


namespace vcg {
namespace tetra {
 /** \addtogroup tetramesh */
/*@{*/

  /**  Class TetraMesh.
 This is class for definition of a mesh.
		@param STL_VERT_CONT (Template Parameter) Specifies the type of the vertices container any the vertex type.
		@param STL_FACE_CONT (Template Parameter) Specifies the type of the faces container any the face type.
 */


template < class STL_VERT_CONT ,class STL_TETRA_CONT >
class Tetramesh{
	public:

/***********************************************/
/** @name Tetramesh Type Definitions **/
//@{
  
  /// The mesh type
	typedef Tetramesh<STL_VERT_CONT,STL_TETRA_CONT> TetraMeshType;

	/// The vertex container
	typedef STL_VERT_CONT VertexContainer;

	/// The tethaedhron container
	typedef STL_TETRA_CONT TetraContainer;

	/// The vertex type 
	typedef typename STL_VERT_CONT::value_type VertexType;
	
	/// The tetrahedron type 
	typedef typename STL_TETRA_CONT::value_type TetraType;

	/// The type of vertex iterator
	typedef typename STL_VERT_CONT::iterator VertexIterator;

	/// The type of tetra iterator
	typedef typename STL_TETRA_CONT::iterator TetraIterator;

	/// The type of constant vertex iterator
	typedef typename STL_VERT_CONT::const_iterator const_VertexIterator;

	/// The type of constant face iterator
	typedef typename STL_TETRA_CONT::const_iterator const_TetraIterator;

	/// The vertex pointer type
	typedef VertexType * VertexPointer;

	/// The tetra pointer type
	typedef TetraType * TetraPointer;

	/// The type of the constant vertex pointer
	typedef const VertexType * const_VertexPointer;

	/// The type of the constant tetrahedron pointer
	typedef const VertexType * const_TetraPointer;

	typedef typename VertexType::ScalarType ScalarType;
//@}

/***********************************************/
/** @Common Attributes of a tetrahedral mesh **/
//@{

	///temporary mark for decimation
	int IMark;

	/// Set of vertices 
	STL_VERT_CONT vert;

	/// Real number of vertices
	int vn;

	/// Set of tetrahedron
	STL_TETRA_CONT tetra;

	/// Real number of tetrahedron
	int tn;
  
  /// Real number of edges
	int en;

  ///Boundingbox della mesh
  Box3<ScalarType> bbox;
//@}
 
/***********************************************/
/** @Default Functions **/
//@{

	/// Default constructor
	Tetramesh()
	{   
		tn = vn = en = 0;
	}
	
	Tetramesh(VertexContainer v,TetraContainer t)
	{
		this->vert=v;
		this->tetra=t;
		vn=v.size();
		tn=t.size();
	}

	inline int MemUsed() const
	{
		return sizeof(Tetramesh)+sizeof(VertexType)*vert.size()+sizeof(TetraType)*tetra.size();
	}

	void Clear(){
		vert.clear();
		tetra.clear();
		tn = 0;
		vn = 0;
		}

	/// Initialize the imark-system of the vertices
	void InitVertexIMark()	
	{
		VertexIterator vi;

		for(vi=vert.begin();vi!=vert.end();++vi)
			if( !(*vi).IsD() && (*vi).IsRW() )
				(*vi).InitIMark();
}
//@}

/***********************************************/
/** @Functions used to retrieve informations**/
//@{
 /// Reflection functions that speak about vertex and face properties.
static bool HasPerVertexNormal()  { return VertexType::HasNormal() ; }
static bool HasPerVertexColor()   { return VertexType::HasColor()  ; }
static bool HasPerVertexMark()    { return VertexType::HasMark()   ; }
static bool HasPerVertexQuality() { return VertexType::HasQuality(); }
static bool HasPerVertexTexCoord(){ return VertexType::HasTexCoord(); }

static bool HasPerTetraNormal()    { return TetraType::HasTetraNormal()  ; }
static bool HasPerTetraMark()      { return TetraType::HasTetraMark()   ; }
static bool HasPerTetraQuality()   { return TetraType::HasTetraQuality(); }

static bool HasTTTopology()       { return TetraType::HasTTAdjacency();  }
static bool HasVTTopology()       { return TetraType::HasVTAdjacency(); }
static bool HasTopology()         { return HasTTTopology() || HasVTTopology(); }

int & SimplexNumber(){ return tn;}
int & VertexNumber(){ return vn;}
/***********************************************/

/** @Functions used for handle the temporany mark of a tetrahedron used in decimation**/
//@{

///Increase the current mark.
	void UnMarkAll()
	{	
		++IMark;
	}

///Mark the vertex with current value
	void Mark(VertexType *v)
	{
		 v->IMark()=IMark;
	}

  ///return the current mark
	int GetMark()
	{
		return (IMark);
	}

///Initialize the mark of all vertices
	void InitIMark()
	{
	VertexIterator vi;
	IMark=0;
	for(vi=vert.begin();vi!=vert.end();vi++)
	{
		(*vi).InitIMark();
	}
	}

//@}
};//End class

/*@}*/


};//end namespace
};//end namespace
#endif

