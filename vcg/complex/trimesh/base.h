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
#ifndef __GNUC
#pragma warning( disable : 4804 )
#endif
#include <vcg/space/box3.h>
#include <vcg/math/shot.h>

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
template < class VertContainerType, class FaceContainerType >
class TriMesh{
	public:
	typedef FaceContainerType FaceContainer;
	typedef VertContainerType VertContainer;
	typedef typename VertContainer::value_type VertexType;
	typedef typename FaceContainer::value_type FaceType;
	typedef typename VertexType::ScalarType ScalarType;
	typedef typename VertexType::CoordType CoordType;
	typedef typename VertContainer::iterator VertexIterator;
	typedef typename FaceContainer::iterator FaceIterator;
	typedef typename VertContainer::const_iterator ConstVertexIterator;
	typedef typename FaceContainer::const_iterator ConstFaceIterator;
	typedef VertexType * VertexPointer;
	typedef const VertexType * ConstVertexPointer;
	typedef FaceType * FacePointer;
	typedef const FaceType * ConstFacePointer;
	typedef Box3<ScalarType> BoxType;

	/// Set of vertices 
	VertContainer vert;
	/// Real number of vertices
	int vn;
	/// Set of faces
	FaceContainer face;
	/// Real number of faces
	int fn;
	/// Bounding box of the mesh
	Box3<ScalarType> bbox;
	
  /// Nomi di textures
	//vector<string> textures;
	//vector<string> normalmaps;

		/// La camera
	Camera<ScalarType> camera; // intrinsic
	Shot<ScalarType> shot;		// extrinsic

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
	TriMesh():shot(camera)
	{
		fn = vn = 0;
		imark = 0;
	}

	inline int MemUsed() const
	{
		return sizeof(MMTYPE)+sizeof(MVTYPE)*vert.size()+sizeof(MFTYPE)*face.size();
	}

	inline int MemNeeded() const
	{
		return sizeof(MMTYPE)+sizeof(MVTYPE)*vn+sizeof(MFTYPE)*fn;
	}



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
static bool HasPerVertexTexture() { return VertexType::HasTexture(); }

static bool HasPerFaceColor()     { return FaceType::HasFaceColor() ; }
static bool HasPerFaceNormal()    { return FaceType::HasFaceNormal()  ; }
static bool HasPerFaceMark()      { return FaceType::HasFaceMark()   ; }
static bool HasPerFaceQuality()   { return FaceType::HasFaceQuality(); }

static bool HasPerWedgeColor()     { return FaceType::HasWedgeColor()  ; }
static bool HasPerWedgeNormal()    { return FaceType::HasWedgeNormal()  ; }
static bool HasPerWedgeMark()      { return FaceType::HasWedgeMark()   ; }
static bool HasPerWedgeQuality()   { return FaceType::HasWedgeQuality(); }
static bool HasPerWedgeTexture()   { return FaceType::HasWedgeTexture(); }

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
	    T = ( (*fi).P( k     ) - (*fi).P(j) ) ^
	        ( (*fi).P((k+1)%3) - (*fi).P(j) ) ;
	    B.Normalize();
	    N = T ^ B;
     
	    CoordType pj = (*fi).P(j);
	    CoordType pk = (*fi).P(k);
   

	    V +=  (pj*  T )*(pj*N)*(pj*B);
	    V +=  (pk*(-T))*(pk*N)*(pk*B);
     }
  }
	return V/6.0;
}


};	// end class Mesh
/*@}*/
/*@}*/
}	 // end namespace
}	 // end namespace


#endif

