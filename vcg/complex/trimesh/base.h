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
Revision 1.2  2004/03/04 00:08:15  cignoni
First working version!

Revision 1.1  2004/02/19 13:11:06  cignoni
Initial commit


****************************************************************************/

#pragma warning( disable : 4804 )

/*
People should subclass his vertex class from these one...
*/

#ifndef __VCG_MESH
#define __VCG_MESH

namespace vcg {
namespace tri {

/** Class Mesh.
    This is class for definition of a mesh.
		@param VertContainer (Template Parameter) Specifies the type of the vertices container any the vertex type.
		@param STL_FACE_CONT (Template Parameter) Specifies the type of the faces container any the face type.
 */
template < class STL_VERT_CONT, class STL_FACE_CONT >
class TriMesh{
	public:
	typedef STL_FACE_CONT FaceContainer;
	typedef STL_VERT_CONT VertContainer;
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
	TriMesh()
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

static bool HasPerFaceColor()     { return FaceType::HasFaceNormal() ; }
static bool HasPerFaceNormal()    { return FaceType::HasFaceColor()  ; }
static bool HasPerFaceMark()      { return FaceType::HasFaceMark()   ; }
static bool HasPerFaceQuality()   { return FaceType::HasFaceQuality(); }

static bool HasPerWedgeColor()     { return FaceType::HasWedgeNormal() ; }
static bool HasPerWedgeNormal()    { return FaceType::HasWedgeColor()  ; }
static bool HasPerWedgeMark()      { return FaceType::HasWedgeMark()   ; }
static bool HasPerWedgeQuality()   { return FaceType::HasWedgeQuality(); }

static bool HasFFTopology()       { return FaceType::HasFFAdjacency();  }
static bool HasVFTopology()       { return FaceType::HasVFAdjacency(); }
static bool HasTopology()         { return HasFFTopology() || HasVFTopology(); }


/// Initialize the imark-system of the faces
void InitFaceIMark()
{
	face_iterator f;
	
	for(f=face.begin();f!=face.end();++f)
		if( !(*f).IsDeleted() && (*f).IsR() && (*f).IsW() )
			(*f).InitIMark();
}

/// Initialize the imark-system of the vertices
void InitVertexIMark()
{
	vertex_iterator vi;

	for(vi=vert.begin();vi!=vert.end();++vi)
		if( !(*vi).IsDeleted() && (*vi).IsRW() )
			(*vi).InitIMark();
}

/// The incremental mark
int imark;

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
inline void Mark( ConstVertexPointer v ) const { v->IMark() = imark; }
/** Set the face incremental mark of the vertex to the one of the mesh.
	@param v Vertex pointer
*/
inline void Mark( ConstFacePointer f ) const { f->IMark() = imark; }
/// Unmark the mesh
inline void UnMarkAll() { ++imark; }



/// Calcolo del volume di una mesh chiusa
ScalarType Volume()
{
 
  face_iterator f;
  int j,k;
  ScalarType V = 0;
  vectorial_type T,N,B;
 
  for(f = face.begin(); f!=face.end(); ++f)
  {
	for(j = 0; j < 3; ++j)
	{
	  /*calcolo tangente, normale e binormale (6 volte)*/
	  k = (j+1)%3;
	  T = (*f).V(k)->P() - (*f).V(j)->P();
	  T.Normalize();
	  T = ( (*f).V( k     )->P() - (*f).V(j)->P() ) ^
	   ( (*f).V((k+1)%3)->P() - (*f).V(j)->P() ) ;
	  B.Normalize();
	  N = T ^ B;
   
	  vectorial_type pj = (*f).V(j)->P();
	  vectorial_type pk = (*f).V(k)->P();
 

	  V +=  (pj*  T )*(pj*N)*(pj*B);
	  V +=  (pk*(-T))*(pk*N)*(pk*B);
    }
  }
	return V/6;
}


};	// end class Mesh


}	 // end namespace
}	 // end namespace


#endif

