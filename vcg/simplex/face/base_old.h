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
Revision 1.43  2007/03/12 15:37:19  tarini
Texture coord name change!  "TCoord" and "Texture" are BAD. "TexCoord" is GOOD.

Revision 1.42  2007/02/20 14:07:53  ganovelli
added QualityType to comply faceplus type

Revision 1.41  2007/01/13 00:25:36  cignoni
Added #ifdefs to guarantee that ComputeNormal would be defined only once

Revision 1.40  2007/01/11 10:13:11  cignoni
Rewrote the template of ComputeNormal functions to a more readable form.

Revision 1.39  2006/07/12 12:14:31  zifnab1974
changes for compilation on linux. Not sure about using either SphereOfTriangle or SphereOfTetra, please check.

Revision 1.38  2006/07/06 12:45:08  ganovelli
added SmallestEnclosingSphere

Revision 1.37  2006/01/22 10:00:43  cignoni
Very Important Change: Area->DoubleArea (and no more Area function)

Revision 1.36  2005/12/16 11:42:23  corsini
Add some user bit

Revision 1.35  2005/12/01 23:54:29  cignoni
Added HasFlags

Revision 1.34  2005/11/23 14:40:09  pietroni
added cFFi function

Revision 1.33  2005/11/22 15:47:35  cignoni
Moved ComputeNormal and ComputeNormalizedNormal out of the face class (no more a member function!)

Revision 1.32  2005/11/12 18:39:54  cignoni
Added dummy static member for avoiding annoying warning in empty functions...

Revision 1.31  2005/11/01 18:16:36  cignoni
Added intialization of _flags to zero in the default constructor of face

Revision 1.30  2005/10/13 09:25:43  cignoni
Added cFFp and cVFp const member functions

Revision 1.29  2005/09/28 19:32:09  m_di_benedetto
Added const qualifier in GetBBox method.

Revision 1.28  2005/06/17 00:43:34  cignoni
Added a named typedef for the per wedge TexCoordinate

Revision 1.27  2005/03/18 16:35:53  fiorin
minor changes to comply gcc compiler

Revision 1.26  2005/03/11 14:14:14  ganovelli
_ffi was a 4 for positions vector (only 3  used)

Revision 1.25  2005/01/28 17:53:13  pietroni
added HasEdgePlane function

Revision 1.24  2004/10/28 00:50:48  cignoni
Better Doxygen documentation

Revision 1.23  2004/10/25 08:22:40  ganovelli
IsBOrder (typecast on return type)

Revision 1.22  2004/10/20 08:28:31  fiorin
Added constant access function FFp and renamed F1 F2 to FFp1 FFp2

Revision 1.21  2004/10/18 17:13:50  ganovelli
added  ::IsBorder

Revision 1.20  2004/09/15 11:20:15  ganovelli
changed P() to cP()

Revision 1.19  2004/09/14 19:47:02  ganovelli
removed "&" in FFp

Revision 1.18  2004/08/25 15:15:27  ganovelli
minor changes to comply gcc compiler (typename's and stuff)

Revision 1.17  2004/07/15 12:03:07  ganovelli
minor changes

Revision 1.16  2004/07/15 11:31:59  ganovelli
minor changes

Revision 1.15  2004/07/12 12:17:09  pietroni
added function NormalizedNormal

Revision 1.14  2004/05/13 11:01:06  turini
Changed ComputeMormalizedNormal() using Triangle3

Revision 1.13  2004/05/12 18:49:05  ganovelli
dist and coputeRT removed (see distance.h and updateEdges)

Revision 1.12  2004/05/12 14:43:36  cignoni
removed warning of unused variables

Revision 1.11  2004/05/12 12:50:20  turini
include color4

Revision 1.10  2004/05/10 14:01:09  ganovelli
assert(i*0) for using "i" and preventing the compiler warning for unreferenced variable

Revision 1.9  2004/05/10 13:19:38  cignoni
Added mandatory template params for edge and face class names to the face class
Changed type of return face pointer to the one passed by templ params
Changed name of func FV to VF (it stores Vertex-Face Topology)

Revision 1.8  2004/05/06 09:06:59  pietroni
changed names to topology functions

Revision 1.7  2004/05/04 02:46:23  ganovelli
added function Dist

Revision 1.5  2004/04/05 11:51:22  cignoni
wrong define FACE_N instead of FACE_FN

Revision 1.4  2004/03/29 08:37:09  cignoni
missing include

Revision 1.3  2004/03/10 00:52:38  cignoni
Moved geometric stuff to the space/triangle  class

Revision 1.2  2004/03/03 16:08:38  cignoni
First working version

Revision 1.1  2004/02/13 00:44:45  cignoni
First commit...

****************************************************************************/

#pragma message("[VCGLIB Warning]  this way to define the simplex face is DEPRECATED and no more SUPPORTED") 
#pragma message("[VCGLIB Warning]  use vcg/simplex/faceplus instead ") 

#ifndef FACE_TYPE 
#pragma error message("\nYou should never directly include this file\n")
#else


#include <vcg/math/base.h>
#include <vcg/space/box3.h>
#include <vcg/space/texcoord2.h>
#include <vcg/space/triangle3.h>
#include <vcg/space/sphere3.h>
#include <vcg/space/color4.h>
#include <vcg/space/plane3.h>
#include <vcg/space/smallest_enclosing.h>
#include <vcg/simplex/face/topology.h>

namespace vcg {
class DUMMYEDGETYPE;
class DUMMYFACETYPE;
class DUMMYTETRATYPE;

/** \addtogroup face */
//@{
/*!
 * This class represent the generic configurable Face; 
 * Usually you never direclty use this class with this name but you build 
 * your own type by directly including one of the .h files under the face/with 
 * directory. Each file specify a class type with the desired fields. So for example 
 * including 'vcg/simplex/face/with/FCFN.h' allow you to use the class FaceFCFN that has per-face color and normal stored inside.
 */
template <class FVTYPE, class FETYPE, class FFTYPE, class TCTYPE = TexCoord2<float,1> > class FACE_TYPE
{
public:
	///	The base type of the face
	typedef FACE_TYPE BaseFaceType;
	///	The base type of the face itself
	typedef FFTYPE FaceType;
	/// The vertex type
	typedef FVTYPE VertexType;
	/// The vertex type
	typedef FETYPE EdgeType;
	/// The Texture Coordinate type
	typedef TCTYPE TexCoordType;
	/// The type of the scalar field of the vertex coordinate
  typedef typename VertexType::ScalarType ScalarType;
	/// The type of the the vertex coordinate
	typedef Point3< ScalarType > CoordType;
	typedef Point3< ScalarType > NormalType;
	/// The geometric type of the face
	typedef Triangle3<ScalarType> GeometricType;
	/// The type of the quality (same as scalar)
	typedef  ScalarType  QualityType;

  typedef typename FVTYPE::FaceType FaceFromVertType;
	/// The bounding box type
	typedef Box3<ScalarType> BoxType;
	
  /// Default Empty Costructor
  inline FACE_TYPE(){_flags=0;}

	/// This are the _flags of face, the default value is 0
	int  _flags;		

/***********************************************/
/** @name Vertex Pointer
Functions to access to the vertexes of the face;
**/
  //@{
protected:
	/// Vector of vertex pointer incident in the face
	VertexType *v[3];
public:
	/** Return the pointer to the j-th vertex of the face.
		@param j Index of the face vertex.
	 */
	inline VertexType * & V( const int j )
	{	
		assert( (_flags & DELETED) == 0 );
		assert( (_flags & NOTREAD) == 0 ); 
		assert( (_flags & NOTWRITE) == 0 );
		assert(j >= 0);
		assert(j <  3);
		return v[j];
	}

	inline  VertexType * const &  V( const int j ) const
	{
		assert( (_flags & DELETED) == 0 );
		assert( (_flags & NOTREAD) == 0 );
		assert(j>=0);
		assert(j<3);
		return v[j];
	}
	inline  VertexType * const  cV( const int j ) const
	{
		assert( (_flags & DELETED) == 0 );
		assert( (_flags & NOTREAD) == 0 );
		assert(j>=0);
		assert(j<3);
		return v[j];
	}

	// Shortcut per accedere ai punti delle facce
	inline CoordType & P( const int j )
	{	
		assert( (_flags & DELETED) == 0 );
		assert( (_flags & NOTREAD) == 0 ); 
		assert( (_flags & NOTWRITE) == 0 );
		assert(j>=0);
		assert(j<3);
		return v[j]->P();
	}

	inline const CoordType & P( const int j ) const
	{
		assert( (_flags & DELETED) == 0 );
		assert( (_flags & NOTREAD) == 0 );
		assert(j>=0);
		assert(j<3);
		return v[j]->cP();
	}
	inline const CoordType & cP( const int j ) const
	{
		assert( (_flags & DELETED) == 0 );
		assert( (_flags & NOTREAD) == 0 );
		assert(j>=0);
		assert(j<3);
		return v[j]->cP();
	}

	/** Return the pointer to the ((j+1)%3)-th vertex of the face.
		@param j Index of the face vertex.
	 */
	inline VertexType * & V0( const int j ) { return V(j);}
	inline VertexType * & V1( const int j ) { return V((j+1)%3);}
	inline VertexType * & V2( const int j ) { return V((j+2)%3);}
	inline const VertexType * const &  V0( const int j ) const { return V(j);}
	inline const VertexType * const &  V1( const int j ) const { return V((j+1)%3);}
	inline const VertexType * const &  V2( const int j ) const { return V((j+2)%3);}
	inline const VertexType * const & cV0( const int j ) const { return cV(j);}
	inline const VertexType * const & cV1( const int j ) const { return cV((j+1)%3);}
	inline const VertexType * const & cV2( const int j ) const { return cV((j+2)%3);}

	/// Shortcut per accedere ai punti delle facce
	inline CoordType & P0( const int j ) { return V(j)->P();}
	inline CoordType & P1( const int j ) { return V((j+1)%3)->P();}
	inline CoordType & P2( const int j ) { return V((j+2)%3)->P();}
	inline const CoordType &  P0( const int j ) const { return V(j)->P();}
	inline const CoordType &  P1( const int j ) const { return V((j+1)%3)->P();}
	inline const CoordType &  P2( const int j ) const { return V((j+2)%3)->P();}
	inline const CoordType & cP0( const int j ) const { return cV(j)->P();}
	inline const CoordType & cP1( const int j ) const { return cV((j+1)%3)->P();}
	inline const CoordType & cP2( const int j ) const { return cV((j+2)%3)->P();}

	inline VertexType * & UberV( const int j )
	{	
		assert(j>=0);
		assert(j<3);
		return v[j];
	}

	inline const VertexType * const & UberV( const int j ) const
	{
		assert(j>=0);
		assert(j<3);
		return v[j];
	}


  //@}

/***********************************************/
/** @name Normal
    blah
    blah
**/
  //@{

#ifdef __VCGLIB_FACE_FN
	/// This vector indicates the normal of the face (defines if FACE_N is defined)
protected:
	CoordType _n;
public:
#endif

  /// Return the reference of the normal to the face (if __VCGLIB_FACE_FN is defined).
	inline CoordType & N()
	{
#ifdef __VCGLIB_FACE_FN
	return _n;
#else
	assert(0);
	return *(CoordType *)0;
#endif
	}
		/// Return the reference of the normal to the face (if __VCGLIB_FACE_FN is defined).
	inline const CoordType & N() const
	{
#ifdef __VCGLIB_FACE_FN
		return _n;
#else
	return *(CoordType *)0;
#endif
	}
	/// Return the reference of the normal to the face (if __VCGLIB_FACE_FN is defined).
	inline const CoordType cN() const
	{
#ifdef __VCGLIB_FACE_FN
		return _n;
#else
	return *(CoordType *)0;
#endif
	}

/// Return the value of the face normal as it correspond to the current geometry.
/// it is always computed and never stored. 
const CoordType Normal() const
{
	return vcg::Normal(*this);
}

/// Return the value of the face normal as it correspond to the current geometry.
/// it is always computed and never stored. 
const CoordType NormalizedNormal() const
{
	return vcg::NormalizedNormal(*this);
}

#ifdef __VCGLIB_FACE_WN
	/// This vector indicates per wedge normal 
	CoordType _wn[3];
#endif

public:
	CoordType & WN(const int i)
	{
#ifdef __VCGLIB_FACE_WN
		return _wn[i];
#else
		assert(0);
		return *(CoordType *)(&_flags);
#endif
	}

const CoordType & WN(const int i) const
	{
#ifdef __VCGLIB_FACE_WN
		return _wn[i];
#else
		return CoordType();
#endif
	}

  //@}

/***********************************************/
/** @name Quality
    blah
    blah
**/
  //@{

#ifdef __VCGLIB_FACE_FQ
protected:
	float _q;
#endif
public:
	float & Q()
	{
#ifdef __VCGLIB_FACE_FQ
		return _q;
#else
		assert(0);
		return *(float*)(&_flags);
#endif
	}

const float & Q() const
	{
#ifdef __VCGLIB_FACE_FQ
		return _q;
#else
		assert(0);
		return *(float*)(&_flags);
#endif
	}

  //@}

/***********************************************/
/** @name Texture
    blah
    blah
**/
  //@{

// Per Wedge Texture Coords
protected:
#ifdef __VCGLIB_FACE_WT
	TCTYPE _wt[3];
#endif
public:
	TCTYPE & WT(const int i)
	{
#ifdef __VCGLIB_FACE_WT
		return _wt[i];
#else
		assert(0);
		return *(TCTYPE*)(&_flags +i) ;
#endif
	}

	const TCTYPE & WT(const int i) const
	{
#ifdef __VCGLIB_FACE_WT
		return _wt[i];
#else
		assert(0);
		return *(TCTYPE*)(&_flags);
#endif
	}


 //@}

/***********************************************/
/** @name Colors
    blah
    blah
**/
  //@{
protected:
#ifdef __VCGLIB_FACE_FC
	Color4b _c;
#endif

public:
	Color4b & C()
	{
#ifdef __VCGLIB_FACE_FC
		return _c;
#else
		assert(0);
		return *(Color4b*)(&_flags);
#endif
	}

	const Color4b C() const
	{
#ifdef __VCGLIB_FACE_FC
		return _c;
#else
		return Color4b(Color4b::White);
#endif
	}

protected:
#ifdef __VCGLIB_FACE_WC
	Color4b _wc[3];
#endif
public:
	Color4b & WC(const int i)
	{
#ifdef __VCGLIB_FACE_WC
		return _wc[i];
#else
		assert(0);
		return *(Color4b*)(&_flags + i);
#endif
	}

const Color4b WC(const int i) const
	{
#ifdef __VCGLIB_FACE_WC
		return _wc[i];
#else
		assert(0);
		return Color4b(Color4b::White);
#endif
	}




  //@}

/***********************************************/
/** @name Adjacency
    blah
    blah
**/
  //@{

#if (defined(__VCGLIB_FACE_AF) && defined(__VCGLIB_FACE_AS))
	#error Error: You cannot specify face-to-face and shared topology together
#endif

#if (defined(__VCGLIB_FACE_AV) && defined(__VCGLIB_FACE_AS))
	#error Error: You cannot specify vertex-face and shared topology together
#endif

protected:
#if defined(__VCGLIB_FACE_AF)
  /// Vector of face pointer, it's used to indicate the adjacency relations (defines if FACE_A is defined)
	FFTYPE   *_ffp[3];				// Facce adiacenti
	/// Index of the face in the arrival face 
	char _ffi[3];									
#endif

#ifdef __VCGLIB_FACE_AV
	///Vettore di puntatori a faccia, utilizzato per indicare le adiacenze vertice faccia
	FFTYPE *_fvp[3];
	char _fvi[3];
#endif

#ifdef __VCGLIB_FACE_AS
	///Vettore di puntatori a faccia, utilizzato per indicare le adiacenze vertice faccia
	FFTYPE *fs[3];
	char zs[3];
#endif
public:




	/** Return the pointer to the j-th adjacent face.
	    @param j Index of the edge.
	 */
	inline FFTYPE * & FFp( const int j )
	{
		assert( (_flags & DELETED) == 0 );
		assert( (_flags & NOTREAD) == 0 );
		assert( (_flags & NOTWRITE) == 0 );
		assert(j>=0);
	    assert(j<3);
#if defined(__VCGLIB_FACE_AF)
		  return _ffp[j];
#elif defined(__VCGLIB_FACE_AS)
			return fs[j];
#else 
		assert(0);
    static FFTYPE *dum=0; dum+=j;
		return dum;

#endif
	}

	inline const FFTYPE * FFp( const int j ) const 
	{
		assert( (_flags & DELETED) == 0 );
		assert( (_flags & NOTREAD) == 0 );
		assert(j>=0 && j<3);

#if defined(__VCGLIB_FACE_AF)
		  return _ffp[j];
#elif defined(__VCGLIB_FACE_AS)
			return fs[j];
#else 
		assert(0);
    static FFTYPE *dum=0; dum+=j;
		return dum;
#endif
	}
  inline const FFTYPE * cFFp( const int j ) const {return FFp(j);}
	
	inline FFTYPE * & FFp1( const int j ) { return FFp((j+1)%3);}
	inline FFTYPE * & FFp2( const int j ) { return FFp((j+2)%3);}
	inline const FFTYPE * const&  FFp1( const int j ) const { return FFp((j+1)%3);}
	inline const FFTYPE * const&  FFp2( const int j ) const { return FFp((j+2)%3);}



/** Return the pointer to the j-th adjacent face.
	    @param j Index of the edge.
	 */
	inline FFTYPE * & UberF( const int j )
	{
		assert(j>=0);
	  assert(j<3);
#if defined(__VCGLIB_FACE_AF)
		  return _ffp[j];
#elif defined(__VCGLIB_FACE_AS)
			return fs[j];
#else 
		assert(0); // if you stop here you are probably trying to use FF topology in a face without it
		return *((FFTYPE **)(_flags));
#endif
	}

	inline const FFTYPE * const & UberF( const int j ) const
	{
		assert(j>=0);
	  assert(j<3);
#if defined(__VCGLIB_FACE_AF)
		  return _ffp[j];
#elif defined(__VCGLIB_FACE_AS)
			return fs[j];
#else
		assert(0); // if you stop here you are probably trying to use FF topology in a face without it
		return *((FFTYPE **)(_flags));
#endif
	}
	

	inline FFTYPE * & VFp( const int j )
	{
		assert( (_flags & DELETED) == 0 );
		assert( (_flags & NOTREAD) == 0 );
		assert( (_flags & NOTWRITE) == 0 );
		assert(j>=0);
		assert(j<3);
#ifdef __VCGLIB_FACE_AV
		return _fvp[j];
#elif defined(__VCGLIB_FACE_AS)
		return fs[j];
#else
		assert(0); // you are probably trying to use VF topology in a vertex without it
		return *((FFTYPE **)(_flags));
#endif
	}

	inline const FFTYPE * const VFp( const int j ) const
	{
		assert( (_flags & DELETED) == 0 );
		assert( (_flags & NOTREAD) == 0 );
		assert(j>=0);
		assert(j<3);
#ifdef __VCGLIB_FACE_AV
		return _fvp[j];
#elif defined(__VCGLIB_FACE_AS)
		return fs[j];
#else
		assert(0);
    static FFTYPE * const DummyVal=0;
		return DummyVal;
#endif
	}
  inline const FFTYPE * cVFp( const int j ) const {return VFp(j);}


	/** Return the index that the face have in the j-th adjacent face.
	    @param j Index of the edge.
	 */
	inline char & FFi( const int j )
	{
		assert( (_flags & DELETED) == 0 );
		assert( (_flags & NOTREAD) == 0 );
		assert( (_flags & NOTWRITE) == 0 );
		assert(j>=0);
		assert(j<3);
#if defined(__VCGLIB_FACE_AF) 
		return _ffi[j];
#elif defined(__VCGLIB_FACE_AS) 
		return zs[j];
#else
		assert(0);
		return *(char *)&_flags; // tanto per farlo compilare...
#endif
	}

	inline const char & FFi( const int j ) const
	{
		assert( (_flags & DELETED) == 0 );
		assert( (_flags & NOTREAD) == 0 );
		assert(j>=0);
		assert(j<3);
#if defined(__VCGLIB_FACE_AF) 
		return _ffi[j];
#elif defined(__VCGLIB_FACE_AS) 
		return zs[j];
#else
		assert(0);
		return *(char *)&_flags;
#endif
	}

	inline const char & cFFi( const int j ) const {return FFi(j);}

		/** Return the index that the face have in the j-th adjacent face.
	    @param j Index of the edge.
	 */
	inline char & UberZ( const int j )
	{
		assert(j>=0);
		assert(j<3);
#if defined(__VCGLIB_FACE_AF) 
		return _ffi[j];
#elif defined(__VCGLIB_FACE_AS) 
		return zs[j];
#else
		assert(0);
		return *(char *)&_flags;
#endif
	}

	inline const char & UberZ( const int j ) const
	{
		assert(j>=0);
		assert(j<3);
#if defined(__VCGLIB_FACE_AF) 
		return _ffi[j];
#elif defined(__VCGLIB_FACE_AS) 
		return zs[j];
#else
		assert(0);
		return *(char *)&_flags;
#endif
	}


	inline char & VFi( const int j )
	{
		assert( (_flags & DELETED) == 0 );
		assert( (_flags & NOTREAD) == 0 );
		assert( (_flags & NOTWRITE) == 0 );
		assert(j>=0);
		assert(j<3);
#ifdef __VCGLIB_FACE_AV
		return _fvi[j];
#elif defined(__VCGLIB_FACE_AS)
		return zs[j];
#else
		assert(0);
		return *(char *)&_flags;
#endif
	}

	inline const char & VFi( const int j ) const
	{
		assert( (_flags & DELETED) == 0 );
		assert( (_flags & NOTREAD) == 0 );
		assert(j>=0);
		assert(j<3);
#ifdef __VCGLIB_FACE_AV
		return _fvi[j];
#elif defined(__VCGLIB_FACE_AS)
		return zs[j];
#else
		assert(0);
		return *(char *)&_flags;
#endif
	}

  //@}

/***********************************************/
/** @name Mark
    blah
    blah
**/
  //@{


#ifdef __VCGLIB_FACE_FM
	/// Incremental mark (defines if FACE_I is defined)
	int imark;
	inline int & IMark()
	{
		assert( (_flags & DELETED) == 0 );
		assert( (_flags & NOTREAD) == 0 );
		assert( (_flags & NOTWRITE) == 0 );
		return imark;
	}

	inline const int & IMark() const
	{
		assert( (_flags & DELETED) == 0 );
		assert( (_flags & NOTREAD) == 0 );
		return imark;
	}
#endif // Mark

	/// Initialize the imark system of the face
	inline void InitIMark()
	{
#ifdef __VCGLIB_FACE_FM
		imark = 0;
#endif
	}

  
  //@}
/***********************************************/
/** @name Flags
    blah
    blah
**/
  //@{


	enum {
		// This bit indicate that the face is deleted from the mesh
		DELETED     = 0x00000001,		// cancellato
		// This bit indicate that the face of the mesh is not readable
		NOTREAD     = 0x00000002,		// non leggibile (ma forse modificabile)
		// This bit indicate that the face is not modifiable
		NOTWRITE    = 0x00000004,		// non modificabile (ma forse leggibile) 
		// This bit indicate that the face is modified
		SELECTED    = 0x00000020,		// Selection _flags
		// Border _flags, it is assumed that BORDERi = BORDER0<<i 
		BORDER0     = 0x00000040,
		BORDER1     = 0x00000080,
		BORDER2     = 0x00000100,
		// Face Orientation Flags, used efficiently compute point face distance  
		NORMX		= 0x00000200,
		NORMY		= 0x00000400,
		NORMZ		= 0x00000800,
		// Crease _flags,  it is assumed that FEATUREi = FEATURE0<<i 
		FEATURE0    = 0x00008000,
		FEATURE1    = 0x00010000,
		FEATURE2    = 0x00020000,
		// User bits
		USER0       = 0x00040000,
		USER1       = 0x00080000,
		USER2       = 0x00100000,
		USER3       = 0x00200000
		};
public:
	static int &LastBitFlag()
		{
			static int b =USER0;
			return b;
		}
	static inline int NewBitFlag()
		{
			LastBitFlag()=LastBitFlag()<<1;
			return LastBitFlag();
		}
	static inline bool DeleteBitFlag(int bitval)
		{	
			if(LastBitFlag()==bitval) {
					LastBitFlag()= LastBitFlag()>>1;
					return true;
			}
			assert(0);
			return false;
		}

    void ClearFlags() {_flags=0;}

	/// Return the _flags.
	inline int & Flags ()
	{
		assert( (_flags & DELETED) == 0 );
		assert( (_flags & NOTREAD) == 0 );
		return _flags;
	}

	inline const int & Flags () const
	{
		assert( (_flags & DELETED) == 0 );
		assert( (_flags & NOTREAD) == 0 );
		return _flags;
	}
	/// Ritorna il _flags senza effettuare alcun controllo sui relativi bit
	inline int & UberFlags()
	{
		return _flags;
	}

	inline const int UberFlags() const
	{
		return _flags;
	}

	/// This function checks if the face is deleted 
	bool IsD() const {return (_flags & DELETED) != 0;}
	/// This function mark the face as deleted
	void SetD()		{_flags |=DELETED;}
	/// This function mark the face as not deleted
	void ClearD()	{_flags &= (~DELETED);}
	/// This function checks if the face is deleted 
	bool IsDeleted() const {return IsD();}
	
	/// This function checks if the face is readable 
	bool IsR() const {return (_flags & NOTREAD) == 0;}
	/// This function marks the face as readable
	void SetR()		{_flags &= (~NOTREAD);}
	/// This function marks the face as not readable
	void ClearR() {_flags |=NOTREAD;}
	
	/// This function checks if the face is readable 
	bool IsW() const {return (_flags & NOTWRITE)== 0;}
	/// This function marks the vertex as not writable
	void SetW() {_flags &=(~NOTWRITE);}
	/// This function marks the face as not writable
	void ClearW() {_flags |=NOTWRITE;}

	/// This funcion checks whether the face is both readable and modifiable
	bool IsRW() const {return (_flags & (NOTREAD | NOTWRITE)) == 0;}

	
	/// This function checks if the face is selected
	bool IsS() const {return (_flags & SELECTED) != 0;}
	/// This function select the face
	void SetS()		{_flags |=SELECTED;}
	/// This funcion execute the inverse operation of SetS()
	void ClearS()	{_flags &= (~SELECTED);}

	/// This function checks if the face is selected
	bool IsB(int i) const {return (_flags & (BORDER0<<i)) != 0;}
	/// This function select the face
	void SetB(int i)		{_flags |=(BORDER0<<i);}
	/// This funcion execute the inverse operation of SetS()
	void ClearB(int i)	{_flags &= (~(BORDER0<<i));}

	/// This function checks if the face is Crease  on side i
	bool IsF(int i) const {return (_flags & (FEATURE0<<i)) != 0;}
	/// This function select the face flag
	void SetF(int i)		{_flags |=(FEATURE0<<i);}
	/// This funcion execute the inverse operation of Set()
	void ClearF(int i)	{_flags &= (~(FEATURE0<<i));}

	/// This function checks if the given user bit is true
	bool IsUserBit(int userBit){return (_flags & userBit) != 0;}
	/// This function set  the given user bit 
	void SetUserBit(int userBit){_flags |=userBit;}
	/// This function clear the given user bit 
	void ClearUserBit(int userBit){_flags &= (~userBit);}


//@}
/*#*******************	
*  Bounding box *
**********************/

void GetBBox( BoxType & bb ) const
{
	bb.Set( v[0]->P() );
	bb.Add( v[1]->P() );
	bb.Add( v[2]->P() );
}

 /***********************************************/
 /** @name Reflection Functions 
 Static functions  that give information about the current vertex type.
Reflection is a mechanism making it possible to investigate yourself. Reflection is used to investigate format of objects at runtime, invoke methods and access fields of these objects. Here we provide static const functions that are resolved at compile time and they give information about the data (normal, color etc.) supported by the current vertex type.
 **/
 //@{
static bool HasFlags()  { 
  return true;
}
static bool HasFaceNormal()  { 
#ifdef __VCGLIB_FACE_FN 
  return true;
#else
  return false;
#endif
}
static bool HasFaceQuality()  { 
#ifdef __VCGLIB_FACE_FQ
  return true;
#else
  return false;
#endif
}
static bool HasFaceColor()  { 
#ifdef __VCGLIB_FACE_FC 
  return true;
#else
  return false;
#endif
}
static bool HasEdgePlane()  { 
#ifdef __VCGLIB_FACE_RT 
  return true;
#else
  return false;
#endif
}
static bool HasFFAdjacency()  { 
#if (defined(__VCGLIB_FACE_AF) || defined(__VCGLIB_FACE_AS))
  return true;
#else
  return false;
#endif
}
static bool HasVFAdjacency()  { 
#if (defined(__VCGLIB_FACE_AV) || defined(__VCGLIB_FACE_AS))
  return true;
#else
  return false;
#endif
}
static bool HasSharedAdjacency()  { 
#if defined(__VCGLIB_FACE_AS)
  return true;
#else
  return false;
#endif
}
static bool HasFaceMark()  { 
#ifdef __VCGLIB_FACE_FC 
  return true;
#else
  return false;
#endif
}
static bool HasWedgeColor()  { 
#ifdef __VCGLIB_FACE_WC 
  return true;
#else
  return false;
#endif
}
static bool HasWedgeTexCoord()  { 
#ifdef __VCGLIB_FACE_WT 
  return true;
#else
  return false;
#endif
}
static bool HasWedgeNormal()  { 
#ifdef __VCGLIB_FACE_WN 
  return true;
#else
  return false;
#endif
}

//@}

  /// operator to compare two faces
	inline bool operator == ( const FFTYPE & f ) const {
		for(int i=0; i<3; ++i)
			if( (V(i) != f.V(0)) && (V(i) != f.V(1)) && (V(i) != f.V(2)) )
				return false;
		return true;
	}

/** Calcola i coefficienti della combinazione convessa.
	@param bq Punto appartenente alla faccia
	@param a Valore di ritorno per il vertice V(0)
	@param b Valore di ritorno per il vertice V(1)
	@param _c Valore di ritorno per il vertice V(2)
	@return true se bq appartiene alla faccia, false altrimenti
*/
bool InterpolationParameters(const CoordType & bq, ScalarType &a, ScalarType &b, ScalarType &_c ) const
{	
const ScalarType EPSILON = ScalarType(0.000001);


#define x1 (cV(0)->P()[0])
#define y1 (cV(0)->P()[1])
#define z1 (cV(0)->P()[2])
#define x2 (cV(1)->P()[0])
#define y2 (cV(1)->P()[1])
#define z2 (cV(1)->P()[2])
#define x3 (cV(2)->P()[0])
#define y3 (cV(2)->P()[1])
#define z3 (cV(2)->P()[2])
#define px (bq[0])
#define py (bq[1])
#define pz (bq[2])

     ScalarType t1  = px*y2;
     ScalarType t2  = px*y3;
     ScalarType t3  = py*x2;
     ScalarType t4  = py*x3;
     ScalarType t5  = x2*y3;
     ScalarType t6  = x3*y2;
     ScalarType t8  = x1*y2;
     ScalarType t9  = x1*y3;
     ScalarType t10 = y1*x2;
     ScalarType t11 = y1*x3;
     ScalarType t13 = t8-t9-t10+t11+t5-t6;
     if(fabs(t13)>=EPSILON)
	 {
         ScalarType t15 = px*y1;
         ScalarType t16 = py*x1;
         a =  (t1 -t2-t3 +t4+t5-t6 )/t13;
         b = -(t15-t2-t16+t4+t9-t11)/t13;
         _c =  (t15-t1-t16+t3+t8-t10)/t13;
		return true;
     }

     t1  = px*z2;
     t2  = px*z3;
     t3  = pz*x2;
     t4  = pz*x3;
     t5  = x2*z3;
     t6  = x3*z2;
     t8  = x1*z2;
     t9  = x1*z3;
     t10 = z1*x2;
     t11 = z1*x3;
     t13 = t8-t9-t10+t11+t5-t6;
     if(fabs(t13)>=EPSILON)
	 {
		ScalarType t15 = px*z1;
		ScalarType t16 = pz*x1;
		a =  (t1 -t2-t3 +t4+t5-t6 )/t13;
		b = -(t15-t2-t16+t4+t9-t11)/t13;
		_c =  (t15-t1-t16+t3+t8-t10)/t13;
		return true;
     }

     t1  = pz*y2; t2  = pz*y3;
     t3  = py*z2; t4  = py*z3;
     t5  = z2*y3; t6  = z3*y2;
     t8  = z1*y2; t9  = z1*y3;
     t10 = y1*z2; t11 = y1*z3;
     t13 = t8-t9-t10+t11+t5-t6;
     if(fabs(t13)>=EPSILON)
	 {
         ScalarType t15 = pz*y1;
         ScalarType t16 = py*z1;
         a =  (t1 -t2-t3 +t4+t5-t6 )/t13;
         b = -(t15-t2-t16+t4+t9-t11)/t13;
         _c =  (t15-t1-t16+t3+t8-t10)/t13;
		return true;
     }
	 
#undef x1
#undef y1
#undef z1
#undef x2
#undef y2
#undef z2
#undef x3
#undef y3
#undef z3
#undef px
#undef py
#undef pz

     return false;
}



/// Return the DOUBLE of the area of the face
// NOTE the old Area function has been removed to intentionally 
// cause compiling error that will help people to check their code...
// A some  people used Area assuming that it returns the double and some not. 
// So please check your codes!!!
// And please DO NOT Insert any Area named function here!
ScalarType DoubleArea() const
{
	return ( (V(1)->cP() - V(0)->cP()) ^ (V(2)->cP() - V(0)->P()) ).Norm();
}

CoordType Barycenter() const
{
	return (V(0)->P()+V(1)->P()+V(2)->P())/ScalarType(3.0);
}

Sphere3<ScalarType> SmallestEnclosingSphere() const
{
		return SmallestEnclosing::SphereOfTriangle<ScalarType>(*this);
}

ScalarType Perimeter() const
{
	return Distance(V(0)->P(),V(1)->P())+
		     Distance(V(1)->P(),V(2)->P())+
				 Distance(V(2)->P(),V(0)->P());
}

/// Return the _q of the face, the return value is in [0,sqrt(3)/2] = [0 - 0.866.. ]
ScalarType QualityFace( ) const
{
	
	return Quality(V(0)->cP(), V(1)->cP(), V(2)->cP());
	/*
	CoordType d10 = V(1)->P() - V(0)->P();
	CoordType d20 = V(2)->P() - V(0)->P();
	CoordType d12 = V(1)->P() - V(2)->P();

	CoordType x = d10^d20;

	ScalarType a = Norm( x );		// doppio dell' Area
	ScalarType b;
	
	b = Norm2( d10 );
	ScalarType t = b; 
	t = Norm2( d20 ); if( b<t ) b = t;
	t = Norm2( d12 ); if( b<t ) b = t;

	assert(b!=0.0);

	return a/b;*/

}

// Funzione di supporto
inline void Nexts( BaseFaceType *&f,int &z )
{
    int t;
    t = z;
    z = (*f).Z(z);
    f = (*f).F(t);
}

	// Sezione dist e ray
#ifdef __VCGLIB_FACE_RT

	CoordType & Edge(const int& i){return edges[i];};
	CoordType  cEdge(const int& i)const{return  edges[i];};

	Plane3<ScalarType>  & Plane(){return plane;};
	Plane3<ScalarType>   cPlane()const{return plane;};


	CoordType edges[3];;
	Plane3<ScalarType> plane;
#endif

	/// return the index [0..2] of a vertex in a face
	inline int VertexIndex( const VertexType * w ) const
	{
			 if( v[0]==w ) return  0;
		else if( v[1]==w ) return  1;
		else if( v[2]==w ) return  2;
		else               return -1;
	}


}; //end Class


#ifndef __VCG_FACE_BASE_SINGLE
#define __VCG_FACE_BASE_SINGLE
// Note that while the whole file can be included more than once, this portion of the file MUST be included once

/// Calculate the normal to the face, the value is store in the field _n of the face
namespace face
{

template <class FaceType>
void ComputeNormal(FaceType &f) {	f.N() = vcg::Normal< FaceType >(f); }

template <class FaceType>
void ComputeNormalizedNormal(FaceType &f) {	f.N() = vcg::NormalizedNormal< FaceType >(f); }

} // end namespace face
//@}
#endif

}	 // end namespace vcg


#endif

