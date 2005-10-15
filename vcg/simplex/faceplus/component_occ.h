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

****************************************************************************/

/* 
Note
OCC = Optional Component Compact
compare with OCF(Optional Component Fast)

****************************************************************************/
#ifndef __VCG_FACE_PLUS_COMPONENT_OCC
#define __VCG_FACE_PLUS_COMPONENT_OCC

#include <vcg/simplex/faceplus/component.h>
#include <vcg/container/vector_occ.h>


namespace vcg {
  namespace face {

///*-------------------------- NORMAL ----------------------------------------*/ 

	template <class A, class T> class NormalOcc: public T {
	public:
		typedef A NormalType;
		NormalType &N() {return CAT< vector_occ<FaceType>,NormalType>::Instance()->Get((FaceType*)this);}
	  static bool HasNormal()   { return true; }
		static bool HasColorOcc()   { return true; }
	};

	template <class T> class Normal3sOcc: public NormalOcc<vcg::Point3s, T> {};
	template <class T> class Normal3fOcc: public NormalOcc<vcg::Point3f, T> {};
	template <class T> class Normal3dOcc: public NormalOcc<vcg::Point3d, T> {};

///*-------------------------- COLOR ----------------------------------------*/ 

template <class A, class T> class ColorOcc: public T {
public:
  typedef A ColorType;
  ColorType &C() { return CAT< vector_occ<FaceType>,ColorType>::Instance()->Get((FaceType*)this); }
  static bool HasColor()   { return true; }
  static bool HasColorOcc()   { return true; }
};

template <class T> class Color4bOcc: public ColorOcc<vcg::Color4b, T> {};

/*----------------------------- VFADJ ---------------------------------------*/ 

// questo tipo serve per tenere tutte le informazioni sull'adiacenza dentro una
// singola classe
template <class FP>
struct VFAdjTypeSup {
		FP  _vfp[3]; 
		char _vfi[3];
	};

template <class A, class T> class VFAdjOccBase: public T {
public:
  typename T::FacePointer &VFp(const int j) {
		return (CAT< vector_occ<FaceType>,VFAdjTypeSup<T::FacePointer> >::Instance()->Get((FaceType*)this))._vfp[j];}

  typename T::FacePointer cVFp(const int j) const {
		return (CAT< vector_occ<FaceType>,VFAdjTypeSup<T::FacePointer> >::Instance()->Get((FaceType*)this))._vfp[j];}

  char &VFi(const int j) { return (CAT< vector_occ<FaceType>,VFAdjTypeSup<T::FacePointer> >::Instance()->Get((FaceType*)this))._vfi[j];}

  static bool HasVFAdjacency()   {   return true; }
  static bool HasVFAdjacencyOcc()   { return true; }
};

template <class T> class VFAdjOcc : public VFAdjOccBase<VFAdjTypeSup<typename T::FacePointer>,T>{};

/*----------------------------- FFADJ -----------------------------------*/ 

// questo tipo serve per tenere tutte le informazioni sull'adiacenza dentro una
// singola classe
template <class FP>
struct FFAdjTypeSup {
		FP  _ffp[3]; 
		char _ffi[3];
	};

template <class A, class T> class FFAdjOccBase: public T {
public:
	
	typedef  A FFAdjType;

  typename T::FacePointer &FFp(const int j) { 
		return (CAT< vector_occ<FaceType>,FFAdjTypeSup<T::FacePointer> >::Instance()->Get((FaceType*)this))._ffp[j];}

  typename T::FacePointer const  FFp(const int j) const { 
		return (CAT< vector_occ<FaceType>,FFAdjTypeSup<T::FacePointer> >::Instance()->Get((FaceType*)this))._ffp[j];}

  typename T::FacePointer const cFFp(const int j) const {
    return (CAT< vector_occ<FaceType>,FFAdjTypeSup<T::FacePointer> >::Instance()->Get((FaceType*)this))._ffp[j];}
 
  char &FFi(const int j) {
   return (CAT< vector_occ<FaceType>,FFAdjTypeSup<T::FacePointer> >::Instance()->Get((FaceType*)this))._ffi[j];}  

  static bool HasFFAdjacency()   {   return true; }
  static bool HasFFAdjacencyOcc()   { return true; }

};

template <class T> class FFAdjOcc : public FFAdjOccBase<FFAdjTypeSup<typename T::FacePointer>,T>{};

template <class T> class VertexRefOcc: public T {
public:
 // typedef typename T::VertexType VertexType;
 // typedef typename T::VertexType::CoordType CoordType;

  inline typename T::VertexType *       & V( const int j ) 	     { assert(j>=0 && j<3); 
		return (CAT< vector_occ<FaceType>,VertexRef<T> >::Instance()->Get((FaceType*)this)).V(j); }

  inline typename T::VertexType * const & V( const int j ) const { assert(j>=0 && j<3); 
		return (CAT< vector_occ<FaceType>,VertexRef<T> >::Instance()->Get((FaceType*)this)).V(j); }

	inline typename T::VertexType * const  cV( const int j ) const { assert(j>=0 && j<3);	
		return (CAT< vector_occ<FaceType>,VertexRef<T> >::Instance()->Get((FaceType*)this)).V(j); }

	// Shortcut per accedere ai punti delle facce
	inline       typename T::CoordType & P( const int j ) 	    {	assert(j>=0 && j<3);	return	V(j)->P();	}
	inline const typename T::CoordType & P( const int j ) const	{	assert(j>=0 && j<3);	return  V(j)->cP(); }
	inline const typename T::CoordType &cP( const int j ) const	{	assert(j>=0 && j<3);	return  V(j)->cP(); }

	/** Return the pointer to the ((j+1)%3)-th vertex of the face.
		@param j Index of the face vertex.
	 */
	inline       typename T::VertexType *       &  V0( const int j )       { return V(j);}
	inline       typename T::VertexType *       &  V1( const int j )       { return V((j+1)%3);}
	inline       typename T::VertexType *       &  V2( const int j )       { return V((j+2)%3);}
	inline const typename T::VertexType * const &  V0( const int j ) const { return V(j);}
	inline const typename T::VertexType * const &  V1( const int j ) const { return V((j+1)%3);}
	inline const typename T::VertexType * const &  V2( const int j ) const { return V((j+2)%3);}
	inline const typename T::VertexType * const & cV0( const int j ) const { return cV(j);}
	inline const typename T::VertexType * const & cV1( const int j ) const { return cV((j+1)%3);}
	inline const typename T::VertexType * const & cV2( const int j ) const { return cV((j+2)%3);}

	/// Shortcut per accedere ai punti delle facce
	inline       typename T::CoordType &  P0( const int j )       { return V(j)->P();}
	inline       typename T::CoordType &  P1( const int j )       { return V((j+1)%3)->P();}
	inline       typename T::CoordType &  P2( const int j )       { return V((j+2)%3)->P();}
	inline const typename T::CoordType &  P0( const int j ) const { return V(j)->P();}
	inline const typename T::CoordType &  P1( const int j ) const { return V((j+1)%3)->P();}
	inline const typename T::CoordType &  P2( const int j ) const { return V((j+2)%3)->P();}
	inline const typename T::CoordType & cP0( const int j ) const { return cV(j)->P();}
	inline const typename T::CoordType & cP1( const int j ) const { return cV((j+1)%3)->P();}
	inline const typename T::CoordType & cP2( const int j ) const { return cV((j+2)%3)->P();}

	inline       typename T::VertexType *       & UberV( const int j )	      { assert(j>=0 && j<3); return v[j]; }
	inline const typename T::VertexType * const & UberV( const int j ) const	{ assert(j>=0 && j<3);	return v[j];	}
  static bool HasVertexRef()   { return true; }
};


  } // end namespace face
}// end namespace vcg
#endif
