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
Revision 1.6  2004/04/23 14:55:06  pietroni
conversion funtion

Revision 1.5  2004/03/10 00:59:06  cignoni
minor changes

Revision 1.4  2004/03/03 16:07:57  cignoni
Yet another cr lf mismatch

Revision 1.3  2004/02/24 21:36:39  cignoni
grouped documentation, changed typenames and reflection mechanism

Revision 1.2  2004/02/13 02:09:39  cignoni
First working release, with doxygen comment structure

Revision 1.1  2004/02/10 01:11:28  cignoni
Edited Comments and GPL license

****************************************************************************/


#ifndef VERTEX_TYPE 
#pragma message("\nYou should never directly include this file\_n")
#else

#include<vcg/space/point3.h>
#include<vcg/space/color4.h>
#include<vcg/space/tcoord2.h>

class DUMMYFACETYPE;

namespace vcg {

/**
    \ingroup vertex
    @name Vertex
    Class Vertex.
    This is the base class for definition of a vertex of the mesh.
	@param FLTYPE (Template Parameter) Specifies the scalar field of the vertex coordinate type.
	@param VFTYPE (Template Parameter) Specifies the type for the face, needed only for VF adjacency.
 */
template <class FLTYPE, class VFTYPE = DUMMYFACETYPE, class TCTYPE = TCoord2<float,1> > class VERTEX_TYPE
{
public:

	/// The scalar type used to represent coords (i.e. float, double, ...)
	typedef FLTYPE         ScalarType;
	/// The coordinate type used to represent the point (i.e. Point3f, Point3d, ...)
	typedef Point3<ScalarType> CoordType;
	typedef Point3<ScalarType> NormalType;
	/// The type base of the vertex, useful for recovering the original typename after user subclassing
	typedef VERTEX_TYPE    BaseVertexType;
	/// The type base of the vertex, useful for recovering the original typename after user subclassing
  typedef VFTYPE         FaceType;


/***********************************************/
/** @name Vertex Coords
    blah
    blah
**/
  //@{
protected:
	/// Spatial coordinates of the vertex
	CoordType _p;

public:
	/// Return the spatial coordinate of the vertex
	inline CoordType & P()
	{
	  assert( (_flags & DELETED) == 0 );
		assert( (_flags & NOTREAD) == 0 );
		assert( (_flags & NOTWRITE) == 0 );
		return _p;
	}

	/// Return the constant spatial coordinate of the vertex
	inline const CoordType & P() const
	{
		assert( (_flags & DELETED) == 0 );
		assert( (_flags & NOTREAD) == 0 );
		return _p;
	}

	/// Return the constant spatial coordinate of the vertex
	inline const CoordType & cP() const
	{
		assert( (_flags & DELETED) == 0 );
		assert( (_flags & NOTREAD) == 0 );
		return _p;
	}

	/// Return the spatial coordinate of the vertex, senza effettuare controlli sul flag
	inline CoordType & UberP()
	{
		return _p;
	}

	/// Return the constant spatial coordinate of the vertex, senza effettuare controlli sul flag
	inline const CoordType & UberP() const
	{
		return _p;
	}

  //@}

/***********************************************/
/** @name Vertex Flags
For each vertex we store a set of boolean values packed in a int. 
The default value for each flag is 0. Most commonly used flags are the \a deleted and the \a selected ones. 
Users can ask and dispose for a bit for their own purposes with the  vcg::VertexFull::NewUserBit() and vcg::VertexFull::DeleteUserBit() functions. 
The value returned by these functions has to be passed to the 
vcg::VertexFull::SetUserBit() vcg::VertexFull::ClearUserBit() and vcg::VertexFull::IsUserBit() functions to check and modify the obtained bit flag.

**/
//@{

protected:
	/// This are the flags of vertex, the default (reasonable) value is 0
	int _flags;		

public:
	/// Return the vector of _flags
	inline int & Flags ()
	{
			assert( (_flags & DELETED) == 0 );
			assert( (_flags & NOTREAD) == 0 );
			return _flags;
	}

	/// Return the vector of _flags, senza effettuare controlli sui bit
	inline int & UberFlags ()
	{
			return _flags;
	}
	inline const int UberFlags() const
	{
		return _flags;
	}

 	///  checks if the vertex is deleted
	bool IsD() const {return (_flags & DELETED) != 0;}
	///  checks if the vertex is readable
	bool IsR() const {return (_flags & NOTREAD) == 0;}
	///  checks if the vertex is modifiable
	bool IsW() const {return (_flags & NOTWRITE)== 0;}
	/// This funcion checks whether the vertex is both readable and modifiable
	bool IsRW() const {return (_flags & (NOTREAD | NOTWRITE)) == 0;}
	///  checks if the vertex is Modified
	bool IsS() const {return (_flags & SELECTED) != 0;}
	///  checks if the vertex is readable
	bool IsB() const {return (_flags & BORDER) != 0;}

	/** Set the flag value
		@param flagp Valore da inserire nel flag
	*/
	void SetFlags(int flagp) {_flags=flagp;}

	/** Set the flag value
		@param flagp Valore da inserire nel flag
	*/
	void ClearFlags() {_flags=0;}

	///  deletes the vertex from the mesh
	void SetD() {_flags |=DELETED;}
	///  un-delete a vertex
	void ClearD() {_flags &=(~DELETED);}
	///  marks the vertex as readable
	void SetR() {_flags &=(~NOTREAD);}
	///  marks the vertex as not readable
	void ClearR() {_flags |=NOTREAD;}
	///  marks the vertex as writable
	void ClearW() {_flags |=NOTWRITE;}
	///  marks the vertex as not writable
	void SetW() {_flags &=(~NOTWRITE);}
	///  select the vertex
	void SetS()		{_flags |=SELECTED;}
	/// Un-select a vertex
	void ClearS()	{_flags &= ~SELECTED;}
	void SetB()		{_flags |=BORDER;}
	void ClearB()	{_flags &=~BORDER;}
	
///  Return the first bit that is not still used
static int &LastBitFlag()
		{
			static int b =USER0;
			return b;
		}

/// allocate a bit among the flags that can be used by user.
static inline int NewUserBit()
		{
			LastBitFlag()=LastBitFlag()<<1;
			return LastBitFlag();
		}
// de-allocate a bit among the flags that can be used by user.
static inline bool DeleteUserBit(int bitval)
		{	
			if(LastBitFlag()==bitval) {
					LastBitFlag()= LastBitFlag()>>1;
					return true;
			}
			assert(0);
			return false;
		}
	/// This function checks if the given user bit is true
	bool IsUserBit(int userBit){return (_flags & userBit) != 0;}
	/// This function set  the given user bit 
	void SetUserBit(int userBit){_flags |=userBit;}
	/// This function clear the given user bit 
	void ClearUserBit(int userBit){_flags &= (~userBit);}
 
 
 //@}
  
  
/***********************************************/
/** @name Vertex Texture Coords
   blah
   blah
   **/
//@{

#ifdef __VCGLIB_VERTEX_T
protected:
	TCTYPE _t;
#endif

public:
	TCTYPE & T()
	{
#ifdef __VCGLIB_VERTEX_T
		return _t;
#else
		assert(0);
		return *(TCTYPE*)(&_flags);
#endif
	}

	const TCTYPE & T() const
	{
#ifdef __VCGLIB_VERTEX_T
		return _t;
#else
		assert(0);
		return *(TCTYPE*)(&_flags);
#endif
	}

//@}

/***********************************************/
/** @name Per vertex Color
   blah
   blah
   **/
//@{

#ifdef __VCGLIB_VERTEX_C
protected:
	Color4b _c;
#endif

public:
	Color4b & C()
	{
#ifdef __VCGLIB_VERTEX_C
		return _c;
#else
		assert(0);
		return *(Color4b*)(&_flags);
#endif
	}

	const Color4b & C() const
	{
#ifdef __VCGLIB_VERTEX_C
		return _c;
#else
		return Color4b(Color4b::White);
#endif
	}
 //@}

/***********************************************/
/** @name Vertex Quality
   blah
   blah
  **/
  //@{

#ifdef __VCGLIB_VERTEX_Q
protected:
	float _q;
#endif

public:
	float & Q()
	{
#ifdef __VCGLIB_VERTEX_Q
		return _q;
#else
		assert(0);
		return *(float*)(&_flags);
#endif
	}

	const float & Q() const
	{
#ifdef __VCGLIB_VERTEX_Q
		return _q;
#else
		return 1;
#endif
	}
 //@}

/***********************************************/
/** @name Vertex-Face Adjacency
   blah
   blah
 **/
 //@{

#if ((defined __VCGLIB_VERTEX_A) || (defined __VCGLIB_VERTEX_AS)) 
	// Puntatore ad una faccia appartenente alla stella del vertice, implementa l'adiacenza vertice-faccia
protected:
	VFTYPE *_fp;
	int _zp;
#endif

public:
inline VFTYPE * & Fp()
	{
#if ((defined __VCGLIB_VERTEX_A) || (defined __VCGLIB_VERTEX_AS))
		  return _fp;
#else
    assert(0);// you are probably trying to use VF topology in a vertex without it
		return *((VFTYPE **)(_flags));  
#endif
	}

inline const VFTYPE * & Fp() const
	{
#if ((defined __VCGLIB_VERTEX_A) || (defined __VCGLIB_VERTEX_AS))
		  return _fp;
#else
		assert(0);// you are probably trying to use VF topology in a vertex without it
		return (VFTYPE *)this;
#endif
	}

inline int & Zp()
	{
#if ((defined __VCGLIB_VERTEX_A) || (defined __VCGLIB_VERTEX_AS))

		  return _zp;
#else
    assert(0);// you are probably trying to use VF topology in a vertex without it
		return _flags;
#endif
	}

inline const int & Zp() const
	{
#if ((defined __VCGLIB_VERTEX_A) || (defined __VCGLIB_VERTEX_AS))
		  return _zp;
#else
		assert(0);// you are probably trying to use VF topology in a vertex without it
		return (void *)this;
#endif
	}



 //@}

/***********************************************/
/** @name Vertex Incremental Mark
   blah
   blah
 **/
 //@{

#ifdef __VCGLIB_VERTEX_M
protected:
  /// The incremental vertex mark
	int _imark;
#endif // Mark
public:
#ifdef __VCGLIB_VERTEX_M
	/// This function return the vertex incremental mark
	inline int & IMark()
	{
		assert( (_flags & DELETED) == 0 );
		assert( (_flags & NOTREAD) == 0 );
		assert( (_flags & NOTWRITE) == 0 );
		return _imark;
	}

	/// This function return the constant vertex incremental mark
	inline const int & IMark() const
	{
		assert( (_flags & DELETED) == 0 );
		assert( (_flags & NOTREAD) == 0 );
		return _imark;
	}
#endif

	/// Initialize the _imark system of the vertex
	inline void InitIMark()
	{
#ifdef __VCGLIB_VERTEX_M
		_imark = 0;
#endif
	}

 //@}
 
 /***********************************************/
 /** @name Vertex Normal 
   blah
 blah
 **/
 //@{

#ifdef __VCGLIB_VERTEX_N
protected:
    CoordType _n;
#endif 

public:
  /// Return the vertex normal
	inline CoordType & N()
	{
		assert( (_flags & DELETED) == 0 );
		assert( (_flags & NOTREAD) == 0 );
		assert( (_flags & NOTWRITE) == 0 );
#ifdef __VCGLIB_VERTEX_N
		return _n;
#else
		assert(0);
		return *(CoordType *)this;
#endif
	}

	/// Return the constant vertex normal
	inline const CoordType & N() const
	{
		assert( (_flags & DELETED) == 0 );
		assert( (_flags & NOTREAD) == 0 );
#ifdef __VCGLIB_VERTEX_N
		return _n;
#else
		assert(0);
		return *(CoordType *)this;
#endif
	}

	inline const CoordType  cN() const
	{
		assert( (_flags & DELETED) == 0 );
		assert( (_flags & NOTREAD) == 0 );
#ifdef __VCGLIB_VERTEX_N
		return _n;
#else
		return CoordType(0,0,0);
#endif
	}
   /// Return the Normal of the vertex
	inline CoordType & UberN()
	{
#ifdef __VCGLIB_VERTEX_N
		return _n;
#else
		assert(0);
		return *(CoordType *)this;
#endif
	}

	/// Return the constant normal of the vertex
	inline const CoordType & UberN() const
	{
#ifdef __VCGLIB_VERTEX_N
		return _n;
#else
		assert(0);
		return *(CoordType *)this;
#endif
	}
 //@}

 /***********************************************/
 /** @name Reflection Functions 
 Static functions  that give information about the current vertex type.
Reflection is a mechanism making it possible to investigate yourself. Reflection is used to investigate format of objects at runtime, invoke methods and access fields of these objects. Here we provide static const functions that are resolved at compile time and they give information about the data (normal, color etc.) supported by the current vertex type.
 **/
 //@{

static bool HasNormal()  { 
#ifdef __VCGLIB_VERTEX_N 
  return true;
#else
  return false;
#endif
}
static bool HasColor()  { 
#ifdef __VCGLIB_VERTEX_C 
  return true;
#else
  return false;
#endif
}
static bool HasMark()     { 
#ifdef __VCGLIB_VERTEX_M 
  return true;
#else
  return false;
#endif
}
static bool HasQuality()   { 
#ifdef __VCGLIB_VERTEX_Q 
  return true;
#else
  return false;
#endif
}
static bool HasTexture()   { 
#ifdef __VCGLIB_VERTEX_T 
  return true;
#else
  return false;
#endif
}
static bool HasVFAdjacency()   { 
#ifdef __VCGLIB_VERTEX_A 
  return true;
#else
  return false;
#endif
}
 //@}

/***********************************************/
 /** @Conversion to other vertex
 **/
 //@{

template <class VERT_TYPE>
inline Convert( const VERT_TYPE &v )
{
  this->P()=v->P();
  this._flags=v._flags;
	if (this->HasNormal())&&(v.HasNormal())
    this->N()=v->N();
  if (this->HasColor())&&(v.HasColor())
    this->C()=v->C();
  if (this->HasMark())&&(v.HasMark())
    this.IMark()=v.IMark();
  if (this->HasQuality())&&(v.HasQuality())
    this->Q()=v->Q();
  if (this->HasTexture())&&(v.HasTexture())
    this->T()=v->T();
}
	
 //@}

	enum { 
		// This bit indicate that the vertex is deleted from the mesh
		DELETED    = 0x0001,		// cancellato
		// This bit indicate that the vertex of the mesh is not readable
		NOTREAD    = 0x0002,		// non leggibile (ma forse modificabile) 
		// This bit indicate that the vertex is not modifiable
		NOTWRITE   = 0x0004,		// non modificabile (ma forse leggibile) 
		// This bit indicate that the vertex is modified
		MODIFIED   = 0x0008,		// modificato 
		// This bit can be used to mark the visited vertex
		VISITED    = 0x0010,		// Visited  
		// This bit can be used to select 
		SELECTED   = 0x0020,		// Selection flag
		// Border Flag
		BORDER     = 0x0100,
		// First user bit
		USER0      = 0x0200			// Fisrt user bit
			};


	/** Return the i-th spatial value of the vertex coordinate.
	    @param i Index of the spatial vertex coordinate (x=0 y=1 z=2).
	 */
	inline ScalarType & operator [] ( const int i ){
			assert(i>=0 && i<3);
			return P().V(i);
	}
	/** Return the i-th spatial value of the const vertex coordinate.
	    @param i Index of the spatial vertex coordinate (x=0 y=1 z=2).
	 */
	inline const FLTYPE & operator [] ( const int i ) const {
			assert(i>=0 && i<3);
			return P().V(i);
	}
	/// Operator to compare two vertices using lexicographic order
	inline bool operator < ( const VERTEX_TYPE & ve) const {
		return _p < ve._p;
		}
	inline VERTEX_TYPE() {
#ifdef _DEBUG 
		_flags=0;
#endif
	};

};


}	 // end namespace
#endif

