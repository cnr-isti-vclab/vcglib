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

/** @name Vertex
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
	/// The type base of the vertex, useful for recovering the original typename after user subclassing
	typedef VERTEX_TYPE    BaseVertexType;
	/// The type base of the vertex, useful for recovering the original typename after user subclassing
  typedef VFTYPE         face_type;


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
  blah
  blah
**/
//@{

protected:
	/// This are the _flags of vertex, the default value is 0
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



	enum {
		OBJ_TYPE_N =  0x0001,
		OBJ_TYPE_M =  0x0002,
		OBJ_TYPE_A =  0x0004,
		OBJ_TYPE_AS = 0x0008,
		OBJ_TYPE_C  = 0x0010,
		OBJ_TYPE_T  = 0x0020,
		OBJ_TYPE_Q  = 0x0040,
	};

	enum {
		OBJ_TYPE = 
#ifdef __VCGLIB_VERTEX_N
		OBJ_TYPE_N |
#endif
#ifdef __VCGLIB_VERTEX_M
		OBJ_TYPE_M |
#endif
#ifdef __VCGLIB_VERTEX_A
		OBJ_TYPE_A |
#endif
#ifdef __VCGLIB_VERTEX_AS
		OBJ_TYPE_AS |
#endif
#ifdef __VCGLIB_VERTEX_C
		OBJ_TYPE_C |
#endif
#ifdef __VCGLIB_VERTEX_T
		OBJ_TYPE_T |
#endif
#ifdef __VCGLIB_VERTEX_Q
		OBJ_TYPE_Q |
#endif
		0
	};
	

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

/*
Queste funzioni servono per ottenere a runtime un bit per i flag

*/	
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

	/// This function checks if the vertex is deleted
	bool IsD() const {return (_flags & DELETED) != 0;}
	/// This function checks if the vertex is readable
	bool IsR() const {return (_flags & NOTREAD) == 0;}
	/// This function checks if the vertex is modifiable
	bool IsW() const {return (_flags & NOTWRITE)== 0;}
	/// This funcion checks whether the vertex is both readable and modifiable
	bool IsRW() const {return (_flags & (NOTREAD | NOTWRITE)) == 0;}
	/// This function checks if the vertex is Modified
	bool IsM() const {return (_flags & MODIFIED)!= 0;}
	/// This function checks if the vertex is marked as visited
	bool IsV() const {return (_flags & VISITED) != 0;}
	/// This function checks if the vertex is selected
	bool IsS() const {return (_flags & SELECTED) != 0;}
	/// This function checks if the vertex is readable
	bool IsB() const {return (_flags & BORDER) != 0;}
//	bool IsMF() const {return (_flags & NOTMANIFOLD) == 0;}

	/// This function checks if the vertex is deleted from the mesh
	bool IsDeleted() const {return IsD();}
	/// This function checks if the vertex is readable
	bool IsReadable() const {return IsR();}
	/** Set the flag value
		@param flagp Valore da inserire nel flag
	*/
	void SetFlags(int flagp) {_flags=flagp;}

	/// This function deletes the vertex from the mesh
	void SetD() {_flags |=DELETED;}
	/// This funcion execute the inverse operation of SetD()
	void ClearD() {_flags &=(~DELETED);}
	/// This function marks the vertex as modified. It's necessary to mark all modified vertex to have a consistent mesh
	void SetM() {_flags |=MODIFIED;}
	/// This function marks the vertex as not modified
	void ClearM() {_flags &=(~MODIFIED);}
	/// This function marks the vertex as readable
	void SetR() {_flags &=(~NOTREAD);}
	/// This function marks the vertex as not readable
	void ClearR() {_flags |=NOTREAD;}
	/// This function marks the vertex as writable
	void ClearW() {_flags |=NOTWRITE;}
	/// This function marks the vertex as not writable
	void SetW() {_flags &=(~NOTWRITE);}
	/// This funcion marks the vertex as visited
	void SetV() {_flags |=VISITED;}
	/// This function marks the vertex as not visited. This flag, initially, is setted to random value, therefore, to the beginnig of every function it is necessary to clean up the flag
	void ClearV() {_flags &=(~VISITED);}
	/// This function select the vertex
	void SetS()		{_flags |=SELECTED;}
	/// This funcion execute the inverse operation of SetS()
	void ClearS()	{_flags &= ~SELECTED;}
	void SetB()		{_flags |=BORDER;}
	void ClearB()	{_flags &=~BORDER;}
	
	/// This function checks if the given user bit is true
	bool IsUserBit(int userBit){return (_flags & userBit) != 0;}
	/// This function set  the given user bit 
	void SetUserBit(int userBit){_flags |=userBit;}
	/// This function clear the given user bit 
	void ClearUserBit(int userBit){_flags &= (~userBit);}
};


}	 // end namespace
#endif

