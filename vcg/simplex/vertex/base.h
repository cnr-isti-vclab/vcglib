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
Revision 1.30  2007/02/26 14:21:44  turini
VTb moved to VTp

Revision 1.29  2007/02/20 14:08:34  ganovelli
added QualityType to comply vertexplus type

Revision 1.28  2006/08/23 15:34:20  marfr960
added minimal comments

Revision 1.26  2005/11/12 18:41:14  cignoni
Added HasFlags and initialization of flags at construction.

Revision 1.25  2005/10/14 13:25:50  cignoni
Added cVFp member

Revision 1.24  2005/10/06 14:26:39  pietroni
added getBBox method

Revision 1.23  2005/03/18 16:38:36  fiorin
Minor changes

Revision 1.22  2005/03/18 00:13:45  cignoni
Removed NormalizedNormalV  (out of standard and wrong) and
added the member functions Normal and NormalizedNormal() (just like for faces)

Revision 1.21  2004/10/28 00:50:49  cignoni
Better Doxygen documentation

Revision 1.20  2004/10/11 17:45:05  ganovelli
added template on corrdinate type (default Point3)

Revision 1.19  2004/09/28 15:24:56  fiorin
DUMMY classes definition moved into vcg namespace

Revision 1.18  2004/08/25 15:15:27  ganovelli
minor changes to comply gcc compiler (typename's and stuff)

Revision 1.17  2004/07/20 15:24:53  pietroni
corrected NormalizedNormalV function...

Revision 1.16  2004/07/15 11:25:01  ganovelli
VFb moved to VFp, userbit to bitflag,setV, inclusion of pos.h

Revision 1.15  2004/07/15 10:13:48  pietroni
adde NormalizedNormalV funtion to compute the normal on a vertex

Revision 1.14  2004/05/13 22:44:40  ganovelli
syntax error (typo)

Revision 1.13  2004/05/13 22:40:02  ganovelli
default template parameters

Revision 1.12  2004/05/13 12:49:22  pietroni
no default template parameters... each one must be specified

Revision 1.12  2004/05/10 13:31:13  ganovelli
function for edge adjacency added

$Log: not supported by cvs2svn $
Revision 1.30  2007/02/26 14:21:44  turini
VTb moved to VTp

Revision 1.29  2007/02/20 14:08:34  ganovelli
added QualityType to comply vertexplus type

Revision 1.28  2006/08/23 15:34:20  marfr960
added minimal comments

Revision 1.26  2005/11/12 18:41:14  cignoni
Added HasFlags and initialization of flags at construction.

Revision 1.25  2005/10/14 13:25:50  cignoni
Added cVFp member

Revision 1.24  2005/10/06 14:26:39  pietroni
added getBBox method

Revision 1.23  2005/03/18 16:38:36  fiorin
Minor changes

Revision 1.22  2005/03/18 00:13:45  cignoni
Removed NormalizedNormalV  (out of standard and wrong) and
added the member functions Normal and NormalizedNormal() (just like for faces)

Revision 1.21  2004/10/28 00:50:49  cignoni
Better Doxygen documentation

Revision 1.20  2004/10/11 17:45:05  ganovelli
added template on corrdinate type (default Point3)

Revision 1.19  2004/09/28 15:24:56  fiorin
DUMMY classes definition moved into vcg namespace

Revision 1.18  2004/08/25 15:15:27  ganovelli
minor changes to comply gcc compiler (typename's and stuff)

Revision 1.17  2004/07/20 15:24:53  pietroni
corrected NormalizedNormalV function...

Revision 1.16  2004/07/15 11:25:01  ganovelli
VFb moved to VFp, userbit to bitflag,setV, inclusion of pos.h

Revision 1.15  2004/07/15 10:13:48  pietroni
adde NormalizedNormalV funtion to compute the normal on a vertex

Revision 1.14  2004/05/13 22:44:40  ganovelli
syntax error (typo)

Revision 1.13  2004/05/13 22:40:02  ganovelli
default template parameters

Revision 1.12  2004/05/13 12:49:22  pietroni
no default template parameters... each one must be specified

Revision 1.11  2004/05/10 13:31:13  ganovelli
function for edge adjacency added

Revision 1.10  2004/05/10 13:13:17  cignoni
added void to Convert, corrected return object in VFp

Revision 1.9  2004/05/06 15:28:10  pietroni
changed names to VF topology function (was missed)

Revision 1.8  2004/05/05 17:03:25  pietroni
changed name to topology functions

Revision 1.7  2004/04/28 11:37:14  pietroni
*** empty log message ***

Revision 1.6  2004/04/26 09:40:15  pietroni
*** empty log message ***

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

#pragma message("VCGLIB Warning: this way to define the simplex vertex is DEPRECATED  and no more SUPPORTED")
#pragma message("VCGLIB Warning: use vcg/simplex/vertexplus instead ")

#include<vcg/space/point3.h>
#include<vcg/space/color4.h>
#include<vcg/space/texcoord2.h>
#include<vcg/simplex/face/pos.h>
#include<vcg/space/box3.h>


namespace vcg {

	class DUMMYFACETYPE;
	class DUMMYEDGETYPE;
	class DUMMYTETRATYPE;
/** \addtogroup vertex */
//@{
/*!
 * This class represent the generic configurable Vertex; 
 * Usually you never direclty use this class with this name but you build 
 * your own type by directly including one of the .h files under the face/with 
 * directory. Each file specify a class type with the desired fields. So for example 
 * including 'vcg/simplex/vertex/with/VCVN.h' allow you to use the class VertVCVN that has per-vertex color and normal stored inside.
 */
template <class FLTYPE, class VETYPE = DUMMYEDGETYPE, class VFTYPE = DUMMYFACETYPE, class VTTYPE = DUMMYTETRATYPE,class TCTYPE = TexCoord2<float,1>, class CoordTYPE= Point3<FLTYPE> > 
class VERTEX_TYPE
{
public:

	/// The scalar type used to represent coords (i.e. float, double, ...)
	typedef FLTYPE         ScalarType;
	/// The coordinate type used to represent the point (i.e. Point3f, Point3d, ...)
	typedef CoordTYPE CoordType;
	typedef Point3<ScalarType> NormalType;
	/// The type base of the vertex, useful for recovering the original typename after user subclassing
	typedef VERTEX_TYPE    BaseVertexType;
	/// The type of the face pointed by the vertex if vertex edge topology is present
	typedef VETYPE         EdgeType;
	/// The type base of the vertex, useful for recovering the original typename after user subclassing
	typedef VFTYPE         FaceType;
	/// The type of the quality (same as scalar)
	typedef  ScalarType  QualityType;


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
	///  checks if the vertex is visited
	bool IsV() const {return (_flags & VISITED) != 0;}


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
	/// Set vertex as ob border
	void SetB()		{_flags |=BORDER;}
	void ClearB()	{_flags &=~BORDER;}
	///  checks if the vertex is visited
	void ClearV()	{_flags &= ~VISITED;}
	///  checks if the vertex is visited
	void SetV()		{_flags |=VISITED;}
	///  Return the first bit that is not still used
static int &LastBitFlag()
		{
			static int b =USER0;
			return b;
		}

/// allocate a bit among the flags that can be used by user.
static inline int NewBitFlag()
		{
			LastBitFlag()=LastBitFlag()<<1;
			return LastBitFlag();
		}
// de-allocate a bit among the flags that can be used by user.
static inline bool DeleteBitFlag(int bitval)
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
  
  
/*#*******************	
*  Bounding box *
**********************/

void GetBBox( Box3<ScalarType> & bb ) const
{
	bb.Set( cP() );
}

/***********************************************/
/** @name Vertex Texture Coords
   blah
   blah
   **/
//@{

#ifdef __VCGLIB_VERTEX_VT
protected:
	TCTYPE _t;
#endif

public:
	TCTYPE & T()
	{
#ifdef __VCGLIB_VERTEX_VT
		return _t;
#else
		assert(0);
		return *(TCTYPE*)(&_flags);
#endif
	}

	const TCTYPE & T() const
	{
#ifdef __VCGLIB_VERTEX_VT
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

#ifdef __VCGLIB_VERTEX_VC
protected:
	Color4b _c;
#endif

public:
	Color4b & C()
	{
#ifdef __VCGLIB_VERTEX_VC
		return _c;
#else
		assert(0);
		return *(Color4b*)(&_flags);
#endif
	}

	const Color4b & C() const
	{
#ifdef __VCGLIB_VERTEX_VC
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

#ifdef __VCGLIB_VERTEX_VQ
protected:
	float _q;
#endif

public:
	float & Q()
	{
#ifdef __VCGLIB_VERTEX_VQ
		return _q;
#else
		assert(0);
		return *(float*)(&_flags);
#endif
	}

	const float & Q() const
	{
#ifdef __VCGLIB_VERTEX_VQ
		return _q;
#else
		return 1;
#endif
	}
 //@}
/** @name Vertex-Edge Adjacency
   blah
   blah
 **/
 //@{

#if ((defined __VCGLIB_VERTEX_EA) || (defined __VCGLIB_VERTEX_EAS)) 
	// Puntatore ad un edge appartenente alla stella del vertice, implementa l'adiacenza vertice-edge
protected:
	EdgeType *_ep;
	int _ei;
#endif

public:
inline EdgeType * & VEp()
	{
#if ((defined __VCGLIB_VERTEX_EA) || (defined __VCGLIB_VERTEX_EAS))
		  return _ep;
#else
    assert(0);// you are probably trying to use VF topology in a vertex without it
		return *((EdgeType **)(_flags));  
#endif
	}

inline const EdgeType * & VEp() const
	{
#if ((defined __VCGLIB_VERTEX_EA) || (defined __VCGLIB_VERTEX_EAS))
		  return _ep;
#else
		assert(0);// you are probably trying to use VF topology in a vertex without it
		return (EdgeType *)this;
#endif
	}

inline int & VEi()
	{
#if ((defined __VCGLIB_VERTEX_EA) || (defined __VCGLIB_VERTEX_EAS))

		  return _ei;
#else
    assert(0);// you are probably trying to use VF topology in a vertex without it
		return _flags;
#endif
	}

inline const int & VEi() const
	{
#if ((defined __VCGLIB_VERTEX_EA) || (defined __VCGLIB_VERTEX_EAS))
		  return _ei;
#else
		assert(0);// you are probably trying to use VF topology in a vertex without it
		return (void *)this;
#endif
	}
 //@}
/***********************************************/
/** @name Vertex-Face Adjacency
   blah
   blah
 **/
 //@{

#if ((defined __VCGLIB_VERTEX_AF) || (defined __VCGLIB_VERTEX_AFS)) 
	// Puntatore ad una faccia appartenente alla stella del vertice, implementa l'adiacenza vertice-faccia
protected:
	VFTYPE *_vfb;
	int _vfi;
#endif

public:
inline VFTYPE * & VFp()
	{
#if ((defined __VCGLIB_VERTEX_AF) || (defined __VCGLIB_VERTEX_AFS))
		  return _vfb;
#else
    assert(0);// you are probably trying to use VF topology in a vertex without it
    static VFTYPE *dum;
    return dum;
#endif
	}

inline const VFTYPE * & VFp() const
	{
#if ((defined __VCGLIB_VERTEX_AF) || (defined __VCGLIB_VERTEX_AFS))
		  return _vfb;
#else
		assert(0);// you are probably trying to use VF topology in a vertex without it
		return (VFTYPE *)0;
#endif
	}

inline const VFTYPE *  cVFp() const
	{
#if ((defined __VCGLIB_VERTEX_AF) || (defined __VCGLIB_VERTEX_AFS))
		  return _vfb;
#else
		assert(0);// you are probably trying to use VF topology in a vertex without it
		return (VFTYPE *)0;
#endif
	}

inline int & VFi()
	{
#if ((defined __VCGLIB_VERTEX_AF) || (defined __VCGLIB_VERTEX_AFS))

		  return _vfi;
#else
    assert(0);// you are probably trying to use VF topology in a vertex without it
		return _flags;
#endif
	}

inline const int & VFi() const
	{
#if ((defined __VCGLIB_VERTEX_AF) || (defined __VCGLIB_VERTEX_AFS))
		  return _vfi;
#else
		assert(0);// you are probably trying to use VF topology in a vertex without it
		return (void *)this;
#endif
	}



 //@}

  /***********************************************/
/** @name Vertex-Tetrahedron Adjacency
   blah
   blah
 **/
 //@{

#if ((defined __VCGLIB_VERTEX_AT) || (defined __VCGLIB_VERTEX_ATS)) 
	// Pointer to first tetrahedron of the start implements the Vertex-Tetrahedron Topology
protected:
	VTTYPE *_vtp;
	int _vti;
#endif

public:
inline VTTYPE * & VTp()
	{
#if ((defined __VCGLIB_VERTEX_AT) || (defined __VCGLIB_VERTEX_ATS))
		  return _vtp;
#else
    assert(0);// you are probably trying to use VF topology in a vertex without it
		return *((VTTYPE **)(_flags));  
#endif
	}

inline const VTTYPE * &  VTp() const
	{
#if ((defined __VCGLIB_VERTEX_AT) || (defined __VCGLIB_VERTEX_ATS))
		  return _vtp;
#else
		assert(0);// you are probably trying to use VF topology in a vertex without it
		return (VTTYPE *)this;
#endif
	}

inline int & VTi()
	{
#if ((defined __VCGLIB_VERTEX_AT) || (defined __VCGLIB_VERTEX_ATS))

		  return _vti;
#else
    assert(0);// you are probably trying to use VF topology in a vertex without it
		return _flags;
#endif
	}

inline const int & VTi() const
	{
#if ((defined __VCGLIB_VERTEX_AT) || (defined __VCGLIB_VERTEX_ATS))
		  return _vti;
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

#ifdef __VCGLIB_VERTEX_VM
protected:
  /// The incremental vertex mark
	int _imark;
#endif // Mark
public:
#ifdef __VCGLIB_VERTEX_VM
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
#ifdef __VCGLIB_VERTEX_VM
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

#ifdef __VCGLIB_VERTEX_VN
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
#ifdef __VCGLIB_VERTEX_VN
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
#ifdef __VCGLIB_VERTEX_VN
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
#ifdef __VCGLIB_VERTEX_VN
		return _n;
#else
		return CoordType(0,0,0);
#endif
	}
   /// Return the Normal of the vertex
	inline CoordType & UberN()
	{
#ifdef __VCGLIB_VERTEX_VN
		return _n;
#else
		assert(0);
		return *(CoordType *)this;
#endif
	}

	/// Return the constant normal of the vertex
	inline const CoordType & UberN() const
	{
#ifdef __VCGLIB_VERTEX_VN
		return _n;
#else
		assert(0);
		return *(CoordType *)this;
#endif
	}
  
template <bool NormalizeFlag> 
const CoordType GenericNormal()
{
  if (!HasVFAdjacency())
  {
    assert(0);
		return (VERTEX_TYPE::CoordType (0,0,0));
  }
  else
  {
    vcg::face::VFIterator<typename VERTEX_TYPE::FaceType> VFi=vcg::face::VFIterator<typename VERTEX_TYPE::FaceType>();
    VFi.f=VFp();
    VFi.z=VFi();
    typename VERTEX_TYPE::CoordType N= typename VERTEX_TYPE::CoordType(0,0,0);
    while (!VFi.End())
    {
      N+=VFi.f->Normal();
      VFi++;
    }
    if(NormalizeFlag) N.Normalize();
    return N;
  }
}

/// Return the un-normalized value of the vertex normal as it correspond to the current geometry.
/// It is always computed and never use any stored value. 
/// REQUIRES vertex-face topology
const CoordType Normal()           { return GenericNormal<false>(); }

/// Return the normalized value of the vertex normal as it correspond to the current geometry.
/// It is always computed and never use any stored value. 
/// REQUIRES vertex-face topology
const CoordType NormalizedNormal() { return GenericNormal<true>(); }


 //@}

 /***********************************************/
 /** @name Reflection Functions 
 Static functions  that give information about the current vertex type.
Reflection is a mechanism making it possible to investigate yourself. Reflection is used to investigate format of objects at runtime, invoke methods and access fields of these objects. Here we provide static const functions that are resolved at compile time and they give information about the data (normal, color etc.) supported by the current vertex type.
 **/
 //@{

static bool HasFlags()  { // Note the plural because ONE vertex has many Flags (but just one color, normal, mark, quality ecc.)
  return true;
}

static bool HasNormal()  { 
#ifdef __VCGLIB_VERTEX_VN 
  return true;
#else
  return false;
#endif
}
static bool HasColor()  { 
#ifdef __VCGLIB_VERTEX_VC 
  return true;
#else
  return false;
#endif
}
static bool HasMark()     { 
#ifdef __VCGLIB_VERTEX_VM 
  return true;
#else
  return false;
#endif
}
static bool HasQuality()   { 
#ifdef __VCGLIB_VERTEX_VQ 
  return true;
#else
  return false;
#endif
}
static bool HasTexCoord()   { 
#ifdef __VCGLIB_VERTEX_VT 
  return true;
#else
  return false;
#endif
}
static bool HasVFAdjacency()   { 
#ifdef  __VCGLIB_VERTEX_AF
  return true;
#else
  return false;
#endif
}
static bool HasVTAdjacency()   { 
#ifdef  __VCGLIB_VERTEX_AT
  return true;
#else
  return false;
#endif
}

static bool HasVEAdjacency()   { 
#ifdef __VCGLIB_VERTEX_EA 
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
inline void Convert( VERT_TYPE &v )
{
  P()=v.P();
  Flags()=v.Flags();
	if ((HasNormal())&&(v.HasNormal()))
    N()=v.N();
  if ((HasColor())&&(v.HasColor()))
    C()=v.C();
#ifdef __VCGLIB_VERTEX_VM
  if ((HasMark())&&(v.HasMark()))
    IMark()=v.IMark();
#endif
  if ((HasQuality())&&(v.HasQuality()))
    Q()=v.Q();
  if ((HasTexCoord())&&(v.HasTexCoord()))
    T()=v.T();
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
//#ifdef _DEBUG 
		_flags=0;
//#endif
  };

};

//@}
}	 // end namespace
#endif

