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
Revision 1.4  2004/04/28 11:37:14  pietroni
*** empty log message ***

Revision 1.3  2004/04/26 09:38:54  pietroni
*** empty log message ***

Revision 1.2  2004/04/20 12:42:37  pietroni
*** empty log message ***

Revision 1.1  2004/04/15 08:54:20  pietroni
*** empty log message ***


****************************************************************************/

#ifndef TETRA_TYPE 
#pragma message("\nYou should never directly include this file\_n")
#else
#define NULL 0
#include<vcg/space/point3.h>
#include<vcg/space/tetra3.h>

namespace vcg {
/**
    \ingroup tetrahedron
    @name Tetrahedron
    Class Tetrahedron.
    This is the base class for definition of a Tetrahedron of the mesh.
	@param VTYPE (Template Parameter) Specifies the type for the vertex.
 */

template < class VTYPE >
class TETRA_TYPE{

public:

  ///	The base type of the face
	typedef TETRA_TYPE BaseTetraType;
	/// The vertex type 
	typedef typename VTYPE VertexType;
  /// The coordinate type used to represent the point (i.e. Point3f, Point3d, ...)
  typedef typename VertexType::CoordType CoordType;
  /// The scalar type used to represent coords (i.e. float, double, ...)
	typedef typename VertexType::ScalarType ScalarType;

 
/***********************************************/
/** @name Tetrahedron Flags
For each Tetrahedron we store a set of boolean values packed in a int. 
The default value for each flag is 0. Most commonly used flags are the \a deleted and the \a selected ones. 
Users can ask and dispose for a bit for their own purposes with the  vcg::TetrahedronFull::NewUserBit() and vcg::TetrahedronFull::DeleteUserBit() functions. 
The value returned by these functions has to be passed to the 
vcg::TetrahedronFull::SetUserBit() vcg::TetrahedronFull::ClearUserBit() and vcg::TetrahedronFull::IsUserBit() functions to check and modify the obtained bit flag.

**/
//@{

/// This are the flags of tetrahedron, the default value is 0
	int  _flags;
 
  enum {
	DELETED     = 0x00000001,	 // deleted tetrahedron flag
	SELECTED		= 0x00000002,	 // Selection flag
	BORDERF0    = 0x00000004,  // Border flag, Face 0
	BORDERF1    = 0x00000008,  // Border flag, Face 1
	BORDERF2    = 0x00000010,  // Border flag, Face 2
  BORDERF3    = 0x00000020,  // Border flag, Face 3
  BORDERE0    = 0x00000040,  // Border flag, Edge 0
	BORDERE1    = 0x00000080,  // Border flag, Edge 1
	BORDERE2    = 0x00000100,  // Border flag, Edge 2
  BORDERE3    = 0x00000200,  // Border flag, Edge 3
  BORDERE4    = 0x00000400,  // Border flag, Edge 4
  BORDERE5    = 0x00000800,  // Border flag, Edge 5
	USER0				= 0x00001000,  // new flag for user
	};

public:
  	/// Return the vector of _flags
	inline int & Flags ()
	{
			assert( (_flags & DELETED) == 0 );
			return _flags;
	}

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

/// This function checks if the given user bit is true.
bool IsUserBit(int userBit){return (_flags & userBit) != 0;}
/// This function set  the given user bit.
void SetUserBit(int userBit){_flags |=userBit;}
/// This function clear the given user bit.
void ClearUserBit(int userBit){_flags &= (~userBit);}
/// This function checks if the tetrahedron is deleted.
bool IsD() const {return (_flags & DELETED) != 0;}
/// This function mark the tetrahedron as deleted.
void SetD()		{_flags |=DELETED;}
/// This function mark the tetrahedron as not deleted.
void ClearD() {_flags &=~DELETED;}
/// This function mark the tetrahedron as selected.
void SetS() {_flags |=SELECTED;}
/// This function mark the tetrahedron as not selected.
void ClearS() {_flags &=~SELECTED;}
/// This function return true if one face is extern.
bool HaveBorderF() {return ((_flags & (BORDERF0 | BORDERF1 | BORDERF2 | BORDERF3)) != 0);}
/// This function return true if the face is extern.
bool IsBorderF(int face) {
  assert ((face<4)&&(face>-1));
  return (this->TTp(face) == this);
}
 //@}

/***********************************************/
/** @name Vertex Pointers
For each Tetrahedron we store 4 pointers to vertex
**/
//@{
 /// The 4 vertices of the tetrahedron
protected:
	VertexType *_v[4];
public:
 
/** Return the pointer to the j-th vertex of the terahedron.
		@param j Index of the tetrahedron's vertex.
	 */
	inline VertexType * & V( const int j )
	{	
		assert( (_flags & DELETED) == 0 );
		assert(j >= 0);
		assert(j <  4);
		return _v[j];
	}

	inline const VertexType * const & V( const int j ) const
	{
		assert( (_flags & DELETED) == 0 );
		assert(j>=0);
		assert(j<4);
		return _v[j];
	}
  
	inline const VertexType * const & cV( const int j ) const
	{
		assert( (_flags & DELETED) == 0 );
		assert(j>=0);
		assert(j<4);
		return _v[j];
	}

/***********************************************/
/** @name Topology Structures
For each Tetrahedron we store 2 array for Tatrahedron - Tetrahedron topology ( sharing Face)
and 2 array to implement the list of Vertex - Tetrahedron Topology (List of Tetrahedron sharing a vertex).
**/
//@{


#ifdef __VCGLIB_TETRA_AT
protected:
  ///pointers to tetrahedron for tetrahedron-tetrahedron topology (sharing same face)
	TETRA_TYPE *_ttp[4];
  ///index of face for tetrahedron-tetrahedron topology (sharing same face)
	int _tti[4]; 
public:
  ///Function to access the Tetrahedron that share the index-face (extern face returns a pointer to himself)
  	TETRA_TYPE *&TTp(const int &index)
	{
		return _ttp[index];
	}
  ///Function to see the index of the face as seen from the other tetrahedron (extern face returns -1)
	int &TTi(const int &index)
	{
		return _tti[index];
	}
#endif 

#ifdef __VCGLIB_TETRA_AV
protected:
	///pointers to tetrahedron for vertex-tetrahedron topology (sharing same vertex)
	TETRA_TYPE *_tvp[4];
  ///index of vertex for vertex-tetrahedron topology (sharing same vertex)
	short int _tvi[4];
public:
  ///Function to access the Next Tetrahedron of the list that share the index-face (end of list is Null)
  	TETRA_TYPE *&TVp(const int &index)
	{
		return _tvp[index];
	}
  ///Function to see the index of the Vertex as seen from the next tetrahedron of the list ( end of list is -1)
	short int &TVi(const int &index)
	{
		return _tvi[index];
	}
#endif
//@}

/***********************************************/
/** @Default Tatrahedron Functions**/
//@{
public:

  ///Constructor
	TETRA_TYPE()
	{	
		_flags=0;
	}
  ///initialize default parameters of tetrahedron
	virtual 	void Init(VertexType * p0,VertexType * p1,VertexType * p2,VertexType * p3)
				{
					_flags = 0;
					_v[0]=p0;
					_v[1]=p1;
					_v[2]=p2;
					_v[3]=p3;

					if(ComputeVolume()<0 )
						swap(_v[1],_v[2]);

#ifdef		__VCGLIB_TETRA_TA
					_z[0]=_z[1]=_z[2]=_z[3]=-1;
					_t[0]=_t[1]=_t[2]=_t[3]=NULL;
#endif
#ifdef		__VCGLIB_TETRA_TV
					_zv[0]=_zv[1]=_zv[2]=_zv[3]=-1;
					_tv[0]=_tv[1]=_tv[2]=_tv[3]=NULL;
#endif					
			}
 ///set border vertices using TT-topology
#ifdef __VCGLIB_TETRA_AT
	void setBorderV()
	{	
		int i;
		for (i=0;i<4;i++)
			if (T(i)==this)
			{
				FV(i,0)->SetB();
				FV(i,1)->SetB();
				FV(i,2)->SetB();
			}
	}
#endif
//@}

/***********************************************/
/** @Generic geometric and quality funtions of a tetrahedron**/
//@{
#ifdef __VCGLIB_TETRA_TN
private:
  CoordType _n[4];
public:
#endif
///return the normal of a face of the tetrahedron
	const CoordType & N(const int &i){
    assert((i>=0)&&(i<4));
#ifdef __VCGLIB_TETRA_TN		
    return _n[i];
#else	  
    Tetra3<ScalarType> T=Tetra3<ScalarType>();
    T.P0(0)=V(0)->P();
    T.P1(0)=V(1)->P();
    T.P2(0)=V(2)->P();
    T.P3(0)=V(3)->P();
    return (Normal<Tetra3<ScalarType> >(T,i));
#endif
		}	

 /// Calculate the normal to all the faces of a tetrahedron, the value is store in a position of vecton _n for each face
void ComputeNormal() 
{
#ifdef __VCGLIB_TETRA_TN
   Tetra3<ScalarType> T=Tetra3<ScalarType>();
   T.P0(0)=V(0)->P();
   T.P1(0)=V(1)->P();
   T.P2(0)=V(2)->P();
   T.P3(0)=V(3)->P();
   
	for (int i=0;i<4;i++)
      _n[i]=(Normal<Tetra3<ScalarType> >(T,i));
#else
	assert(0);
#endif
}
//@}

/***********************************************/
/** @Generic geometric and quality funtions of a tetrahedron**/
//@{

#ifdef __VCGLIB_TETRA_TQ
		ScalarType _volume;	
		ScalarType _aspect_ratio; 
#endif


 	ScalarType ComputeVolume(){
      Tetra3<ScalarType> T=Tetra3<ScalarType>();
       T.P0(0)=V(0)->cP();
       T.P1(0)=V(1)->cP();
       T.P2(0)=V(2)->cP();
       T.P3(0)=V(3)->cP();
      #ifdef __VCGLIB_TETRA_TQ
			_volume = T.ComputeVolume();
      return _volume;
      #else
       return (T.ComputeVolume());
      #endif
		}

	///return the volume of the tetrahedron
	const double & Volume(){
#ifdef __VCGLIB_TETRA_TQ
		return _volume;
#else
		return (( V(2)->cP()-V(0)->cP())^(V(1)->cP()-V(0)->cP() ))*(V(3)->cP()-V(0)->cP())/6.0;
#endif
		}	
  ///return aspect ratio of the tetrahedron
	double AspectRatio(){
#ifdef __VCGLIB_TETRA_TQ
		return _aspect_ratio;
#else
		return ComputeAspectRatio();
#endif
		}
//@}


 /***********************************************/
 /** @name Reflection Functions 
 Static functions  that give information about the current tetra type.
Reflection is a mechanism making it possible to investigate yourself. Reflection is used to investigate format of objects at runtime, invoke methods and access fields of these objects. Here we provide static const functions that are resolved at compile time and they give information about the data supported by the current tetra type.
 **/
 //@{

static bool HasTetraNormal()  { 
#ifdef __VCGLIB_TETRA_TN 
  return true;
#else
  return false;
#endif
}
static bool HasTetraMark()  { 
#ifdef __VCGLIB_TETRA_TM
  return true;
#else
  return false;
#endif
}
static bool HasTetraQuality()  { 
#ifdef __VCGLIB_TETRA_TQ 
  return true;
#else
  return false;
#endif
}

static bool HasTTAdjacency()  { 
#if (defined(__VCGLIB_TETRA_AT) || defined(__VCGLIB_TETRA_SAT))
  return true;
#else
  return false;
#endif
}
static bool HasVTAdjacency()  { 
#if (defined(__VCGLIB_TETRA_AV) || defined(__VCGLIB_TETRA_SAT))
  return true;
#else
  return false;
#endif
}

//@}
};//end class

}//end namespace
#endif