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

#ifndef TETRA_TYPE 
#pragma message("\nYou should never directly include this file\_n")
#else
#define TETRA_TYPE
#include<vcg/space/point3.h>
#include<vcg/space/tetra4.h>

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
bool HaveBorderF() 
{
  {return ((_flags & (BORDERF0 | BORDERF1 | BORDERF2 | BORDERF3)) != 0);}
{
/// This function return true if the face is extern.
bool IsBorderF(int face) {
  assert ((face<4)&&(face>-1));
  switch (face) {
  case 0:
   {return ((_flags & BORDERF0) != 0);}
    break;
  case 1:
    {return ((_flags & BORDERF1) != 0);}
    break;
  case 2:
    {return ((_flags & BORDERF2) != 0);}
    break;
   case 3:
    {return ((_flags & BORDERF3) != 0);}
    break;
  }
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
  /// The Functions to access a vertex
  	inline  MVTYPE * &V(int index) 
	{
		return _v[index];
	}
  /// The Functions to access a vertex
	inline  const MVTYPE * &V(int index) const
	{
		return _v[index];
	}
 //@}


/***********************************************/
/** @name Topology Structures
For each Tetrahedron we store 2 array for Tatrahedron - Tetrahedron topology ( sharing Face)
and 2 array to implement the list of Vertex - Tetrahedron Topology (List of Tetrahedron sharing a vertex).
**/
//@{


#ifdef __VCGLIB_TETRA_A
protected:
  ///pointers to tetrahedron for tetrahedron-tetrahedron topology (sharing same face)
	TETRA_TYPE *_t[4];
  ///index of face for tetrahedron-tetrahedron topology (sharing same face)
	int _z[4]; 
public:
  ///Function to access the Tetrahedron that share the index-face (extern face returns a pointer to himself)
  	TETRA_TYPE *&T(const int &index)
	{
		return t[index];
	}
  ///Function to see the index of the face as seen from the other tetrahedron (extern face returns -1)
	int &Z(const int &index)
	{
		return z[index];
	}
#endif 

#ifdef __VCGLIB_TETRA_V
protected:
	///pointers to tetrahedron for vertex-tetrahedron topology (sharing same vertex)
	TETRA_TYPE *_tv[4];
  ///index of vertex for vertex-tetrahedron topology (sharing same vertex)
	short int _zv[4];
public:
  ///Function to access the Next Tetrahedron of the list that share the index-face (end of list is Null)
  	TETRA_TYPE *&TV(const int &index)
	{
		return t[index];
	}
  ///Function to see the index of the Vertex as seen from the next tetrahedron of the list ( end of list is -1)
	int &ZV(const int &index)
	{
		return z[index];
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
	virtual 	void Init(MVTYPE * p0,MVTYPE * p1,MVTYPE * p2,MVTYPE * p3)
				{
					_flags = 0;
					_v[0]=p0;
					_v[1]=p1;
					_v[2]=p2;
					_v[3]=p3;

					if(ComputeVolume()<0 )
						swap(_v[1],_v[2]);

#ifdef		__VCGLIB_TETRA_A
					z[0]=z[1]=z[2]=z[3]=-1;
					t[0]=t[1]=t[2]=t[3]=NULL;
#endif
#ifdef		__VCGLIB_TETRA_V
					zv[0]=zv[1]=zv[2]=zv[3]=-1;
					tv[0]=tv[1]=tv[2]=tv[3]=NULL;
#endif					
			}
 ///set border vertices using TT-topology
#ifdef __VCGLIB_TETRA_A
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

#ifdef __VCGLIB_TETRA_Q
		scalar_type _volume;	
		scalar_type _aspect_ratio; 
#endif

#ifdef __VCGLIB_TETRA_Q
 	scalar_type ComputeVolume(){
      Tetra4<scalar_type> T(V(0)->cP(),V(1)->cP(),V(2)->cP(),V(3)->cP());
			_volume = T.ComputeVolume();
		}
#endif

	///return the volume of the tetrahedron
	const double & Volume(){
#ifdef __VCGLIB_TETRA_Q
		return _volume;
#else
		return (( V(2)->cP()-V(0)->cP())^(V(1)->cP()-V(0)->cP() ))*(V(3)->cP()-V(0)->cP())/6.0;
#endif
		}	
  ///return aspect ratio of the tetrahedron
	double AspectRatio(){
#ifdef __VCGLIB_TETRA_Q
		return _aspect_ratio;
#else
		return ComputeAspectRatio();
#endif
		}
//@}

/***********************************************/
/** @Tatrahedron Functions to retrieve information about relation between faces of tetrahedron
(faces,adges,vertices).**/
//@{

	MVTYPE *FV(const int &indexF,const int &indexV)
	{	int facevert[4][3]={{0,1,2},
					{0,3,1},
					{0,2,3},
					{1,3,2}};
		assert ((indexF<4)&&(indexV<3));
		return _v[facevert[indexF][indexV]];
	}

	int iFV(const int &indexF,const int &indexV)
	{	int facevert[4][3]={{0,1,2},
					{0,3,1},
					{0,2,3},
					{1,3,2}};

		assert ((indexF<4)&&(indexV<3));
		return facevert[indexF][indexV];
	}

	int VF(const int &indexV,const int &indexF)
	{	int vertface[4][3]={{0,1,2},
					{0,3,1},
					{0,2,3},
					{1,3,2}};

		assert ((indexV<4)&&(indexF<3));
		return vertface[indexV][indexF];
	}

	int FVE(const int &indexF,const int &indexV)
	{	
	
		int facevertedge[4][4]={{1,0,3,-1},
						{2,0,-1,4},
						{-1,3,5,4},
						{1,-1,5,2}
						};
		assert ((indexF<4)&&(indexV<4));
		return facevertedge[indexF][indexV];
	}

	inline  MVTYPE * VE(const int &indexE,const int &indexV) 
	{	int edgevert[6][2]={{0,1},
					{0,2},
					{0,3},
					{1,2},
					{1,3},
					{2,3}};
		assert ((indexE<6)&&(indexV<2));
		return _v[edgevert[indexE][indexV]];
	}
// Tetrahedron Vertex (Couple) To Edge Conversion Function
	inline int VEC(const int &indexV1,const int &indexV2) 
		{	int ve[4][4]={{-1,  0,  1,  2},
						  { 0, -1,  3,  4},
						  { 1,  3, -1,  5},
						  { 2,  4,  5, -1}};

			assert ((indexV1<4)&&(indexV2<4));
			return ve[indexV1][indexV2];
		}

	inline int EV(const int &indexV,const int &indexE) 
	{	
		int vertedge[4][3]={{0,1,2},
											{0,3,4},
											{5,1,3},
											{4,5,2}};
		assert ((indexE<3)&&(indexV<4));
		return vertedge[indexV][indexE];
	}
	
	inline int iVE(const int &indexE,const int &indexV) 
	{	int edgevert[6][2]={{0,1},
					{0,2},
					{0,3},
					{1,2},
					{1,3},
					{2,3}};
		assert ((indexE<6)&&(indexV<2));
		return edgevert[indexE][indexV];
	}

	inline  int FE(const int &indexE,const int &indexSide) 
	{	int edgeface[6][2]={{0,1},
					{0,2},
					{1,2},
					{0,3},
					{1,3},
					{2,3}};

		assert ((indexE<6)&&(indexSide<2));
		return edgeface [indexE][indexSide];
	}

	inline  int EF(const int &indexF,const int &faceindexEdge) 
	{	int faceedge[4][3]={{0,3,1},
					{2,4,0},
					{1,5,2},
					{4,5,3}
					};

		assert ((indexF<4)&&(faceindexEdge<3));
		return faceedge [indexF][faceindexEdge];
	}

//@}




	

};//end class

}//end namespace
#endif