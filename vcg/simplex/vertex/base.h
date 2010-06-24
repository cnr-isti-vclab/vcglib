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
Revision 1.12  2008/03/17 11:39:14  ganovelli
added curvature and curvatruredir (compiled .net 2005 and gcc)

Revision 1.11  2008/02/04 21:26:49  ganovelli
added ImportData which imports all local attributes into vertexplus and faceplus.
A local attribute is everything (N(), C(), Q()....) except pointers to other simplices
(i.e. FFAdj, VFAdj, VertexRef) which are set to NULL.
Added some function for const attributes

Revision 1.10  2007/03/12 15:37:21  tarini
Texture coord name change!  "TCoord" and "Texture" are BAD. "TexCoord" is GOOD.

Revision 1.9  2007/02/12 19:00:56  ganovelli
added Name(std:vector<std::string>& n) that fills n with the names of the attribute of the vertex type

Revision 1.8  2006/09/28 17:34:11  cignoni
Added Missing GetBBox function

Revision 1.7  2006/02/27 17:42:43  ponchio
Added some documentation.

Revision 1.6  2005/12/05 15:58:10  cignoni
Removed spurious definition of flags in Aritymax that was overriding the correct definition in EmplyBitFlags and BitFlags classes

Revision 1.5  2005/12/02 00:44:41  cignoni
Reformatted and compacted flags code.

Revision 1.4  2005/11/16 22:59:35  cignoni
Standardized name of flags. It is plural becouse each simplex has many flag.

Revision 1.3  2005/11/12 18:36:51  cignoni
Added 'Visited' flag functions

Revision 1.2  lags2004/04/03 13:33:55  cignoni
Missing include

Revision 1.1  2004/03/29 08:36:26  cignoni
First working version!


****************************************************************************/
#ifndef __VCG_VERTEX_PLUS
#define __VCG_VERTEX_PLUS

//#include <vcg/space/point3.h>
#include <vcg/space/texcoord2.h>
#include <vcg/space/color4.h>
#include <vcg/complex/all_types.h>
#include <vcg/simplex/vertex/component.h>
//#include <vcg/complex/used_types.h>
#include <vcg/container/derivation_chain.h>

namespace vcg {

/*------------------------------------------------------------------*/ 
/* 
The base class of all the recusive definition chain. It is just a container of the typenames of the various simplexes.
These typenames must be known form all the derived classes.
*/
//
//template <class BVT, class BET, class BFT, class BTT>
//class VertexTypeHolder{
//  public:
//  typedef BVT VertType;
//  typedef BET EdgeType;
//  typedef BFT FaceType;
//  typedef BTT TetraType;
//  typedef BVT *VertPointer;
//  typedef BET *EdgePointer;
//  typedef BFT *FacePointer;
//  typedef BTT *TetraPointer;
//	template < class LeftV>
//	void ImportData(const LeftV  & /* left */ ) { }
//    static void Name(std::vector<std::string> & /* name */){}
//
//};

/* The base class form which we start to add our components.
it has the empty definition for all the standard members (coords, color flags)
Note:
in order to avoid both virtual classes and ambiguous definitions all 
the subsequent overrides must be done in a sequence of derivation.

In other words we cannot derive and add in a single derivation step 
(with multiple ancestor), both the real (non-empty) normal and color but 
we have to build the type a step a time (deriving from a single ancestor at a time). 


*/ 
//template <class UserUsedTypes>
//class VertexBase: public vertex::EmptyCore< UserUsedTypes >{
//};


/* The Real Big Vertex class;

The class __VertexArityMax__ is the one that is the Last to be derived,
and therefore is the only one to know the real members 
(after the many overrides) so all the functions with common behaviour 
using the members defined in the various Empty/nonEmpty component classes 
MUST be defined here. 

I.e. IsD() that uses the overridden Flags() member must be defined here.

*/

template <class UserTypes,
          template <typename> class A, template <typename> class B,
          template <typename> class C, template <typename> class D,
          template <typename> class E, template <typename> class F,
          template <typename> class G, template <typename> class H,
					template <typename> class I, template <typename> class J,
					template <typename> class K, template <typename> class L> 
class VertexArityMax: public Arity12<vertex::EmptyCore<UserTypes>, A, B, C, D, E, F, G, H, I, J, K, L> {

// ----- Flags stuff -----
public:



 	enum { 
		
		DELETED    = 0x0001,		// This bit indicate that the vertex is deleted from the mesh
		NOTREAD    = 0x0002,		// This bit indicate that the vertex of the mesh is not readable
		NOTWRITE   = 0x0004,		// This bit indicate that the vertex is not modifiable
		MODIFIED   = 0x0008,		// This bit indicate that the vertex is modified
		VISITED    = 0x0010,		// This bit can be used to mark the visited vertex
		SELECTED   = 0x0020,		// This bit can be used to select 
		BORDER     = 0x0100,    // Border Flag
		USER0      = 0x0200			// First user bit
  };

	inline int & UberFlags () { return this->Flags();	}
        inline int UberFlags() const 	{		return this->Flags();	}
 	
	bool IsD() const {return (this->Flags() & DELETED) != 0;} ///  checks if the vertex is deleted
	bool IsR() const {return (this->Flags() & NOTREAD) == 0;} ///  checks if the vertex is readable
	bool IsW() const {return (this->Flags() & NOTWRITE)== 0;}///  checks if the vertex is modifiable
	bool IsRW() const {return (this->Flags() & (NOTREAD | NOTWRITE)) == 0;}/// This funcion checks whether the vertex is both readable and modifiable
	bool IsS() const {return (this->Flags() & SELECTED) != 0;}///  checks if the vertex is Selected
	bool IsB() const {return (this->Flags() & BORDER) != 0;}///  checks if the vertex is a border one
	bool IsV() const {return (this->Flags() & VISITED) != 0;}///  checks if the vertex Has been visited
	

	/** Set the flag value
		@param flagp Valore da inserire nel flag
	*/
	void SetFlags(int flagp) {this->Flags()=flagp;}

	/** Set the flag value
		@param flagp Valore da inserire nel flag
	*/
	void ClearFlags() {this->Flags()=0;}
	void SetD() {this->Flags() |=DELETED;}///  deletes the vertex from the mesh
	void ClearD() {this->Flags() &=(~DELETED);}///  un-delete a vertex
	void SetR() {this->Flags() &=(~NOTREAD);}///  marks the vertex as readable
	void ClearR() {this->Flags() |=NOTREAD;}///  marks the vertex as not readable
	void ClearW() {this->Flags() |=NOTWRITE;}///  marks the vertex as writable
	void SetW() {this->Flags() &=(~NOTWRITE);}///  marks the vertex as not writable
	void SetS()		{this->Flags() |=SELECTED;}///  select the vertex
	void ClearS()	{this->Flags() &= ~SELECTED;}/// Un-select a vertex
	void SetB()		{this->Flags() |=BORDER;}
	void ClearB()	{this->Flags() &=~BORDER;}
	void SetV()		{this->Flags() |=VISITED;}
	void ClearV()	{this->Flags() &=~VISITED;}
	
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
	bool IsUserBit(int userBit){return (this->Flags() & userBit) != 0;}
	/// This function set  the given user bit 
	void SetUserBit(int userBit){this->Flags() |=userBit;}
	/// This function clear the given user bit 
	void ClearUserBit(int userBit){this->Flags() &= (~userBit);}

 template<class BoxType>
  void GetBBox( BoxType & bb ) const
  {	  bb.Set(this->P());  }

          };

          
/*

These are the three main classes that are used by the library user to define its own vertexes.
The user MUST specify the names of all the type involved in a generic complex.
so for example when defining a vertex of a trimesh you must know the name of the type of the edge and of the face.
Typical usage example:

A vertex with coords, flags and normal for use in a standard trimesh:

class VertexNf   : public VertexSimp2< VertexNf, EdgeProto, FaceProto, vert::Coord3d, vert::Flag, vert::Normal3f  > {};


A vertex with coords, and normal for use in a tetrahedral mesh AND in a standard trimesh:

class TetraVertex   : public VertexSimp3< TetraVertex, EdgeProto, FaceProto, TetraProto, vert::Coord3d, vert::Normal3f  > {};


A summary of the available vertex attributes (see component.h for more details):
          
Coord3f,  Coord3d, 
Normal3s,  Normal3f,  Normal3d
Mark                              //a int component (incremental mark)
BitFlags
TexCoord2s,  TexCoord2f,  TexCoord2d
Color4b
Qualitys, Qualityf, Qualityd
VFAdj                             //topology (vertex->face adjacency)
*/

template <class UserTypes,
          template <typename> class A = DefaultDeriver, template <typename> class B = DefaultDeriver,
          template <typename> class C = DefaultDeriver, template <typename> class D = DefaultDeriver,
          template <typename> class E = DefaultDeriver, template <typename> class F = DefaultDeriver,
          template <typename> class G = DefaultDeriver, template <typename> class H = DefaultDeriver,
					template <typename> class I = DefaultDeriver, template <typename> class J = DefaultDeriver,
					template <typename> class K = DefaultDeriver, template <typename> class L = DefaultDeriver>
							class Vertex: public VertexArityMax<UserTypes, A, B, C, D, E, F, G, H, I, J, K, L>  {
			 public: typedef AllTypes::AVertexType IAm; typedef UserTypes TypesPool;};

}// end namespace
#endif
