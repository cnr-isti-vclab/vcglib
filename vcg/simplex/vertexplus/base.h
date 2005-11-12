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
Revision 1.2  2004/04/03 13:33:55  cignoni
Missing include

Revision 1.1  2004/03/29 08:36:26  cignoni
First working version!


****************************************************************************/
#ifndef __VCG_VERTEX_PLUS
#define __VCG_VERTEX_PLUS

#include <vcg/space/point3.h>
#include <vcg/space/tcoord2.h>
#include <vcg/space/color4.h>
#include <vcg/simplex/vertexplus/component.h>

namespace vcg {

  class DumET {};
  class DumFT {};
  class DumTT {};

/*------------------------------------------------------------------*/ 
/* 
The base class of all the recusive definition chain. It is just a container of the typenames of the various simplexes.
These typenames must be known form all the derived classes.
*/

template <class BVT, class BET, class BFT, class BTT>
class VertexTypeHolder{
  public:
  typedef BVT VertType;
  typedef BET EdgeType;
  typedef BFT FaceType;
  typedef BTT TetraType;
  typedef BVT *VertPointer;
  typedef BET *EdgePointer;
  typedef BFT *FacePointer;
  typedef BTT *TetraPointer;
};

/* The base class form which we start to add our components.
it has the empty definition for all the standard members (coords, color flags)
Note:
in order to avoid both virtual classes and ambiguous definitions all 
the subsequent overrides must be done in a sequence of derivation.

In other words we cannot derive and add in a single derivation step 
(with multiple ancestor), both the real (non-empty) normal and color but 
we have to build the type a step a time (deriving from a single ancestor at a time). 


*/ 
template <class BVT, class BET=DumET, class BFT=DumFT, class BTT=DumTT>
class VertexBase: public vert::EmptyTexture<
                         vert::EmptyVFAdj<
                         vert::EmptyColor<
                         vert::EmptyQuality<
                         vert::EmptyNormal<
                         vert::EmptyFlag<
                         vert::EmptyCoord<
                            VertexTypeHolder <BVT, BET, BFT, BTT> > > > > > > >{
};

// Metaprogramming Core

template <class BVT, class BET, class BFT,class BTT,
          template <typename> class A> 
          class VertexArity1: public A<VertexBase<BVT,BET,BFT,BTT> > {
          };

template <class BVT, class BET, typename BFT, class BTT,
          template <typename> class A, template <typename> class B> 
          class VertexArity2: public B<VertexArity1<BVT,BET,BFT,BTT, A> > {};

template <class BVT, class BET, typename BFT,class BTT,
          template <typename> class A, template <typename> class B, 
          template <typename> class C > 
          class VertexArity3: public C<VertexArity2<BVT,BET,BFT,BTT, A, B> > {};

template <class BVT, class BET, typename BFT,class BTT,
          template <typename> class A, template <typename> class B, 
          template <typename> class C, template <typename> class D> 
          class VertexArity4: public D<VertexArity3<BVT,BET,BFT,BTT, A, B, C> > {};

template <class BVT, class BET, typename BFT,class BTT,
          template <typename> class A, template <typename> class B, 
          template <typename> class C, template <typename> class D,
          template <typename> class E > 
          class VertexArity5: public E<VertexArity4<BVT,BET,BFT,BTT, A, B, C, D> > {};
template <class BVT, class BET, typename BFT,class BTT,
          template <typename> class A, template <typename> class B, 
          template <typename> class C, template <typename> class D,
          template <typename> class E, template <typename> class F > 
          class VertexArity6: public F<VertexArity5<BVT,BET,BFT,BTT, A, B, C, D, E> > {};

/* The Real Big Vertex class;

The class __VertexArityMax__ is the one that is the Last to be derived,
and therefore is the only one to know the real members 
(after the many overrides) so all the functions with common behaviour 
using the members defined in the various Empty/nonEmpty component classes 
MUST be defined here. 

I.e. IsD() that uses the overridden Flags() member must be defined here.

*/

template <class BVT, class BET, typename BFT,class BTT,
          template <typename> class A, template <typename> class B, 
          template <typename> class C, template <typename> class D, 
          template <typename> class E, template <typename> class F,
          template <typename> class G> 
class VertexArityMax: public G<VertexArity6<BVT,BET,BFT,BTT, A, B, C, D, E, F> > {

// ----- Flags stuff -----
public:
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


	inline int & UberFlags ()
	{
			return Flags();
	}
	inline const int UberFlags() const
	{
		return Flags();
	}

 	///  checks if the vertex is deleted
	bool IsD() const {return (Flags() & DELETED) != 0;}
	///  checks if the vertex is readable
	bool IsR() const {return (Flags() & NOTREAD) == 0;}
	///  checks if the vertex is modifiable
	bool IsW() const {return (Flags() & NOTWRITE)== 0;}
	/// This funcion checks whether the vertex is both readable and modifiable
	bool IsRW() const {return (Flags() & (NOTREAD | NOTWRITE)) == 0;}
	///  checks if the vertex is Selected
	bool IsS() const {return (Flags() & SELECTED) != 0;}
	///  checks if the vertex is a border one
	bool IsB() const {return (Flags() & BORDER) != 0;}
	///  checks if the vertex Has been visited
	bool IsV() const {return (Flags() & VISITED) != 0;}

	/** Set the flag value
		@param flagp Valore da inserire nel flag
	*/
	void SetFlags(int flagp) {Flags()=flagp;}

	/** Set the flag value
		@param flagp Valore da inserire nel flag
	*/
	void ClearFlags() {Flags()=0;}

	///  deletes the vertex from the mesh
	void SetD() {Flags() |=DELETED;}
	///  un-delete a vertex
	void ClearD() {Flags() &=(~DELETED);}
	///  marks the vertex as readable
	void SetR() {Flags() &=(~NOTREAD);}
	///  marks the vertex as not readable
	void ClearR() {Flags() |=NOTREAD;}
	///  marks the vertex as writable
	void ClearW() {Flags() |=NOTWRITE;}
	///  marks the vertex as not writable
	void SetW() {Flags() &=(~NOTWRITE);}
	///  select the vertex
	void SetS()		{Flags() |=SELECTED;}
	/// Un-select a vertex
	void ClearS()	{Flags() &= ~SELECTED;}
	void SetB()		{Flags() |=BORDER;}
	void ClearB()	{Flags() &=~BORDER;}
	void SetV()		{Flags() |=VISITED;}
	void ClearV()	{Flags() &=~VISITED;}
	
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
	bool IsUserBit(int userBit){return (Flags() & userBit) != 0;}
	/// This function set  the given user bit 
	void SetUserBit(int userBit){Flags() |=userBit;}
	/// This function clear the given user bit 
	void ClearUserBit(int userBit){Flags() &= (~userBit);}

          };


template < typename T=int>
class DefaultDeriver : public T {};
          
/*

These are the three main classes that are used by the library user to define its own vertexes.
The user MUST specify the names of all the type involved in a generic complex.
so for example when defining a vertex of a trimesh you must know the name of the type of the edge and of the face.
Typical usage example:

A vertex with coords, flags and normal for use in a standard trimesh:

class VertexNf   : public VertexSimp2< VertexNf, EdgeProto, FaceProto, vert::Coord3d, vert::Flag, vert::Normal3f  > {};


A vertex with coords, and normal for use in a tetrahedral mesh AND in a standard trimesh:

class TetraVertex   : public VertexSimp3< TetraVertex, EdgeProto, FaceProto, TetraProto, vert::Coord3d, vert::Normal3f  > {};

*/

template <class BVT, class BET, class BFT, class BTT,
          template <typename> class A = DefaultDeriver, template <typename> class B = DefaultDeriver,
          template <typename> class C = DefaultDeriver, template <typename> class D = DefaultDeriver,
          template <typename> class E = DefaultDeriver, template <typename> class F = DefaultDeriver,
          template <typename> class G = DefaultDeriver > 
              class VertexSimp3: public VertexArityMax<BVT,BET,BFT,BTT, A, B, C, D, E, F, G>  {};

template <class BVT, class BET, class BFT, 
          template <typename> class A = DefaultDeriver, template <typename> class B = DefaultDeriver,
          template <typename> class C = DefaultDeriver, template <typename> class D = DefaultDeriver,
          template <typename> class E = DefaultDeriver, template <typename> class F = DefaultDeriver,
          template <typename> class G = DefaultDeriver > 
              class VertexSimp2: public VertexArityMax<BVT,BET,BFT,DumTT, A, B, C, D, E, F, G>  {};

template <class BVT, class BET, 
          template <typename> class A = DefaultDeriver, template <typename> class B = DefaultDeriver,
          template <typename> class C = DefaultDeriver, template <typename> class D = DefaultDeriver,
          template <typename> class E = DefaultDeriver, template <typename> class F = DefaultDeriver,
          template <typename> class G = DefaultDeriver > 
                class VertexSimp1: public VertexArityMax<BVT,BET,DumFT,DumTT, A, B, C, D, E, F, G>  {};

}// end namespace
#endif
