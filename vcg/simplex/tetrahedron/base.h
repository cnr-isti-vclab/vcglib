/****************************************************************************
* VCGLib                                                            o o     *
* Visual and Computer Graphics Library                            o     o   *
*                                                                _   O  _   *
* Copyright(C) 2004-2016                                           \/)\/    *
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
Revision 1.1  2007/05/09 10:31:53  ganovelli
added



****************************************************************************/
// #ifndef __VCG_TETRA_MESH
// #error "This file should not be included alone. It is automatically included by complex.h"
// #endif
#ifndef __VCG_TETRA_PLUS
#define __VCG_TETRA_PLUS

//#include <vcg/space/point3.h>
//#include <vcg/space/texcoord2.h>
//#include <vcg/space/color4.h>
//#include <vcg/simplex/tetrahedron/component.h>

namespace vcg {

/*------------------------------------------------------------------*/ 

// /* The base class form which we start to add our components.
// it has the empty definition for all the standard members (coords, color flags)
// Note:
// in order to avoid both virtual classes and ambiguous definitions all 
// the subsequent overrides must be done in a sequence of derivation.

// In other words we cannot derive and add in a single derivation step 
// (with multiple ancestor), both the real (non-empty) normal and color but 
// we have to build the type a step a time (deriving from a single ancestor at a time). 




template <class UserTypes>
                class TetraTypeHolder: public UserTypes {
  public:

    template <class RightT>
    void ImportData(const RightT & ){}
    static void Name(std::vector<std::string> & /* name */){}


 // prot
    inline int VN()  const { return 4;}
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
template <class UserTypes>
class TetraSimpBase: public
            tetrahedron::EmptyCore< TetraTypeHolder <UserTypes> > {
};

template <class UserTypes,
          template <typename> class A, template <typename> class B,
          template <typename> class C, template <typename> class D,
          template <typename> class E, template <typename> class F,
          template <typename> class G, template <typename> class H,
          template <typename> class I, template <typename> class J,
          template <typename> class K, template <typename> class L>
class TetraArityMax: public Arity12<TetraSimpBase<UserTypes>, A, B, C, D, E, F, G, H, I, J, K, L> {

// ----- Flags stuff -----
public:
  
   	enum { 
		
		DELETED     = 0x00000001,		// Tet is deleted from the mesh
		NOTREAD     = 0x00000002,		// Tet of the mesh is not readable
		NOTWRITE    = 0x00000004,		// Tet of the mesh is not writable
        VISITED     = 0x00000010,		// Tet has been visited. Usualy this is a per-algorithm used bit.
		SELECTED    = 0x00000020,		// Tet is selected. Algorithms should try to work only on selected face (if explicitly requested)
		// Border _flags, it is assumed that BORDERi = BORDER0<<i 
		BORDER0     = 0x00000040,
		BORDER1     = 0x00000080,
		BORDER2     = 0x00000100,
		BORDER3     = 0x00000200,
		BORDER0123  = BORDER0 | BORDER1 | BORDER2 | BORDER3,
		// Crease _flags,  it is assumed that FEATUREi = FEATURE0<<i 
		// First user bit
		USER0       = 0x00004000
			};

 
 	///  checks if the Tet is deleted
        bool IsD() const {return (this->cFlags() & DELETED) != 0;}
	///  checks if the Tet is readable
        bool IsR() const {return (this->cFlags() & NOTREAD) == 0;}
	///  checks if the Tet is modifiable
        bool IsW() const {return (this->cFlags() & NOTWRITE)== 0;}
	/// This funcion checks whether the Tet is both readable and modifiable
        bool IsRW() const {return (this->cFlags() & (NOTREAD | NOTWRITE)) == 0;}
	///  checks if the Tet is Modified
        bool IsS() const {return (this->cFlags() & SELECTED) != 0;}
	///  checks if the Tet is Modified
        bool IsV() const {return (this->cFlags() & VISITED) != 0;}
	
	/** Set the flag value
		@param flagp Valore da inserire nel flag
	*/
	void SetFlags(int flagp) {this->Flags()=flagp;}

	/** Set the flag value
		@param flagp Valore da inserire nel flag
	*/
	void ClearFlags() {this->Flags()=0;}

	///  deletes the Tet from the mesh
	void SetD() {this->Flags() |=DELETED;}
	///  un-delete a Tet
	void ClearD() {this->Flags() &=(~DELETED);}
	///  marks the Tet as readable
	void SetR() {this->Flags() &=(~NOTREAD);}
	///  marks the Tet as not readable
	void ClearR() {this->Flags() |=NOTREAD;}
	///  marks the Tet as writable
	void SetW() {this->Flags() &=(~NOTWRITE);}
	///  marks the Tet as notwritable
	void ClearW() {this->Flags() |=NOTWRITE;}
	///  select the Tet
	void SetS()		{this->Flags() |=SELECTED;}
	/// Un-select a Tet
  void ClearS()	{this->Flags() &= ~SELECTED;}
	///  select the Tet
	void SetV()		{this->Flags() |=VISITED;}
	/// Un-select a Tet
  void ClearV()	{this->Flags() &= ~VISITED;}
	
	/// This function checks if the face is border
	bool IsB(int i) const {return (this->cFlags() & (BORDER0<<i)) != 0;}
	bool IsAnyB() const {return (this->cFlags() & (BORDER0123)) != 0;}
	/// This function select the face
  void SetB(int i)		{this->Flags() |=(BORDER0<<i);}
	/// This funcion execute the inverse operation of SetS()
	void ClearB(int i)	{this->Flags() &= (~(BORDER0<<i));}
	
	///  Return the first bit that is not still used
	static int &FirstUnusedBitFlag()
	{
	  static int b =USER0;
	  return b;
	}

	/// Allocate a bit among the flags that can be used by user. It updates the FirstUnusedBitFlag.
	static inline int NewBitFlag()
	{
	  int bitForTheUser = FirstUnusedBitFlag();
	  FirstUnusedBitFlag()=FirstUnusedBitFlag()<<1;
	  return bitForTheUser;
	}

	/// De-allocate a pre allocated bit. It updates the FirstUnusedBitFlag.
	// Note you must deallocate bit in the inverse order of the allocation (as in a stack)
	static inline bool DeleteBitFlag(int bitval)
	{
	  if(FirstUnusedBitFlag()>>1==bitval) {
		FirstUnusedBitFlag() = FirstUnusedBitFlag()>>1;
		return true;
	  }
	  assert(0);
	  return false;
	}

	/// This function checks if the given user bit is true
	bool IsUserBit(int userBit){return (this->Flags() & userBit) != 0;}

	/// This function set the given user bit
	void SetUserBit(int userBit){this->Flags() |=userBit;}

	/// This function clear the given user bit
	void ClearUserBit(int userBit){this->Flags() &= (~userBit);}

 template<class BoxType>
  void GetBBox( BoxType & bb ) const
  {
	  bb.Set(this->cP(0));
	  bb.Add(this->cP(1));
	  bb.Add(this->cP(2));
	  bb.Add(this->cP(3));
  }


};

// template < typename T=int>
// class TetraDefaultDeriver : public T {};
          
/*

These are the three main classes that are used by the library user to define its own Facees.
The user MUST specify the names of all the type involved in a generic complex.
so for example when defining a Face of a trimesh you must know the name of the type of the edge and of the face.
Typical usage example:

A Face with coords, flags and normal for use in a standard trimesh:

class MyFaceNf   : public FaceSimp2< VertProto, EdgeProto, MyFaceNf, face::Flag, face::Normal3f  > {};


A Face with coords, and normal for use in a tetrahedral mesh AND in a standard trimesh:

class TetraFace   : public FaceSimp3< VertProto, EdgeProto, TetraFace, TetraProto, face::Coord3d, face::Normal3f  > {};


A summary of the components that can be added to a face (see components.h for details):
          
VertexRef
Mark                                            //Incremental mark (int)
VTAdj                                           //Topology vertex face adjacency
                                                 (pointers to next face in the ring of the vertex
TTAdj                                           //topology: face face adj
                                                  pointers to adjacent faces

*/

template <class UserTypes,
          template <typename> class A = DefaultDeriver, template <typename> class B = DefaultDeriver,
          template <typename> class C = DefaultDeriver, template <typename> class D = DefaultDeriver,
          template <typename> class E = DefaultDeriver, template <typename> class F = DefaultDeriver,
          template <typename> class G = DefaultDeriver, template <typename> class H = DefaultDeriver,
          template <typename> class I = DefaultDeriver, template <typename> class J = DefaultDeriver,
          template <typename> class K = DefaultDeriver, template <typename> class L = DefaultDeriver>
class TetraSimp : public TetraArityMax<UserTypes, A, B, C, D, E, F, G, H, I, J, K, L>  {
			 public: typedef AllTypes::ATetraType IAm; typedef UserTypes TypesPool;};

}// end namespace
#endif

