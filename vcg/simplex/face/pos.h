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
Revision 1.32  2007/10/17 19:46:50  cignoni
Added I() access function for the z member to the pos

Revision 1.31  2007/05/28 14:09:41  fiorin
Added Set method which takes a face pointer and a vertex pointer.

Revision 1.30  2007/05/16 15:11:32  fiorin
Replaced ambigous StarSize method with NumberOfIncidentVertices and NumberOfIncidentFaces

Revision 1.29  2007/04/20 12:40:31  cignoni
Corrected  V() operator. It was plainly wrong. Luckly enough it was not very used

Revision 1.28  2007/01/11 10:37:08  cignoni
Added include assert.h

Revision 1.27  2007/01/02 10:06:53  giec
Added access functions F()

Revision 1.26  2006/12/29 13:13:00  giec
Corrected wrong assert in V(i) access function

Revision 1.25  2006/12/04 16:06:12  cignoni
Added FFlip() and const VFlip() operators

Revision 1.24  2006/11/13 01:57:23  cignoni
Added a missing prototype to ismanifold

Revision 1.23  2006/11/09 17:22:56  cignoni
Added ismanifold

Revision 1.22  2006/10/07 14:24:26  cignoni
Explained the use of V() operator of a pos

Revision 1.21  2006/09/25 09:57:49  cignoni
Better comment on usage of VF iterators

Revision 1.20  2005/12/15 11:57:48  corsini
Replace Pos<FaceType> with PosType

Revision 1.19  2005/12/15 11:19:00  corsini
Fix operators

Revision 1.18  2005/12/15 10:53:16  corsini
Add constructor which takes as input a face and a vertex

Revision 1.17  2005/10/16 23:30:39  ponchio
IsBorder(...) declaration needed.

Revision 1.16  2005/10/13 09:29:10  cignoni
Removed the reference to Deprecated f->IsBorder(i) now everyone should use IsBorder(*f,i);

Revision 1.15  2005/01/03 11:22:31  cignoni
Added better documentation (with an example and the V0 V1 V2 access members

Revision 1.14  2004/10/28 00:50:48  cignoni
Better Doxygen documentation

Revision 1.13  2004/10/18 17:14:42  ganovelli
error FFP -> FFp

Revision 1.12  2004/09/14 19:46:10  ganovelli
constructor added

Revision 1.11  2004/08/25 15:15:27  ganovelli
minor changes to comply gcc compiler (typename's and stuff)

Revision 1.10  2004/07/27 09:47:49  cignoni
Added V() access function instead of V(0)

Revision 1.9  2004/07/18 07:45:30  cignoni
Removed two const modifiers from the VFIterator

Revision 1.8  2004/07/15 12:03:07  ganovelli
minor changes

Revision 1.7  2004/07/15 11:28:44  ganovelli
basefacetype to facetype

Revision 1.6  2004/07/06 06:25:44  cignoni
changed the VFIterator ++ to return a facepointer instead of a bool

Revision 1.5  2004/06/02 16:25:45  ganovelli
changed F(.. to FFp
changed Z(   to FFi(

Revision 1.4  2004/05/10 15:21:47  cignoni
Added a constructor without vertex pointer

Revision 1.3  2004/05/10 13:41:57  cignoni
Added VFIterator

Revision 1.2  2004/03/12 15:22:28  cignoni
Written some documentation and added to the trimes doxygen module

Revision 1.1  2004/03/10 08:32:30  cignoni
Initial commit


****************************************************************************/

/** \file face/pos.h
 * Definition of vcg:face::Pos class.
 * This file contain the definition of vcg::face::Pos class and the derived vcg::face::PosN class.
 */

#ifndef __VCG_FACE_POS
#define __VCG_FACE_POS

#include <assert.h>

namespace vcg {
namespace face {

/** \addtogroup face */
/*@{*/

// Needed Prototypes (pos is include before topology)
template <class FaceType>
bool IsBorder(FaceType const & f,  const int j );
template <class FaceType>
bool IsManifold(FaceType const & f,  const int j );

/**  Templated over the class face, it stores a \em position over a face in a mesh.
	It contain a pointer to the current face, 
	the index of one edge and a pointer to one of the vertices of the edge.
	See also the JumpingPos in jumping_pos.h for an iterator that loops
	around the faces of a vertex without requiring the VF topology.
 */
 

template <class FaceType> 
class Pos
{
public:

	/// The vertex type
	typedef typename FaceType::VertexType VertexType;
	///The Pos type
	typedef Pos<FaceType> PosType;
	/// The scalar type
	typedef typename VertexType::ScalarType ScalarType;
	
	/// Pointer to the face of the half-edge
	typename FaceType::FaceType *f;
	/// Index of the edge
	int z;
	/// Pointer to the vertex
	VertexType *v;

	/// Default constructor
	Pos(){}
	/// Constructor which associates the half-edge element with a face, its edge and its vertex
	Pos(FaceType * const fp, int const zp, VertexType * const vp){f=fp; z=zp; v=vp;}
	Pos(FaceType * const fp, int const zp){f=fp; z=zp; v=f->V(zp);}
	Pos(FaceType * const fp, VertexType * const vp)
	{
		f = fp;
		v = vp;
		for(int i = 0; i < f->VN(); ++i)
			if (f->V(i) == v) { z = f->Prev(i); break;}
	}

	// Official Access functions functions 
   VertexType *& V(){ return v; }
	 int         & E(){ return z; }
	 FaceType   *& F(){ return f; }


// Returns the face index of the vertex inside the face.
// Note that this is DIFFERENT from using the z member that denotes the edge index inside the face. 
// It should holds that Vind != (z+1)%3   &&   Vind == z || Vind = z+2%3
   int VInd()
	 {
		 for(int i = 0; i < f->VN(); ++i) if(v==f->V(i)) return i;
		 assert(0);
		 return -1;
	 }

  
	/// Operator to compare two half-edge
	inline bool operator == ( PosType const & p ) const {
			return (f==p.f && z==p.z && v==p.v);
	} 

	/// Operator to compare two half-edge
	inline bool operator != ( PosType const & p ) const {
			return (f!=p.f || z!=p.z || v!=p.v);
	} 
	/// Operator to order half-edge; it's compare at the first the face pointers, then the index of the edge and finally the vertex pointers
	inline bool operator <= ( PosType const & p) const {
		return	(f!=p.f)?(f<f.p):
						(z!=p.z)?(z<p.z):
						(v<=p.v);
	}	

	/// Assignment operator
	inline FaceType & operator = ( const FaceType & h ){
		f=h.f;
		z=h.z;
		v=h.v;
		return *this;
		}
	/// Set to null the half-edge
	void SetNull(){
		f=0;
		v=0;
		z=-1;
	}
	/// Check if the half-edge is null
	bool IsNull() const {
		return f==0 || v==0 || z<0;
	}

	//Cambia Faccia lungo z
	// e' uguale a FlipF solo che funziona anche per non manifold.
	/// Change face via z
	void NextF()
	{
		FaceType * t = f;
		f = t->FFp(z);
		z = t->FFi(z);
	}
  
		// Paolo Cignoni 19/6/99
		// Si muove sulla faccia adiacente a f, lungo uno spigolo che
		// NON e' j, e che e' adiacente a v 
		// in questo modo si scandiscono tutte le facce incidenti in un 
		// vertice f facendo Next() finche' non si ritorna all'inizio
		// Nota che sul bordo rimbalza, cioe' se lo spigolo !=j e' di bordo
		// restituisce sempre la faccia f ma con nj che e' il nuovo spigolo di bordo 
		// vecchi parametri:     	FaceType * & f, VertexType * v, int & j

	/// It moves on the adjacent face incident to v, via a different edge that j
	void NextE()
	{
		assert( f->V(z)==v || f->V(f->Next(z))==v ); // L'edge j deve contenere v
		FlipE();
		FlipF();
		assert( f->V(z)==v || f->V(f->Next(z))==v ); 
	}
	// Cambia edge mantenendo la stessa faccia e lo stesso vertice
	/// Changes edge maintaining the same face and the same vertex
	void FlipE()
	{
		assert(f->V(f->Prev(z))!=v && (f->V(f->Next(z))==v || f->V((z+0)%f->VN())==v));
		if(f->V(f->Next(z))==v) z=f->Next(z);
		else z= f->Prev(z);
		assert(f->V(f->Prev(z))!=v && (f->V(f->Next(z))==v || f->V((z))==v));
	}

	// Cambia Faccia mantenendo lo stesso vertice e lo stesso edge
	// Vale che he.flipf.flipf= he
	// Se l'he e' di bordo he.flipf()==he
	// Si puo' usare SOLO se l'edge e' 2manifold altrimenti 
	// si deve usare nextf

	/// Changes face maintaining the same vertex and the same edge
	void FlipF()
	{
		assert( f->FFp(z)->FFp(f->FFi(z))==f );  // two manifoldness check
		assert(f->V(f->Prev(z))!=v && (f->V(f->Next(z))==v || f->V((z))==v));
		FaceType *nf=f->FFp(z);
		int nz=f->FFi(z);
		assert(nf->V(f->Prev(nz))!=v && (nf->V(f->Next(nz))==v || nf->V((nz))==v));
		f=nf;
		z=nz;
		assert(f->V(f->Prev(z))!=v && (f->V(f->Next(z))==v || f->V(z)==v));
	}

	/// Changes vertex maintaining the same face and the same edge
	void FlipV()
	{
		assert(f->V(f->Prev(z))!=v && (f->V(f->Next(z))==v || f->V(z)==v));
		
		if(f->V(f->Next(z))==v)
			v=f->V(z);
		else
			v=f->V(f->Next(z));

		assert(f->V(f->Prev(z))!=v && (f->V(f->Next(z))==v || f->V(z)==v));
	}
	
	// return the vertex that it should have if we make FlipV;
	VertexType *VFlip()
	{
		assert(f->V(f->Prev(z))!=v && (f->V(f->Next(z))==v || f->V(z)==v));
		if(f->V(f->Next(z))==v)	return f->V(z);
								else			return f->V(f->Next(z));
	}

  // return the vertex that it should have if we make FlipV;
	const VertexType *VFlip() const 
	{
		assert(f->cV(f->Prev(z))!=v && (f->cV(f->Next(z))==v || f->cV(z)==v));
		if(f->cV(f->Next(z))==v)	return f->cV(z);
								else			return f->cV(f->Next(z));
	}

  // return the face that it should have if we make FlipF;
	const FaceType *FFlip() const 
	{
		assert( f->FFp(z)->FFp(f->FFi(z))==f );
		assert(f->V(f->Prev(z))!=v && (f->V(f->Next(z))==v || f->V((z+0)%f->VN())==v));
		FaceType *nf=f->FFp(z);
		return nf;
  }


	// Trova il prossimo half-edge di bordo (nhe)
	// tale che 
	// --nhe.f adiacente per vertice a he.f
	// --nhe.v adiacente per edge di bordo a he.v 
	// l'idea e' che se he e' un half edge di bordo 
	// si puo scorrere tutto un bordo facendo 
	//
	//		hei=he;
	//		do
	//			hei.Nextb()
	//		while(hei!=he);
	
	/// Finds the next half-edge border 
	void NextB( )
	{
		assert(f->V(f->Prev(z))!=v && (f->V(f->Next(z))==v || f->V(z)==v));
		assert(f->FFp(z)==f); // f is border along j
	// Si deve cambiare faccia intorno allo stesso vertice v
	//finche' non si trova una faccia di bordo.
		do
			NextE();
      while(!IsBorder());
		
		// L'edge j e' di bordo e deve contenere v
		assert(IsBorder() &&( f->V(z)==v || f->V(f->Next(z))==v )); 
		
		FlipV();
		assert(f->V(f->Prev(z))!=v && (f->V(f->Next(z))==v || f->V(z)==v));
		assert(f->FFp(z)==f); // f is border along j
	}

	/// Checks if the half-edge is of border
	bool IsBorder()
	{
    return face::IsBorder(*f,z);
	}

  bool IsManifold()
	{
    return face::IsManifold(*f,z);
	}

	/*!
	 * Returns the number of vertices incident on the vertex pos is currently pointing to.
	 */
	int NumberOfIncidentVertices()
	{
		int  count		 = 0;
		bool on_border = false;
		CheckIncidentFaces(count, on_border);
		if(on_border) return (count/2)+1;
		else					return count;
	}

	/*!
	* Returns the number of faces incident on the vertex pos is currently pointing to.
	*/
	int NumberOfIncidentFaces()
	{
		int  count		 = 0;
		bool on_border = false;
		CheckIncidentFaces(count, on_border);
		if(on_border) return count/2;
		else					return count;
	}

	/** Function to inizialize an half-edge.
		@param fp Puntatore alla faccia
		@param zp Indice dell'edge
		@param vp Puntatore al vertice
	*/
	void Set(FaceType  * const fp, int const zp,  VertexType  * const vp)
	{
		f=fp;z=zp;v=vp;
		assert(f->V(f->Prev(z))!=v && (f->V(f->Next(z))==v || f->V(z)==v));
	}
	
	void Set(FaceType  * const pFace, VertexType  * const pVertex)
	{
		f = pFace;
		v = pVertex;
		for(int i  = 0; i < f->VN(); ++i) if(f->V(i) == v ) {z = f->Prev(i);break;}
	}

	void Assert()
	#ifdef _DEBUG
	{
		FaceType ht=*this;
		ht.FlipF();
		ht.FlipF();
		assert(ht==*this);

		ht.FlipE();
		ht.FlipE();
		assert(ht==*this);

		ht.FlipV();
		ht.FlipV();
		assert(ht==*this);
	}
	#else
	{}
	#endif


	protected:
		void CheckIncidentFaces(int & count, bool & on_border)
		{
			PosType ht = *this;
			do
			{
				++count;
				ht.NextE();
				if(ht.IsBorder()) on_border=true;
			} while (ht != *this);
		}
};

template <class FaceType>
/** Class PosN.
	This structure is equivalent to a Pos, but it contains a normal.
	@param FaceType (Template-Parameter) Specifies the type of the faces
 */ 
class PosN : public Pos<FaceType>
{
public:
	typedef typename FaceType::CoordType CoordType;
	//normale per visualizzazione creaseangle
	CoordType normal;
};


/** Class VFIterator.
	This class is used as an iterator over the VF adjacency. 
  It allow to easily traverse all the faces around a given vertex v;
  The faces are traversed in no particular order. No Manifoldness requirement.

  typical example:

    VertexPointer v;
    vcg::face::VFIterator<FaceType> vfi(v);	
    for (;!vfi.End();++vfi)
			vfi.F()->ClearV();
			
		// Alternative 

    vcg::face::VFIterator<FaceType> vfi(f, 1);
		while (!vfi.End()){
			vfi.F()->ClearV();
			++vfi;
		}


	See also the JumpingPos in jumping_pos.h for an iterator that loops
	around the faces of a vertex using FF topology and without requiring the VF topology.

 */ 

template <typename FaceType> 
class VFIterator
{
public:

	/// The vertex type
	typedef typename FaceType::VertexType VertexType;
	/// The Base face type
	typedef  FaceType  VFIFaceType;
	/// The vector type
	typedef typename VertexType::CoordType CoordType;
	/// The scalar type
	typedef typename VertexType::ScalarType ScalarType;

	/// Pointer to the face of the half-edge
	FaceType *f;
	/// Index of the vertex
	int z;
	
	/// Default constructor
	VFIterator(){}
	/// Constructor which associates the half-edge elementet with a face and its vertex
	VFIterator(FaceType * _f,  const int &  _z){f = _f; z = _z;}

	/// Constructor which takes a pointer to vertex 
	VFIterator(VertexType * _v){f = _v->VFp(); z = _v->VFi();}

	VFIFaceType *&	F() { return f;}
	int	&					  I() { return z;}
  
  // Access to the vertex. Having a VFIterator vfi, it corresponds to 
  // vfi.V() = vfi.F()->V(vfi.I())
  inline VertexType *V() const { return f->V(z);}

  inline VertexType * const & V0() const { return f->V0(z);} 
  inline VertexType * const & V1() const { return f->V1(z);}
  inline VertexType * const & V2() const { return f->V2(z);}
	
  bool End() const {return f==0;}
  VFIFaceType *operator++() {
    FaceType* t = f;
		f = f->VFp(z);
		z = t->VFi(z);
    return f;
  }
		
};

/*@}*/
}	 // end namespace
}	 // end namespace
#endif
