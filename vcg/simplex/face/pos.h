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

/** \file pos.h
 * Definition of vcg:face::Pos class.
 * This file contain the definition of vcg::face::Pos class and the derived vcg::face::PosN class.
 */

#ifndef __VCG_FACE_POS
#define __VCG_FACE_POS

namespace vcg {
namespace face {

/** \addtogroup face */
/*@{*/

/**  Templated over the class face, it stores a \em position over a face in a mesh.
	It contain a pointer to the current face, 
	the index of one edge and a edge's incident vertex.
 */
template <class FaceType> 
class Pos
{
public:

	/// The vertex type
	typedef typename FaceType::VertexType VertexType;
	///The HEdgePos type
	typedef Pos<FaceType> PosType;
	/// The vector type
	typedef typename VertexType::CoordType CoordType;
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
	/// Constructor which associates the half-edge elementet with a face, its edge and its vertex
	Pos(FaceType  * const fp, int const zp,	VertexType  * const vp){f=fp; z=zp; v=vp;}
	Pos(FaceType  * const fp, int const zp){f=fp; z=zp; v=f->V(zp);}

	// access functions
	VertexType *& V(){return f->UberV(z);}
	VertexType *& V(const int & i){assert( (i>=0) && (i<2)); return f->UberV( (z +i) %3);}

	/// Operator to compare two half-edge
	inline bool operator == ( FaceType const & p ) const {
			return (f==p.f && z==p.z && v==p.v);
	} 

	/// Operator to compare two half-edge
	inline bool operator != ( FaceType const & p ) const {
			return (f!=p.f || z!=p.z || v!=p.v);
	} 
	/// Operator to order half-edge; it's compare at the first the face pointers, then the index of the edge and finally the vertex pointers
	inline bool operator <= ( FaceType const & p) const {
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
		assert( f->V(z)==v || f->V((z+1)%3)==v ); // L'edge j deve contenere v
		FlipE();
		FlipF();
		assert( f->V(z)==v || f->V((z+1)%3)==v ); 
	}
	// Cambia edge mantenendo la stessa faccia e lo stesso vertice
	/// Changes edge maintaining the same face and the same vertex
	void FlipE()
	{
		assert(f->V((z+2)%3)!=v && (f->V((z+1)%3)==v || f->V((z+0)%3)==v));
		if(f->V((z+1)%3)==v) z=(z+1)%3;
		else z=(z-1+3)%3;
		assert(f->V((z+2)%3)!=v && (f->V((z+1)%3)==v || f->V((z+0)%3)==v));
	}

	// Cambia Faccia mantenendo lo stesso vertice e lo stesso edge
	// Vale che he.flipf.flipf= he
	// Se l'he e' di bordo he.flipf()==he
	// Si puo' usare SOLO se l'edge e' 2manifold altrimenti 
	// si deve usare nextf

	/// Changes face maintaining the same vertex and the same edge
	void FlipF()
	{
		assert( f->FFp(z)->FFp(f->FFi(z))==f );
		assert(f->V((z+2)%3)!=v && (f->V((z+1)%3)==v || f->V((z+0)%3)==v));
		FaceType *nf=f->FFp(z);
		int nz=f->FFi(z);
		assert(nf->V((nz+2)%3)!=v && (nf->V((nz+1)%3)==v || nf->V((nz+0)%3)==v));
		f=nf;
		z=nz;
		assert(f->V((z+2)%3)!=v && (f->V((z+1)%3)==v || f->V((z+0)%3)==v));
	}

	/// Changes vertex maintaining the same face and the same edge
	void FlipV()
	{
		assert(f->V((z+2)%3)!=v && (f->V((z+1)%3)==v || f->V((z+0)%3)==v));
		
		if(f->V((z+1)%3)==v)
			v=f->V((z+0)%3);
		else
			v=f->V((z+1)%3);

		assert(f->V((z+2)%3)!=v && (f->V((z+1)%3)==v || f->V((z+0)%3)==v));
	}
	
	// return the vertex that it should have if we make FlipV;
	VertexType *VFlip()
	{
		assert(f->V((z+2)%3)!=v && (f->V((z+1)%3)==v || f->V((z+0)%3)==v));
		if(f->V((z+1)%3)==v)	return f->V((z+0)%3);
								else			return f->V((z+1)%3);
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
		assert(f->V((z+2)%3)!=v && (f->V((z+1)%3)==v || f->V((z+0)%3)==v));
		assert(f->FFp(z)==f); // f is border along j
	// Si deve cambiare faccia intorno allo stesso vertice v
	//finche' non si trova una faccia di bordo.
		do
			NextE();
		while(!f->IsBorder(z));
		
		// L'edge j e' di bordo e deve contenere v
		assert(f->IsBorder(z) &&( f->V(z)==v || f->V((z+1)%3)==v )); 
		
		FlipV();
		assert(f->V((z+2)%3)!=v && (f->V((z+1)%3)==v || f->V((z+0)%3)==v));
		assert(f->FFp(z)==f); // f is border along j
	}

	/// Checks if the half-edge is of border
	bool IsBorder()
	{
		return f->IsBorder(z);
	}

	/// Return the dimension of the star
	int StarSize()
	{
		int n=0;
		FaceType ht=*this;
		bool bf=false;
		do
		{
			++n;
			ht.NextE();
			if(ht.IsBorder()) bf=true;
		} while(ht!=*this);

		if(bf) return n/2;
		else return n;
	}

	/** Function to inizialize an half-edge.
		@param fp Puntatore alla faccia
		@param zp Indice dell'edge
		@param vp Puntatore al vertice
	*/
	void Set(FaceType  * const fp, int const zp,  VertexType  * const vp)
	{
		f=fp;z=zp;v=vp;
		assert(f->V((z+2)%3)!=v && (f->V((z+1)%3)==v || f->V((z+0)%3)==v));
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

	// Controlla la coerenza di orientamento di un hpos con la relativa faccia
	/// Checks the orientation coherence of a half-edge with the face
	inline bool Coerent() const
	{
		return v == f->V(z);	// e^(ip)+1=0 ovvero E=mc^2
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
	int	&					I() { return z;}

  bool End() const {return f==0;}
  VFIFaceType *operator++() {
    FaceType* t = f;
		f = t->VFp(z);
		z = t->VFi(z);
    return f;
  }
		
};

/*@}*/
}	 // end namespace
}	 // end namespace
#endif
