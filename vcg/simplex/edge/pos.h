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
Revision 1.2  2004/05/10 14:40:47  ganovelli
name of adhacency function updated

Revision 1.1  2004/05/10 14:01:56  ganovelli
created

****************************************************************************/


#ifndef __VCG_EDGE_POS
#define __VCG_EDGE_POS

namespace vcg {
namespace edge {

/*
 Vertex_Edge: run over the fan of a vertex (no order is specified)
*/
/** Class VertexStar
	@param EDGETYPE Specifies the type of the faces
 */
template <class EDGETYPE> 
class VertexStar
{
public:
	/// Pointer to an edge
	EDGETYPE *e;
	/// Local index of the vertex
	int z;
	/// Default Constructor
	VertexStar() {}
	/// Constructor which associates the EdgePos elementet with a face and its edge
	VertexStar(EDGETYPE  * const ep, int const zp)
	{
		e=ep;
		z=zp;
	}

	/// Function to jump on the next face of the list of vertex z
	void NextF()
	{
		EDGETYPE * t = e;
		e = (EDGETYPE *)t->VEp(z);
		z = t->VEi(z);
	}
};


/*
 
*/
/** Class Pos.
	This structure is equivalent to a half-edge. 
	@param MFTYPE (Template-Parameter) Specifies the type of the edges
 */
template <class EDGETYPE> 
class Pos
{
public:

	/// The vertex type
	typedef	typename EDGETYPE::VertexType VertexType;
	/////The HEdgePos type
	typedef Pos< EDGETYPE> POSTYPE;
	///// The vector type
	//typedef typename MVTYPE::coord_type vectorial_type;
	///// The scalar type
	//typedef typename MVTYPE::scalar_type scalar_type;

	/// Pointer to the face of the half-edge
	EDGETYPE *e;
	/// Pointer to the vertex
	VertexType *v;

	/// Default constructor
	Pos(){}
	/// Constructor which associates the half-edge elementet with a face, its edge and its vertex
	Pos(EDGETYPE  * const ep, int const zp,
				VertexType  * const vp){e=ep;v=vp;}

	/// Operator to compare two half-edge
	inline bool operator == ( POSTYPE const & p ) const {
			return (e==p.e &&v==p.v);
	} 

	/// Operator to compare two half-edge
	inline bool operator != ( POSTYPE const & p ) const {
			return (e!=p.e || v!=p.v);
	} 
	/// Operator to order half-edge; it's compare at the first the face pointers, then the index of the edge and finally the vertex pointers
	inline bool operator <= ( POSTYPE const & p) const {
		return	(e!=p.e)?(e<e.p):
						(v<=p.v);
	}	

	/// Assignment operator
	inline POSTYPE & operator = ( const POSTYPE & h ){
		e=h.e;
		v=h.v;
		return *this;
		}
	/// Set to null the half-edge
	void SetNull(){
		e=0;
		v=0;
	}
	/// Check if the half-edge is null
	bool IsNull() const {
		return e==0 || v==0 ;
	}

	//Cambia Faccia lungo z
	// e' uguale a FlipF solo che funziona anche per non manifold.
	/// Change face via z
	void NextE()
	{
		FlipV();
		FlipE();
	}
  
		// Paolo Cignoni 19/6/99
		// Si muove sulla faccia adiacente a f, lungo uno spigolo che
		// NON e' j, e che e' adiacente a v 
		// in questo modo si scandiscono tutte le facce incidenti in un 
		// vertice f facendo Next() finche' non si ritorna all'inizio
		// Nota che sul bordo rimbalza, cioe' se lo spigolo !=j e' di bordo
		// restituisce sempre la faccia f ma con nj che e' il nuovo spigolo di bordo 
		// vecchi parametri:     	MFTYPE * & f, MVTYPE * v, int & j

	// Cambia edge mantenendo la stessa faccia e lo stesso vertice
	/// Changes edge maintaining the same face and the same vertex
	void FlipV()
	{
		v = (e->V(0)==v)?e->V(1):e->V(0);
	}
	void FlipE()
	{
		assert( (e->V(0)==v) ||(e->V(1)==v));
		e = (e->V(0)==v)?e->EEp(0):e->EEp(1);
	}
	int Z(){
		return (e->V(0)==v)?0:1;
		}
	// return the vertex that it should have if we make FlipV;
	VertexType *VFlip()
	{
		return (e->V(0)==v)?e->V(1):e->V(0);
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
	//	assert(f->V((z+2)%3)!=v && (f->V((z+1)%3)==v || f->V((z+0)%3)==v));
	//	assert(f->F(z)==f); // f is border along j
	//// Si deve cambiare faccia intorno allo stesso vertice v
	////finche' non si trova una faccia di bordo.
	//	do
	//		NextE();
	//	while(!f->IsBorder(z));
	//	
	//	// L'edge j e' di bordo e deve contenere v
	//	assert(f->IsBorder(z) &&( f->V(z)==v || f->V((z+1)%3)==v )); 
	//	
	//	FlipV();
	//	assert(f->V((z+2)%3)!=v && (f->V((z+1)%3)==v || f->V((z+0)%3)==v));
	//	assert(f->F(z)==f); // f is border along j
	}

	/// Checks if the half-edge is of border
	//bool IsBorder()
	//{
		//return f->IsBorder(z);
	//}

	/// Return the dimension of the star
	//int StarSize()
	//{
		//int n=0;
		//POSTYPE ht=*this;
		//bool bf=false;
		//do
		//{
		//	++n;
		//	ht.NextE();
		//	if(ht.IsBorder()) bf=true;
		//} while(ht!=*this);

		//if(bf) return n/2;
		//else return n;
	//}

	/** Function to inizialize an half-edge.
		@param fp Puntatore alla faccia
		@param zp Indice dell'edge
		@param vp Puntatore al vertice
	*/
	void Set(EDGETYPE  * const ep,VertexType  * const vp)
	{
		e=ep;v=vp;
	}

	void Assert()
	#ifdef _DEBUG
	{/*
		POSTYPE ht=*this;
		ht.FlipE();
		ht.FlipE();
		assert(ht==*this);

		ht.FlipE();
		ht.FlipE();
		assert(ht==*this);

		ht.FlipV();
		ht.FlipV();
		assert(ht==*this);*/
	}
	#else
	{}
	#endif

	// Controlla la coerenza di orientamento di un hpos con la relativa faccia
	/// Checks the orientation coherence of a half-edge with the face
	//inline bool Coherent() const
	//{
	//	return v == f->V(z);	// e^(ip)+1=0 ovvero E=mc^2
	//}

};
	}	 // end namespace
}	 // end namespace
#endif
