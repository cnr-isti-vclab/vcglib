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



****************************************************************************/

#ifndef __VCG_TETRA_POS
#define __VCG_TETRA_POS

namespace vcg {
namespace tetra {

 /** \addtogroup tetra */
/*@{*/


  /**  Class VTIterator.
 This is a vertex - tetrahedron iterator
		@param MTTYPE (Template Parameter) Specifies the type of the tetrahedron.
 */
template < class MTTYPE> 
class VTIterator
{
public:
	/// The tetrahedron type 
	typedef typename MTTYPE TetraType;
private:
	/// Pointer to a tetrahedron
	TetraType *_vt;
	/// Index of one vertex
	int _vi;
	/// Default Constructor
public:
	VTIterator() {}
	/// Constructor which associates the EdgePos elementet with a face and its edge
	VTIterator(TetraType  * const tp, int const zp)
	{	
		_vt=tp;
		_vi=zp;
	}

	~VTIterator(){};

  	/// Return the tetrahedron stored in the half edge
	inline TetraType & Vt()
	{
		return _vt;
	} 

  	/// Return the tetrahedron stored in the half edge
	inline const TetraType & Vt() const
	{
		return _vt;
	}

  	/// Return the index of vertex as seen from the tetrahedron
	inline int & Vi()
	{
		return _vi;
	}

  	/// Return the index of vertex as seen from the tetrahedron
	inline const int & Vi() const
	{
		return _vi;
	}

	/// move on the next tetrahedron that share the vertex
	bool NextT()
	{	
		int vi=Vi();
		TetraType * tw = Vt();
		Vt() = tw->TVp[vi];
		Vi() = tw->TVi[vi];
		assert(((tw->V(vi))==(Vt()->V(Vi())))||(t==NULL));
    return (Vt()!=NULL);
  }
};

/** \addtogroup tetra */
/*@{*/

/**  Templated over the class tetrahedron, it stores a \em position over a tetrahedron in a mesh.
	It contain a pointer to the current tetrahedron, 
	the index of one face,edge and a edge's incident vertex.
 */
template < class MTTYPE> 
class Pos
{
public:

  /// The tetrahedron type
	typedef	typename MTTYPE TetraType;
	/// The vertex type
	typedef	typename TetraType::MVTYPE VertexType;
  /// The coordinate type
	typedef	typename TetraType::MVTYPE::CoordType CoordType;
  ///The HEdgePos type
	typedef Pos<TetraType> BasePosType;

private:
	/// Pointer to the tetrahedron of the half-edge
	TetraType *_t;
	/// Index of the face
	char _f;
	/// Index of the edge
	char _e;
	/// Pointer to the vertex
	char _v;

public:
	/// Default constructor
	Pos(){SetNull();};
	/// Constructor which associates the half-edge elementet with a face, its edge and its vertex
	Pos(TetraType * const tp, char const fap,char const ep,
		char const const vp){_t=tp;_f=fap;_e=ep;_v=vp;}

	~Pos(){};
  
  	/// Return the tetrahedron stored in the half edge
	inline TetraType & T()
	{
		return _t;
	} 

  	/// Return the tetrahedron stored in the half edge
	inline const TetraType & T() const
	{
		return _t;
	}

  	/// Return the index of face as seen from the tetrahedron
	inline char & F()
	{
		return _f;
	}

  	/// Return the index of face as seen from the tetrahedron
	inline const char & F() const
	{
		return _f;
	}

	/// Return the index of face as seen from the tetrahedron
	inline char & E()
	{
		return _e;
	}

  	/// Return the index of face as seen from the tetrahedron
	inline const char & E() const
	{
		return _e;
	}

/// Return the index of vertex as seen from the tetrahedron
	inline char & V()
	{
		return _v;
	}

  	/// Return the index of vertex as seen from the tetrahedron
	inline const char & V() const
	{
		return _v;
	}

	/// Operator to compare two half-edge
	inline bool operator == ( BasePosType const & p ) const {
			return (T()==p.T() && F()==p.F() && E==p.E() && V==p.V());
	} 

	/// Operator to compare two half-edge
	inline bool operator != ( BasePosType const & p ) const {
			return (!((*this)==p));
	} 

	/// Assignment operator
	inline BasePosType & operator = ( const BasePosType & h ){
		T()=h.T();
		F()=h.F();
		E()=h.E();
		V()=h.V();
		return *this;
		}

	/// Set to null the half-edge
	void SetNull(){
		T()=0;
		F()=-1;
		E()=-1;
		V()=-1;
	}

	/// Check if the half-edge is null
	bool IsNull() const {
		return ((T()==0) || (F()<0) || (E()<0) || (V()<0));
	}


	/// Changes edge maintaining the same face and the same vertex
	void FlipE()
	{
		
		//take the absolute index of the tree edges of the faces
    char e0=vcg::Tetra::EofF(fa,0);
		char e1=vcg::Tetra::EofF(fa,1);
		char e2=vcg::Tetra::EofF(fa,2);
		//eliminate the same as himself
		if (e0==E())
			{
			 e0=e1;
			 e1=e2;
			}
		else
		if (e1==E())
			{
			 e1=e2;
			}

		//now choose the one that preserve the same vertex
     if ((vcg::Tetra::VofE(e1,0)==V())||(vcg::Tetra::VofE(e1,1)==V()))
			E()=e1;
		else
			E()=e0;
	}


	/// Changes vertex maintaining the same face and the same edge
	void FlipV()
	{
		// in the same edge choose the one that change
		char v0=vcg::Tetra::VofE(E(),0);
		char v1=vcg::Tetra::VofE(E(),1);
		if (v0!=V())
			V()=v0;
		else
			V()=v1;
	}

	/// Changes face maintaining the same vertex and the same edge
	void FlipF()
	{
    char f0=vcg::Tetra::FofE(z,0);
		char f1=vcg::Tetra::FofE(z,1);
		if (f0!=F())
			F()=f0;
		else
			F()=f1;
	}

	/// Changes tetrahedron maintaining the same face edge and vertex'... to finish
	void FlipT()
	{
		
		//save the two vertices of the old edge
		char *v0=vcg::Tetra::VofE(z,0);
		char *v1=vcg::Tetra::VofE(z,1);

		//get new tetrahedron according to faceto face topology
		TetraType *nt=T()->TTp(F());
		char nfa=T()->TTi(F());
		if (nfa!=-1)
		{
			//find the right edge
      char ne0=vcg::Tetra::EofF(nfa,0);
			char ne1=vcg::Tetra::EofF(nfa,1);
			char ne2=vcg::Tetra::EofF(nfa,2);

			//verify that the two vertices of tetrahedron are identical
			if (((nt->VE(ne0,0)==v0)&&(nt->VE(ne0,1)==v1))||
				((nt->VE(ne0,1)==v0)&&(nt->VE(ne0,0)==v1)))
			z=ne0;
			else
			if (((nt->VE(ne1,0)==v0)&&(nt->VE(ne1,1)==v1))||
				((nt->VE(ne1,1)==v0)&&(nt->VE(ne1,0)==v1)))	
			z=ne1;
			else
			z=ne2;
			t=nt;
			fa=nfa;
		}
			}

	
	
	
	void NextE( )
	{
		//assert(t->V((z+2)%4)!=v && (t->V((z+1)%4)==v || t->V((z+0)%4)==v));
		#ifdef _DEBUG
		vertex_type *v0old=t->VE(z,0);
		vertex_type *v1old=t->VE(z,1);
		#endif
		FlipT();
		FlipF();
		#ifdef _DEBUG
		vertex_type *v0=t->VE(z,0);
		vertex_type *v1=t->VE(z,1);
		assert(v1!=v0);
		assert(((v0==v0old)&&(v1==v1old))||((v1==v0old)&&(v0==v1old)));
		#endif
		
	}
	
	void NextV( )
	{
		//assert(t->V((z+2)%4)!=v && (t->V((z+1)%4)==v || t->V((z+0)%4)==v));
		int j;
		int indexv;

		// find the index of the current vertex
		for (j=0;j<4;j++)
			{
			if (v==t->V(j))
					indexv=j;
			}
		//increase the iterator
		EdgePosT <MTTYPE> e(t,indexv);
		e.NextT();
		t=e.t;

		//assert(t->V((z+2)%4)!=v && (t->V((z+1)%4)==v || t->V((z+0)%4)==v));
	}

	void NextF( )
	{
		assert(t->V((z+2)%4)!=v && (t->V((z+1)%4)==v || t->V((z+0)%4)==v));
		FlipT();
		assert(t->V((z+2)%4)!=v && (t->V((z+1)%4)==v || t->V((z+0)%4)==v));
	}

	void NextT( )
	{
		assert(t->V((z+2)%4)!=v && (t->V((z+1)%4)==v || t->V((z+0)%4)==v));
		/*fa=(fa+1)%4;
		t=T(*/
		assert(t->V((z+2)%4)!=v && (t->V((z+1)%4)==v || t->V((z+0)%4)==v));
	}

	
	/** Function to inizialize an half-edge.
		@param fp Puntatore alla faccia
		@param zp Indice dell'edge
		@param vp Puntatore al vertice
	*/
	void Set(MTTYPE * const tp, int const fap,int const zp,vertex_type  * const vp)
	{	t=tp;fa=fap;z=zp;v=vp;
		assert(t->V((z+2)%4)!=v && (t->V((z+1)%4)==v || t->V((z+0)%4)==v));
	}

	void Assert()
	#ifdef _DEBUG
	{	
		HETYPE ht=*this;
		ht.FlipT();
		ht.FlipT();
		assert(ht==*this);

		ht=*this;
		ht.FlipF();
		ht.FlipF();
		assert(ht==*this);

		ht=*this;
		ht.FlipE();
		ht.FlipE();
		assert(ht==*this);

		ht=*this;
		ht.FlipV();
		ht.FlipV();
		assert(ht==*this);
	}
	#else
	{}
	#endif

	/*// Controlla la coerenza di orientamento di un hpos con la relativa faccia
	/// Checks the orientation coherence of a half-edge with the face
	inline bool Coerent() const
	{
		return v == t->V(z);	// e^(ip)+1=0 ovvero E=mc^2
	}*/

};

template < class MTTYPE> 
class HEdgePosTEdge:public HEdgePosT<MTTYPE>
{
	public :
	MTTYPE *t_initial;
	short int fa_initial;
	short int back;

	/// Constructor which associates the half-edge elementet with a face, its edge and its vertex
	HEdgePosTEdge(){}

	HEdgePosTEdge(MTTYPE * const tp,const int  fap,const int  zp,
		vertex_type  *  vp){t=tp;fa=fap;fa_initial=fap;z=zp;v=vp;t_initial=tp;back=0;}

	void NextE()
	{	
#ifdef _DEBUG
		int cont=0;
#endif
		MTTYPE *tpred=t;
		HEdgePosT<MTTYPE>::NextE();
		//rimbalzo
		if (tpred==t)
		{
			while (t!=t_initial)
			{
				HEdgePosT<MTTYPE>::NextE();
				#ifdef _DEBUG
				 cont++;
				 assert (cont<500);
				#endif
			}
			back++;
			if (back==1)
			{
				HEdgePosT<MTTYPE>::NextE();
			}
			
		}
		
	}

//change tetrahedron endreturn the number of the face to put on the fan
int NextFaceOnFan()
	{
	HEdgePosTEdge::NextE();
	//get the faces that are not on the edge 
	int fa0=t->FE(z,0);
	int fa1=t->FE(z,1);
	//they are the 2 faces that remain
	int fa2=(fa0+1)%4;
	while ((fa2==fa0)||(fa2==fa1))
	{
		fa2=(fa2+1)%4;
	}
	int fa3=(fa2+1)%4;
	while ((fa3==fa0)||(fa3==fa1)||(fa3==fa2))
	{
		fa3=(fa3+1)%4;
	}
	bool first=false;
	for (int i=0;i<3;i++)
		if (t->FV(fa2,i)==v)
			first=true;
	if (first)
		return fa2;
	else 
		return fa3;
	}
};


#endif