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

$LOG$

****************************************************************************/

#ifndef _VCG_FACE_TOPOLOGY
#define _VCG_FACE_TOPOLOGY

namespace vcg {
  namespace face {
/*#*******************	
*  Adjacency Members *
**********************/


/** Return a boolean that indicate if the face is complex.
    @param j Index of the edge
	@return true se la faccia e' manifold, false altrimenti
*/
template <class FaceType>
inline bool IsManifold( FaceType const & f, const int j ) 
{
  if(FaceType::HasFFAdjacency())
	  return ( f.F(j) == &f || &f == f.F(j)->F(f.Z(j)) );
  else 
    return true;
}

/** Return a boolean that indicate if the j-th edge of the face is a border.
	@param j Index of the edge
	@return true if j is an edge of border, false otherwise
*/
template <class FaceType>
inline bool IsBorder(FaceType const & f,  const int j ) 
{
  if(FaceType::HasFFAdjacency())
	  return f.F(j) == &f ;
  else 
    return true;
}


/// This function counts the boreders of the face
template <class FaceType>
inline int BorderCount(FaceType const & f) 
{
  if(FaceType::HasFFAdjacency())
  {
    int t = 0;
	  if( f.IsBorder(0) ) ++t;
	  if( f.IsBorder(1) ) ++t;
	  if( f.IsBorder(2) ) ++t;
	  return t;
  }
	else 	return 3;
}


/// This function counts the number of incident faces in a complex edge
template <class FaceType>
inline int ComplexSize(FaceType const & f, const int e)
{
if(FaceType::HasFFAdjacency())
{
  int cnt=0;
	FACE_TYPE *fi=(FACE_TYPE *)this;
	int nzi,zi=e;
	do
	{
		nzi=fi->Z(zi);
		fi=fi->F(zi);
		zi=nzi;
		++cnt;
	}
	while(fi!=this);
	return cnt;
}
  assert(0);
	return 2;
}

/*Funzione di detach che scollega una faccia da un ciclo 
(eventualmente costituito da due soli elementi) incidente su un edge*/
/** This function detach the face from the adjacent face via the edge e. It's possible to use it also in non-two manifold situation.
		The function cannot be applicated if the adjacencies among faces aren't define.
		@param e Index of the edge
*/
template <class FaceType>
void Detach(FaceType & f, const int e)
{
	typedef FEdgePosB< FACE_TYPE > ETYPE;

	assert(!IsBorder(e));
	ETYPE EPB(this,e);  // la faccia dall'altra parte
	EPB.NextF();
	int cnt=0;
	while ( EPB.f->F(EPB.z) != this)
	{ 
		assert(!IsManifold(e));   // Si entra in questo loop solo se siamo in una situazione non manifold.
		assert(!EPB.f->IsBorder(EPB.z));
		EPB.NextF();
		cnt++;
	}
	assert(EPB.f->F(EPB.z)==this);

	EPB.f->F(EPB.z) = F(e);
	EPB.f->Z(EPB.z) = Z(e);
	
	F(e) = this;
	Z(e) = e;

	EPB.f->SetM();
	this->SetM();
}


/** This function attach the face (via the edge z1) to another face (via the edge z2). It's possible to use it also in non-two manifold situation.
		The function cannot be applicated if the adjacencies among faces aren't define.
		@param z1 Index of the edge
		@param f2 Pointer to the face
		@param z2 The edge of the face f2 
*/
template <class FaceType>
void Attach(int z1, FaceType *&f2, int z2)
{
	typedef FEdgePosB< FACE_TYPE > ETYPE;
	ETYPE EPB(f2,z2);
	ETYPE TEPB;
	TEPB = EPB;
	EPB.NextF();
	while( EPB.f != f2)  //Alla fine del ciclo TEPB contiene la faccia che precede f2
	{
		TEPB = EPB;
		EPB.NextF();
	}
	//Salvo i dati di f1 prima di sovrascrivere
	face_base *f1prec = this->F(z1);  
	int z1prec = this->Z(z1);
	//Aggiorno f1
	this->F(z1) = TEPB.f->F(TEPB.z);  
	this->Z(z1) = TEPB.f->Z(TEPB.z);
	//Aggiorno la faccia che precede f2
	TEPB.f->F(TEPB.z) = f1prec;
	TEPB.f->Z(TEPB.z) = z1prec;
}


template <class FaceType>
void AssertAdj()
{
	assert(F(0)->F(Z(0))==this);
	assert(F(1)->F(Z(1))==this);
	assert(F(2)->F(Z(2))==this);

	assert(F(0)->Z(Z(0))==0);
	assert(F(1)->Z(Z(1))==1);
	assert(F(2)->Z(Z(2))==2); 
}
// Funzione di supporto usata da swap?
//template <class FaceType>
//inline void Nexts(  *&f, int &z )
//{
//    int t;
//    t = z;
//    z = (*f).Z(z);
//    f = (*f).F(t);
//}

/** This function change the orientation of the face. Inverting the index of two vertex 
@param z Index of the edge
*/
template <class SwapFaceType>
void Swap (SwapFaceType &f, const int z )
{
  int i;
  face_base *tmp, *prec;
  int t, precz;

  swap ( f.V((z  )%3),f.V((z+1)%3));

  if(f.HasFFAdjacency() )
  {
	swap ( f.F((z+1)%3),f.F((z+2)%3));
	swap ( f.Z((z+1)%3),f.Z((z+2)%3));

	for(i = 1; i < 3; i++)
	{

      tmp = this;
      t = (z+i)%3;
      do {
					prec = tmp;
					precz = t;
					Nexts(tmp,t);
      }
      while (tmp != this);
  
      (*prec).Z(precz) = (z+i)%3;
    }
  }
}



// Stacca la faccia corrente dalla catena di facce incidenti sul vertice z, 
// NOTA funziona SOLO per la topologia VF!!!
// usata nelle classi di collapse
template <class FaceType>
void VFDetach(FaceType & f, int z)
{
	if(f.V(z)->Fp()==this )
	{
		int fz = f.V(z)->Zp();
		f.V(z)->Fp() = (face_from_vert_type *) f.F(fz);
		f.V(z)->Zp() = f.Z(fz);
	}
	else
	{
			VEdgePosB<FACE_TYPE> x,y;

		x.f = V(z)->Fp();
		x.z = V(z)->Zp();

		for(;;)
		{
			y = x;
			x.NextF();
			assert(x.f!=0);
			if(x.f==this)
			{
				y.f->F(y.z) = f.F(z);
				y.f->Z(y.z) = f.Z(z);
				break;
			}
		}
	}
}


}	 // end namespace
}	 // end namespace

#endif

