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
Revision 1.7  2004/05/12 12:23:23  cignoni
Conformed C++ syntax to GCC requirements

Revision 1.6  2004/05/11 16:03:18  ganovelli
changed from "thi" to "&f" in Vfdetach

Revision 1.5  2004/05/10 15:20:49  cignoni
Updated names of POS and adj functions to the new standards for many functions

Revision 1.4  2004/03/18 16:00:10  cignoni
minor changes


****************************************************************************/

#ifndef _VCG_FACE_TOPOLOGY
#define _VCG_FACE_TOPOLOGY
#include <vcg/simplex/face/pos.h>
namespace vcg {
namespace face {
/** \addtogroup face */
/*@{*/

/** Return a boolean that indicate if the face is complex.
    @param j Index of the edge
	@return true se la faccia e' manifold, false altrimenti
*/
template <class FaceType>
inline bool IsManifold( FaceType const & f, const int j ) 
{
  if(FaceType::HasFFAdjacency())
	  return ( f.FFp(j) == &f || &f == f.FFp(j)->FFp(f.FFi(j)) );
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
	  return f.FFp(j) == &f ;
  else 
    return true;
}


/// Count border edges of the face
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


/// Counts the number of incident faces in a complex edge
template <class FaceType>
inline int ComplexSize(FaceType const & f, const int e)
{
if(FaceType::HasFFAdjacency())
{
  Pos< FaceType > fpos(&f,e); 
  int cnt=0;
  do
  {
		fpos.NextF();
		++cnt;
	}
	while(fpos.f=&f);
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
	assert(!IsBorder(e));
	Pos< FaceType > EPB(&f,e);  // la faccia dall'altra parte
	EPB.NextF();
	int cnt=0;
	while ( EPB.f->FFp(EPB.z) != &f)
	{ 
		assert(!IsManifold(e));   // Si entra in questo loop solo se siamo in una situazione non manifold.
		assert(!EPB.f->IsBorder(EPB.z));
		EPB.NextF();
		cnt++;
	}
	assert(EPB.f->FFp(EPB.z)==&f);

	EPB.f->FFp(EPB.z) = F(e);
	EPB.f->FFi(EPB.z) = Z(e);
	
	F(e) = &f;
	Z(e) = e;

	EPB.f->SetM();
	f.SetM();
}


/** This function attach the face (via the edge z1) to another face (via the edge z2). It's possible to use it also in non-two manifold situation.
		The function cannot be applicated if the adjacencies among faces aren't define.
		@param z1 Index of the edge
		@param f2 Pointer to the face
		@param z2 The edge of the face f2 
*/
template <class FaceType>
void Attach(FaceType * &f, int z1, FaceType *&f2, int z2)
{
	//typedef FEdgePosB< FACE_TYPE > ETYPE;
	Pos< FaceType > EPB(f2,z2);
	Pos< FaceType > TEPB;
	TEPB = EPB;
	EPB.NextF();
	while( EPB.f != f2)  //Alla fine del ciclo TEPB contiene la faccia che precede f2
	{
		TEPB = EPB;
		EPB.NextF();
	}
	//Salvo i dati di f1 prima di sovrascrivere
	FaceType *f1prec = f.FFp(z1);  
	int z1prec = f.FFi(z1);
	//Aggiorno f1
	f->FFp(z1) = TEPB.f->FFp(TEPB.z);  
	f->FFi(z1) = TEPB.f->FFi(TEPB.z);
	//Aggiorno la faccia che precede f2
	TEPB.f->FFp(TEPB.z) = f1prec;
	TEPB.f->FFi(TEPB.z) = z1prec;
}


template <class FaceType>
void AssertAdj(FaceType & f)
{
	assert(f.FFp(0)->FFp(f.FFi(0))==&f);
	assert(f.FFp(1)->FFp(f.FFi(1))==&f);
	assert(f.FFp(2)->FFp(f.FFi(2))==&f);

	assert(f.FFp(0)->FFi(f.FFi(0))==0);
	assert(f.FFp(1)->FFi(f.FFi(1))==1);
	assert(f.FFp(2)->FFi(f.FFi(2))==2); 
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
  SwapFaceType *tmp, *prec;
  int t, precz;

  swap ( f.V((z  )%3),f.V((z+1)%3));

  if(f.HasFFAdjacency() )
  {
    // TODO!!!
  }
}



// Stacca la faccia corrente dalla catena di facce incidenti sul vertice z, 
// NOTA funziona SOLO per la topologia VF!!!
// usata nelle classi di collapse
template <class FaceType>
void VFDetach(FaceType & f, int z)
{
	printf("detach %d \n",&f);
	if(f.V(z)->VFp()==&f )  //if it is the first face detach from the begin
	{
		int fz = f.V(z)->VFi();
		f.V(z)->VFp() = f.VFp(fz);
		f.V(z)->VFi() = f.VFi(fz);
	}
	else  // scan the list of faces in order to finde the current face f to be detached
	{
    VFIterator<FaceType> x(f.V(z)->VFp(),f.V(z)->VFi());
    VFIterator<FaceType> y;

		for(;;)
		{
			y = x;
			x++;
			assert(x.f!=0);
			if(x.f==&f) // found!
			{
				y.f->VFp(y.z) = f.VFp(z);
				y.f->VFi(y.z) = f.VFi(z);
				break;
			}
		}
	}
}

/*@}*/
}	 // end namespace
}	 // end namespace

#endif

