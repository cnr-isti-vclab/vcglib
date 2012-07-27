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

#ifndef _VCG_FACE_TOPOLOGY
#define _VCG_FACE_TOPOLOGY

#include <vcg/simplex/face/pos.h>

#include <vector>
#include <algorithm>

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
  assert(f.cFFp(j) != 0); // never try to use this on uncomputed topology
  if(FaceType::HasFFAdjacency())
	  return ( f.cFFp(j) == &f || &f == f.cFFp(j)->cFFp(f.cFFi(j)) );
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
	  return f.cFFp(j)==&f;
    //return f.IsBorder(j);
  
  assert(0);
  return true;
}


/// Count border edges of the face
template <class FaceType>
inline int BorderCount(FaceType const & f) 
{
  if(FaceType::HasFFAdjacency())
  {
    int t = 0;
	  if( IsBorder(f,0) ) ++t;
	  if( IsBorder(f,1) ) ++t;
	  if( IsBorder(f,2) ) ++t;
	  return t;
  }
	else 	return 3;
}


/// Counts the number of incident faces in a complex edge
template <class FaceType>
inline int ComplexSize(FaceType & f, const int e)
{
  if(FaceType::HasFFAdjacency())
  {
    if(face::IsBorder<FaceType>(f,e))  return 1;
    if(face::IsManifold<FaceType>(f,e)) return 2;
                      
    // Non manifold case
    Pos< FaceType > fpos(&f,e); 
    int cnt=0;
    do
    {
		  fpos.NextF();
      assert(!fpos.IsBorder());
      assert(!fpos.IsManifold());
		  ++cnt;
	  }
	  while(fpos.f!=&f);
    assert (cnt>2);
	  return cnt;
  }
  assert(0);
	return 2;
}


/** This function check the FF topology correctness for an edge of a face. 
    It's possible to use it also in non-two manifold situation.
		The function cannot be applicated if the adjacencies among faces aren't defined.
		@param f the face to be checked 
		@param e Index of the edge to be checked 
*/
template <class FaceType>
bool FFCorrectness(FaceType & f, const int e)
{
  if(f.FFp(e)==0) return false;   // Not computed or inconsistent topology

  if(f.FFp(e)==&f) // Border
  {
   if(f.FFi(e)==e) return true;
   else return false;
  }

  if(f.FFp(e)->FFp(f.FFi(e))==&f) // plain two manifold 
  {
    if(f.FFp(e)->FFi(f.FFi(e))==e) return true;
    else return false;
  }

  // Non Manifold Case
  // all the faces must be connected in a loop.

  Pos< FaceType > curFace(&f,e);  // Build the half edge
  	int cnt=0;
  do
	{ 
		if(curFace.IsManifold()) return false;  
		if(curFace.IsBorder()) return false;
		curFace.NextF();
		cnt++;
    assert(cnt<100);
	}
  while ( curFace.f != &f);
  return true;
}


/** This function detach the face from the adjacent face via the edge e.
    It's possible to use  this function it ONLY in non-two manifold situation.
        The function cannot be applicated if the adjacencies among faces aren't defined.
        @param f the face to be detached
        @param e Index of the edge to be detached
*/
template <class FaceType>
void FFDetachManifold(FaceType & f, const int e)
{
    assert(FFCorrectness<FaceType>(f,e));
    assert(!IsBorder<FaceType>(f,e));  // Never try to detach a border edge!
    FaceType *ffp = f.FFp(e);
    //int ffi=f.FFp(e);
	int ffi=f.FFi(e);

    f.FFp(e)=&f;
    f.FFi(e)=e;
    ffp->FFp(ffi)=ffp;
    ffp->FFi(ffi)=ffi;

    f.SetB(e);
    f.ClearF(e);
    ffp->SetB(ffi);
    ffp->ClearF(ffi);

    assert(FFCorrectness<FaceType>(f,e));
    assert(FFCorrectness<FaceType>(*ffp,ffi));
}

/** This function detach the face from the adjacent face via the edge e. 
    It's possible to use it also in non-two manifold situation.
		The function cannot be applicated if the adjacencies among faces aren't defined.
		@param f the face to be detached 
		@param e Index of the edge to be detached 
*/

template <class FaceType>
void FFDetach(FaceType & f, const int e)
{
    assert(FFCorrectness<FaceType>(f,e));
    assert(!IsBorder<FaceType>(f,e));  // Never try to detach a border edge!
    int complexity;
    assert(complexity=ComplexSize(f,e));

    Pos< FaceType > FirstFace(&f,e);  // Build the half edge
	Pos< FaceType > LastFace(&f,e);  // Build the half edge
	FirstFace.NextF(); 
    LastFace.NextF();
	int cnt=0;

    // then in case of non manifold face continue to advance LastFace
    // until I find it become the one that
    // preceed the face I want to erase

	while ( LastFace.f->FFp(LastFace.z) != &f)
	{ 
        assert(ComplexSize(*LastFace.f,LastFace.z)==complexity);
		assert(!LastFace.IsManifold());   // We enter in this loop only if we are on a non manifold edge
		assert(!LastFace.IsBorder());
		LastFace.NextF();
		cnt++;
        assert(cnt<100);
	}

	assert(LastFace.f->FFp(LastFace.z)==&f);
    assert(f.FFp(e)== FirstFace.f);

	// Now we link the last one to the first one, skipping the face to be detached;
    LastFace.f->FFp(LastFace.z) = FirstFace.f;
    LastFace.f->FFi(LastFace.z) = FirstFace.z;
    assert(ComplexSize(*LastFace.f,LastFace.z)==complexity-1);

    // At the end selfconnect the chosen edge to make a border.
    f.FFp(e) = &f;
	f.FFi(e) = e;
    assert(ComplexSize(f,e)==1);

    assert(FFCorrectness<FaceType>(*LastFace.f,LastFace.z));
    assert(FFCorrectness<FaceType>(f,e));
}


/** This function attach the face (via the edge z1) to another face (via the edge z2). It's possible to use it also in non-two manifold situation.
		The function cannot be applicated if the adjacencies among faces aren't define.
		@param z1 Index of the edge
		@param f2 Pointer to the face
		@param z2 The edge of the face f2 
*/
template <class FaceType>
void FFAttach(FaceType * &f, int z1, FaceType *&f2, int z2)
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
  FaceType *f1prec = f->FFp(z1);
  int z1prec = f->FFi(z1);
	//Aggiorno f1
	f->FFp(z1) = TEPB.f->FFp(TEPB.z);  
	f->FFi(z1) = TEPB.f->FFi(TEPB.z);
	//Aggiorno la faccia che precede f2
	TEPB.f->FFp(TEPB.z) = f1prec;
	TEPB.f->FFi(TEPB.z) = z1prec;
}

/** This function attach the face (via the edge z1) to another face (via the edge z2).
		It is not possible to use it also in non-two manifold situation.
		The function cannot be applicated if the adjacencies among faces aren't define.
		@param z1 Index of the edge
		@param f2 Pointer to the face
		@param z2 The edge of the face f2
*/
template <class FaceType>
void FFAttachManifold(FaceType * &f1, int z1, FaceType *&f2, int z2)
{
  assert(IsBorder<FaceType>(*f1,z1));
  assert(IsBorder<FaceType>(*f2,z2));
  assert(f1->V0(z1) == f2->V0(z2) || f1->V0(z1) == f2->V1(z2));
  assert(f1->V1(z1) == f2->V0(z2) || f1->V1(z1) == f2->V1(z2));
  f1->FFp(z1) = f2;
  f1->FFi(z1) = z2;
  f2->FFp(z2) = f1;
  f2->FFi(z2) = z1;
}

// This one should be called only on uniitialized faces.
template <class FaceType>
void FFSetBorder(FaceType * &f1, int z1)
{
  assert(f1->FFp(z1)==0 || IsBorder(*f1,z1));

  f1->FFp(z1)=f1;
  f1->FFi(z1)=z1;
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

/**
 * Check if the given face is oriented as the one adjacent to the specified edge.
 * @param f Face to check the orientation
 * @param z Index of the edge
 */
template <class FaceType>
bool CheckOrientation(FaceType &f, int z)
{
	if (IsBorder(f, z))
		return true;
	else
	{
		FaceType *g = f.FFp(z);
		int gi = f.FFi(z);
		if (f.V0(z) == g->V1(gi))
			return true;
		else
			return false;
	}
}


/** 
 * This function change the orientation of the face by inverting the index of two vertex.
 * @param z Index of the edge
 */
template <class FaceType>
void SwapEdge(FaceType &f, const int z) { SwapEdge<FaceType,true>(f,z); }

template <class FaceType, bool UpdateTopology>
void SwapEdge(FaceType &f, const int z)
{
	// swap V0(z) with V1(z)
	std::swap(f.V0(z), f.V1(z));

	if(f.HasFFAdjacency() && UpdateTopology)
	{
		// store information to preserve topology
		int z1 = (z+1)%3;
		int z2 = (z+2)%3;
		FaceType *g1p = f.FFp(z1);
		FaceType *g2p = f.FFp(z2);
		int g1i = f.FFi(z1);
		int g2i = f.FFi(z2);

		// g0 face topology is not affected by the swap

		if (g1p != &f)
		{
			g1p->FFi(g1i) = z2;
			f.FFi(z2) = g1i;
		}
		else
		{
			f.FFi(z2) = z2;
		}

		if (g2p != &f)
		{
			g2p->FFi(g2i) = z1;
			f.FFi(z1) = g2i;
		}
		else
		{
			f.FFi(z1) = z1;
		}

		// finalize swap
		f.FFp(z1) = g2p;
		f.FFp(z2) = g1p;
	}
}

/*!
* Perform a Geometric Check about the normals of a edge flip.
* return trues if after the flip the normals does not change more than the given threshold angle;
* it assumes that the flip is topologically correct.
*
*	\param f	the face
*	\param z	the edge index
*   \param angleRad the threshold angle
*
*  oldD1 ___________ newD1
*       |\          |
*       |  \        |
*       |    \      |
*       |  f  z\    |
*       |        \  |
*       |__________\|
* newD0               oldD0
*/

template <class FaceType>
static bool CheckFlipEdgeNormal(FaceType &f, const int z, const float angleRad)
{
  typedef typename FaceType::VertexType VertexType;
  typedef typename VertexType::CoordType CoordType;
  typedef typename VertexType::ScalarType ScalarType;

  VertexType *OldDiag0 = f.V0(z);
  VertexType *OldDiag1 = f.V1(z);

  VertexType *NewDiag0 = f.V2(z);
  VertexType *NewDiag1 = f.FFp(z)->V2(f.FFi(z));

  assert((NewDiag1 != NewDiag0) && (NewDiag1 != OldDiag0) && (NewDiag1 != OldDiag1));

  CoordType oldN0 = NormalizedNormal( NewDiag0->cP(),OldDiag0->cP(),OldDiag1->cP());
  CoordType oldN1 = NormalizedNormal( NewDiag1->cP(),OldDiag1->cP(),OldDiag0->cP());
  CoordType newN0 = NormalizedNormal( OldDiag0->cP(),NewDiag1->cP(),NewDiag0->cP());
  CoordType newN1 = NormalizedNormal( OldDiag1->cP(),NewDiag0->cP(),NewDiag1->cP());
  if(AngleN(oldN0,newN0) > angleRad) return false;
  if(AngleN(oldN0,newN1) > angleRad) return false;
  if(AngleN(oldN1,newN0) > angleRad) return false;
  if(AngleN(oldN1,newN1) > angleRad) return false;

  return true;
}

/*!
* Perform a Topological check to see if the z-th edge of the face f can be flipped.
* No Geometric test are done. (see CheckFlipEdgeNormal)
*	\param f	pointer to the face
*	\param z	the edge index
*/
template <class FaceType>
static bool CheckFlipEdge(FaceType &f, int z)
{
  typedef typename FaceType::VertexType VertexType;
  typedef typename vcg::face::Pos< FaceType > PosType;

  if (z<0 || z>2)  return false;

	// boundary edges cannot be flipped
  if (face::IsBorder(f, z)) return false;

	FaceType *g = f.FFp(z);
	int		 w = f.FFi(z);

	// check if the vertices of the edge are the same
  // e.g. the mesh has to be well oriented
	if (g->V(w)!=f.V1(z) || g->V1(w)!=f.V(z) )
		return false;

	// check if the flipped edge is already present in the mesh
  // f_v2 and g_v2 are the vertices of the new edge
  VertexType *f_v2 = f.V2(z);
	VertexType *g_v2 = g->V2(w);

  // just a sanity check. If this happens the mesh is not manifold.
  if (f_v2 == g_v2) return false;

  // Now walk around f_v2, one of the two vertexes of the new edge
  // and check that it does not already exists.

  PosType pos(&f, (z+2)%3, f_v2);
  PosType startPos=pos;
	do
	{
		pos.NextE();
    if (g_v2 == pos.VFlip())
			return false;
	}
  while (pos != startPos);

	return true;
}

/*!
* Flip the z-th edge of the face f.
* Check for topological correctness first using <CODE>CheckFlipFace()</CODE>.
*	\param f	pointer to the face
*	\param z	the edge index
*
* Note: For <em>edge flip</em> we intend the swap of the diagonal of the rectangle 
*       formed by the face \a f and the face adjacent to the specified edge.
*/
template <class FaceType>
static void FlipEdge(FaceType &f, const int z)
{	
	assert(z>=0);
	assert(z<3);
	assert( !IsBorder(f,z) );
	assert( face::IsManifold<FaceType>(f, z));

 	FaceType *g = f.FFp(z);
	int		 w = f.FFi(z);
	
	assert( g->V(w)	== f.V1(z) );
	assert( g->V1(w)== f.V(z) );
	assert( g->V2(w)!= f.V(z) );
	assert( g->V2(w)!= f.V1(z) );
	assert( g->V2(w)!= f.V2(z) );

	f.V1(z) = g->V2(w);
	g->V1(w) = f.V2(z);
	
    f.FFp(z)				= g->FFp((w+1)%3);
	f.FFi(z)				= g->FFi((w+1)%3);
    g->FFp(w)				= f.FFp((z+1)%3);
	g->FFi(w)				= f.FFi((z+1)%3);
    f.FFp((z+1)%3)				= g;
	f.FFi((z+1)%3)	= (w+1)%3;
    g->FFp((w+1)%3)			= &f;
	g->FFi((w+1)%3) = (z+1)%3;

	if(f.FFp(z)==g)
	{
		f.FFp(z) = &f;
		f.FFi(z) = z;
	}
	else
	{
		f.FFp(z)->FFp( f.FFi(z) ) = &f;
		f.FFp(z)->FFi( f.FFi(z) ) = z;
	}
	if(g->FFp(w)==&f)
	{
		g->FFp(w)=g;
		g->FFi(w)=w;
	}
	else
	{
		g->FFp(w)->FFp( g->FFi(w) ) = g;
		g->FFp(w)->FFi( g->FFi(w) ) = w;
	}
}


// Stacca la faccia corrente dalla catena di facce incidenti sul vertice z, 
// NOTA funziona SOLO per la topologia VF!!!
// usata nelle classi di collapse
template <class FaceType>
void VFDetach(FaceType & f, int z)
{
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
			++x;
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

/// Append a face in VF list of vertex f->V(z) 
template <class FaceType>
void VFAppend(FaceType* & f, int z)
{
	typename FaceType::VertexType *v = f->V(z);
	if (v->VFp()!=0)
	{
		FaceType *f0=v->VFp();	
		int z0=v->VFi();
		//append
		f->VFp(z)=f0;
		f->VFi(z)=z0;
	}
	v->VFp()=f;
	v->VFi()=z;
}

/*!
* Compute the set of vertices adjacent to a given vertex using VF adjacency. 
*	\param vp	pointer to the vertex whose star has to be computed.
*	\param starVec a std::vector of Vertex pointer that is filled with the adjacent vertices.
*
*/

template <class FaceType>
void VVStarVF( typename FaceType::VertexType* vp, std::vector<typename FaceType::VertexType *> &starVec)
{
	typedef typename FaceType::VertexType* VertexPointer;
	starVec.clear();
	face::VFIterator<FaceType> vfi(vp);
	while(!vfi.End())
			{
				starVec.push_back(vfi.F()->V1(vfi.I()));
				starVec.push_back(vfi.F()->V2(vfi.I()));
				++vfi;
			}
				
	std::sort(starVec.begin(),starVec.end());
	typename std::vector<VertexPointer>::iterator new_end = std::unique(starVec.begin(),starVec.end());
	starVec.resize(new_end-starVec.begin());
}

/*!
* Compute the set of faces adjacent to a given vertex using VF adjacency. 
*	\param vp	pointer to the vertex whose star has to be computed.
*	\param faceVec a std::vector of Face pointer that is filled with the adjacent faces.
*
*/
template <class FaceType>
void VFStarVF( typename FaceType::VertexType* vp, std::vector<FaceType *> &faceVec)
{
	typedef typename FaceType::VertexType* VertexPointer;
	faceVec.clear();
	face::VFIterator<FaceType> vfi(vp);
	while(!vfi.End())
	{
		faceVec.push_back(vfi.F());
		++vfi;
	}
}

/*!
* Compute the ordered set of faces adjacent to a given vertex using VF adjacency.and FF adiacency 
*	\param vp	pointer to the vertex whose star has to be computed.
*	\param faceVec a std::vector of Face pointer that is filled with the adjacent faces.
*
*/
template <class FaceType>
static void VFOrderedStarVF_FF(typename FaceType::VertexType &vp,
								std::vector<FaceType*> &faceVec)
{

	///check that is not on border..
	assert (!vp.IsB());

	///get first face sharing the edge
	FaceType *f_init=vp.VFp();
	int edge_init=vp.VFi(); 

	///and initialize the pos
	vcg::face::Pos<FaceType> VFI(f_init,edge_init);
	bool complete_turn=false;
	do  
	{
		FaceType *curr_f=VFI.F();
		faceVec.push_back(curr_f);

		int curr_edge=VFI.E();

		///assert that is not a border edge
		assert(curr_f->FFp(curr_edge)!=curr_f);

		///continue moving 
		VFI.FlipF();
		VFI.FlipE();

		FaceType *next_f=VFI.F();

		///test if I've finiseh with the face exploration
		complete_turn=(next_f==f_init);
		/// or if I've just crossed a mismatch
	}while (!complete_turn);
}


/*!
* Check if two faces share and edge through the FF topology.
*	\param f0,f1 the two face to be checked
* \param i0,i1 the index of the shared edge;
*/

template <class FaceType>
bool ShareEdgeFF(FaceType *f0,FaceType *f1, int *i0=0, int *i1=0)
{
  assert((!f0->IsD())&&(!f1->IsD()));
  for (int i=0;i<3;i++)
      if (f0->FFp(i)==f1)
      {
        if((i0!=0) && (i1!=0)) {
          *i0=i;
          *i1=f0->FFi(i);
        }
        return true;
      }
  return false;
}

/*!
* Count the number of vertices shared between two faces.
*	\param f0,f1 the two face to be checked
* ;
*/
template <class FaceType>
int CountSharedVertex(FaceType *f0,FaceType *f1)
{
  int sharedCnt=0;
  for (int i=0;i<3;i++)
      for (int j=0;j<3;j++)
          if (f0->V(i)==f1->V(j)) {
                  sharedCnt++;
              }
  return sharedCnt;
}

/*!
* find the first shared vertex between two faces.
*	\param f0,f1 the two face to be checked
* \param i,j the indexes of the shared vertex in the two faces. Meaningful only if there is one single shared vertex
* ;
*/
template <class FaceType>
bool FindSharedVertex(FaceType *f0,FaceType *f1, int &i, int &j)
{
  for (i=0;i<3;i++)
      for (j=0;j<3;j++)
          if (f0->V(i)==f1->V(j)) return true;

  i=-1;j=-1;
  return false;
}

/*!
* find the first shared edge between two faces.
*	\param f0,f1 the two face to be checked
* \param i,j the indexes of the shared edge in the two faces. Meaningful only if there is a shared edge
*
*/
template <class FaceType>
bool FindSharedEdge(FaceType *f0,FaceType *f1, int &i, int &j)
{
  for (i=0;i<3;i++)
      for (j=0;j<3;j++)
        if( ( f0->V0(i)==f1->V0(j) || f0->V0(i)==f1->V1(j) ) &&
            ( f0->V1(i)==f1->V0(j) || f0->V1(i)==f1->V1(j) ) )
            return true;
  i=-1;j=-1;
  return false;
}

/*@}*/
}	 // end namespace
}	 // end namespace

#endif

