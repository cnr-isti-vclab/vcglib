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
Revision 1.20  2007/05/31 15:24:50  ponchio
FIxed off-by-one error on FaceBorderFromNone.

Revision 1.19  2007/05/22 15:19:42  cignoni
Added VertexClear

Revision 1.18  2007/01/30 18:49:23  tarini
aggiunta la VertexBorderFromNone (flag bordo per vertici senza richiedere nulla)

Revision 1.17  2006/08/31 13:11:12  marfr960
corrected bounds of a vector scan

Revision 1.16  2006/08/30 12:59:49  marfr960
Added missing std:: to swap

Revision 1.15  2006/08/30 06:50:07  cignoni
Reverted to version 1.13. Version 1.14 was done on outdated version.

Revision 1.13  2006/06/18 20:49:30  cignoni
Added missing IsD tests

Revision 1.12  2006/05/03 21:23:25  cignoni
Corrected IsDeleted -> isD

Revision 1.11  2005/12/02 00:09:12  cignoni
Added assert(HasFlags) everywhere..

Revision 1.10  2005/07/06 08:16:34  ganovelli
set VertexBorderFromFace as static

Revision 1.9  2005/06/10 15:07:23  cignoni
Completed FaceBorderFromNone (and added a missing helper class)

Revision 1.8  2005/04/01 13:04:55  fiorin
Minor changes

Revision 1.7  2004/09/14 19:49:43  ganovelli
first compilation version

Revision 1.6  2004/07/15 00:13:39  cignoni
Better doxigen documentation

Revision 1.5  2004/07/06 06:27:02  cignoni
Added  FaceBorderFromVF

Revision 1.4  2004/05/13 15:58:55  ganovelli
function Clear added

Revision 1.3  2004/03/12 15:22:19  cignoni
Written some documentation and added to the trimes doxygen module

Revision 1.2  2004/03/10 00:46:10  cignoni
changed to the face::IsBorder() style

Revision 1.1  2004/03/05 10:59:24  cignoni
Changed name from plural to singular (normals->normal)

Revision 1.1  2004/03/04 00:37:56  cignoni
First working version!


****************************************************************************/
#ifndef __VCG_TRI_UPDATE_FLAGS
#define __VCG_TRI_UPDATE_FLAGS

#include <vcg/simplex/face/pos.h>

namespace vcg {
namespace tri {
/// \ingroup trimesh 

/// \headerfile flag.h vcg/complex/trimesh/update/flag.h

/// \brief Management, updating and computation of per-vertex and per-face flags (like border flags).

/** 
This class is used to compute or update some of the flags that can be stored in the mesh components. For now just Border flags (e.g. the flag that tells if a given edge of a face belong to a border of the mesh or not).
*/

template <class UpdateMeshType>
class UpdateFlags
{

public:
typedef UpdateMeshType MeshType; 
typedef vcg::face::Pos<typename UpdateMeshType::FaceType> PosType;
typedef typename MeshType::ScalarType     ScalarType;
typedef typename MeshType::VertexType     VertexType;
typedef typename MeshType::VertexPointer  VertexPointer;
typedef typename MeshType::VertexIterator VertexIterator;
typedef typename MeshType::FaceType       FaceType;
typedef typename MeshType::FacePointer    FacePointer;
typedef typename MeshType::FaceIterator   FaceIterator;

/// \brief Reset all the mesh flags (both vertexes and faces) setting everithing to zero (the default value for flags)

static void Clear(MeshType &m)
{
  assert(HasPerFaceFlags(m));
	FaceIterator fi;
	VertexIterator vi;
	for(fi=m.face.begin(); fi!=m.face.end(); ++fi)
		(*fi).Flags() = 0;
	for(vi=m.vert.begin(); vi!=m.vert.end(); ++vi)
		(*vi).Flags() = 0;
}

static void VertexClear(MeshType &m, unsigned int FlagMask = 0xffffffff)
{
	VertexIterator vi;
	int andMask = ~FlagMask;
	for(vi=m.vert.begin(); vi!=m.vert.end(); ++vi)
		if(!(*vi).IsD()) (*vi).Flags() &= andMask ;
}

static void FaceClear(MeshType &m, unsigned int FlagMask = 0xffffffff)
{
	FaceIterator fi;
	int andMask = ~FlagMask;
	for(fi=m.face.begin(); fi!=m.face.end(); ++fi)
		if(!(*fi).IsD()) (*fi).Flags() &= andMask ;
}

static void VertexSet(MeshType &m, unsigned int FlagMask)
{
	VertexIterator vi;
	for(vi=m.vert.begin(); vi!=m.vert.end(); ++vi)
		if(!(*vi).IsD()) (*vi).Flags() |= FlagMask ;
}

static void FaceSet(MeshType &m, unsigned int FlagMask)
{
	FaceIterator fi;
  for(fi=m.face.begin(); fi!=m.face.end(); ++fi)
		if(!(*fi).IsD()) (*fi).Flags() |= FlagMask ;
}



static void VertexClearV(MeshType &m) { VertexClear(m,VertexType::VISITED);}
static void VertexClearB(MeshType &m) { VertexClear(m,VertexType::BORDER);}
static void FaceClearV(MeshType &m) { FaceClear(m,FaceType::VISITED);}
static void FaceClearB(MeshType &m) { FaceClear(m,FaceType::BORDER012);}
static void FaceClearF(MeshType &m) { FaceClear(m,FaceType::FAUX012);}

static void VertexSetV(MeshType &m) { VertexSet(m,VertexType::VISITED);}
static void VertexSetB(MeshType &m) { VertexSet(m,VertexType::BORDER);}
static void FaceSetV(MeshType &m) { FaceSet(m,FaceType::VISITED);}
static void FaceSetB(MeshType &m) { FaceSet(m,FaceType::BORDER);}
static void FaceSetF(MeshType &m) { FaceSet(m,FaceType::FAUX012);}

/// \brief Compute the border flags for the faces using the Face-Face Topology. 

/**
 \warning Obviously it assumes that the topology has been correctly computed (see: UpdateTopology::FaceFace )
*/
static void FaceBorderFromFF(MeshType &m)
{
  assert(HasPerFaceFlags(m));
//	const int BORDERFLAG[3]={FaceType::BORDER0,FaceType::BORDER1,FaceType::BORDER2};
	FaceIterator fi;
	for(fi=m.face.begin();fi!=m.face.end();++fi)if(!(*fi).IsD())
		for(int j=0;j<3;++j)
		{
			//if(!(*fi).IsManifold(j)) (*fi).SetCF(j);
			//else 
      if(face::IsBorder(*fi,j)) (*fi).SetB(j);
					 else (*fi).ClearB(j);
		}
}


static void FaceBorderFromVF(MeshType &m)
{
  assert(HasPerFaceFlags(m));
	VertexIterator vi;
  assert(m.HasVFTopology());

  int visitedBit=VertexType::NewBitFlag();

	// Calcolo dei bordi
  // per ogni vertice vi si cercano i vertici adiacenti che sono toccati da una faccia sola
	// (o meglio da un numero dispari di facce)

		const int BORDERFLAG[3]={FaceType::BORDER0, FaceType::BORDER1, FaceType::BORDER2};
		
		for(vi=m.vert.begin();vi!=m.vert.end();++vi)
			if(!(*vi).IsD())
				{
					for(face::VFIterator<FaceType> vfi(&*vi) ; !vfi.End(); ++vfi )
          {
			          vfi.f->V1(vfi.z)->ClearUserBit(visitedBit);
								vfi.f->V2(vfi.z)->ClearUserBit(visitedBit);
          }
					for(face::VFIterator<FaceType> vfi(&*vi) ; !vfi.End(); ++vfi )
          {
            if(vfi.f->V1(vfi.z)->IsUserBit(visitedBit))  vfi.f->V1(vfi.z)->ClearUserBit(visitedBit);
																                    else vfi.f->V1(vfi.z)->SetUserBit(visitedBit);
						if(vfi.f->V2(vfi.z)->IsUserBit(visitedBit))  vfi.f->V2(vfi.z)->ClearUserBit(visitedBit);
																                    else vfi.f->V2(vfi.z)->SetUserBit(visitedBit);
					}
          for(face::VFIterator<FaceType> vfi(&*vi) ; !vfi.End(); ++vfi )
          {
	        		if(vfi.f->V(vfi.z)< vfi.f->V1(vfi.z)  &&  vfi.f->V1(vfi.z)->IsUserBit(visitedBit)) 
					        vfi.f->Flags() |= BORDERFLAG[vfi.z];
	        		if(vfi.f->V(vfi.z)< vfi.f->V2(vfi.z)  &&  vfi.f->V2(vfi.z)->IsUserBit(visitedBit)) 
					        vfi.f->Flags() |= BORDERFLAG[(vfi.z+2)%3];
          }
			}	
		VertexType::DeleteBitFlag(VertexType::LastBitFlag());
}

 
class EdgeSorter
{
public:
	
	VertexPointer v[2];		// Puntatore ai due vertici (Ordinati)
	FacePointer    f;				// Puntatore alla faccia generatrice
	int      z;				// Indice dell'edge nella faccia

  EdgeSorter() {} // Nothing to do


void Set( const FacePointer pf, const int nz )
{
	assert(pf!=0);
	assert(nz>=0);
	assert(nz<3);
	
	v[0] = pf->V(nz);
	v[1] = pf->V((nz+1)%3);
	assert(v[0] != v[1]);

	if( v[0] > v[1] ) std::swap(v[0],v[1]);
	f    = pf;
	z    = nz;
}

inline bool operator <  ( const EdgeSorter & pe ) const {
	if( v[0]<pe.v[0] ) return true;
	else if( v[0]>pe.v[0] ) return false;
	else return v[1] < pe.v[1];
}

inline bool operator == ( const EdgeSorter & pe ) const
{
	return v[0]==pe.v[0] && v[1]==pe.v[1];
}
inline bool operator != ( const EdgeSorter & pe ) const
{
	return v[0]!=pe.v[0] || v[1]!=pe.v[1];
}

};


// versione minimale che non calcola i complex flag.
static void VertexBorderFromNone(MeshType &m)
{
  assert(HasPerVertexFlags(m));
	std::vector<EdgeSorter> e;
	typename UpdateMeshType::FaceIterator pf;
	typename std::vector<EdgeSorter>::iterator p;

	if( m.fn == 0 ) 
		return;

	e.resize(m.fn*3);								// Alloco il vettore ausiliario
	p = e.begin();
	for(pf=m.face.begin();pf!=m.face.end();++pf)			// Lo riempio con i dati delle facce
		if( ! (*pf).IsD() )
			for(int j=0;j<3;++j)
			{
				(*p).Set(&(*pf),j);
				(*pf).ClearB(j);
				++p;
			}
	assert(p==e.end());
	sort(e.begin(), e.end());							// Lo ordino per vertici
	
	typename std::vector<EdgeSorter>::iterator pe,ps;
	for(ps = e.begin(), pe = e.begin(); pe < e.end(); ++pe)	// Scansione vettore ausiliario
	{
		if( pe==e.end() ||  *pe != *ps )					// Trovo blocco di edge uguali
		{
			if(pe-ps==1) 	{	
					ps->v[0]->SetB();
					ps->v[1]->SetB();
			} else
			if(pe-ps!=2)  {  // not twomanyfold!
				for(;ps!=pe;++ps) {
						ps->v[0]->SetB(); // Si settano border anche i complex.
						ps->v[1]->SetB();
        }
			} 
			ps = pe;
		}
	}
}

/// This function fill the flags with the info on what is the best projection direction
/// for a given face. Used by the point-face distance function when do not exploiting pre-computed 
/// per-face data (the so called edge component)  
static void FaceProjection(MeshType &m)
{
	FaceIterator fi;
	for(fi=m.face.begin();fi!=m.face.end();++fi)			// Lo riempio con i dati delle facce
		if( ! (*fi).IsD() )
		{
			ScalarType nx = math::Abs((*fi).cN()[0]);
			ScalarType ny = math::Abs((*fi).cN()[1]);
			ScalarType nz = math::Abs((*fi).cN()[2]);
			if(nx>ny && nx>nz) { (*fi).Flags() |= FaceType::NORMX; }
			else if(ny>nz)     { (*fi).Flags() |= FaceType::NORMY; }
			else               { (*fi).Flags() |= FaceType::NORMZ; }
		}
}

/// Computes per-face border flags without requiring any kind of topology 
/// It has a O(fn log fn) complexity. 
static void FaceBorderFromNone(MeshType &m)
{
  assert(HasPerFaceFlags(m));
	std::vector<EdgeSorter> e;
	typename UpdateMeshType::FaceIterator pf;
	typename std::vector<EdgeSorter>::iterator p;

	for(VertexIterator v=m.vert.begin();v!=m.vert.end();++v)
		(*v).ClearB();
		
	if( m.fn == 0 ) 
		return;

	FaceIterator fi;
	int n_edges = 0;
	for(fi = m.face.begin(); fi != m.face.end(); ++fi) if(! (*fi).IsD()) n_edges+=(*fi).VN();
	e.resize(n_edges);

	p = e.begin();
	for(pf=m.face.begin();pf!=m.face.end();++pf)			// Lo riempio con i dati delle facce
		if( ! (*pf).IsD() )
			for(int j=0;j<(*pf).VN();++j)
			{
				(*p).Set(&(*pf),j);
				(*pf).ClearB(j);
				++p;
			}
	assert(p==e.end());
	sort(e.begin(), e.end());							// Lo ordino per vertici
	 
	typename std::vector<EdgeSorter>::iterator pe,ps;
	ps = e.begin();pe=e.begin();
	do
	{	 
		if( pe==e.end() ||  *pe != *ps )					// Trovo blocco di edge uguali
		{
			if(pe-ps==1) 	{	
					ps->f->SetB(ps->z);
			} else
			if(pe-ps!=2)  {  // Caso complex!!
				for(;ps!=pe;++ps)
						ps->f->SetB(ps->z); // Si settano border anche i complex.
			} 
			ps = pe;
		}
		if(pe==e.end()) break;
		++pe;
	} while(true);
//	TRACE("found %i border (%i complex) on %i edges\n",nborder,ncomplex,ne);
}

/// Compute the PerVertex Border flag deriving it from the faces
static void VertexBorderFromFace(MeshType &m)
{
  assert(HasPerFaceFlags(m));
	typename MeshType::VertexIterator v;
	typename MeshType::FaceIterator f;

	for(v=m.vert.begin();v!=m.vert.end();++v)
		(*v).ClearB();

	for(f=m.face.begin();f!=m.face.end();++f)
    if(!(*f).IsD())
	    {
		    for(int z=0;z<(*f).VN();++z)
			    if( (*f).IsB(z) )
			    {
				    (*f).V(z)->SetB();
				    (*f).V((*f).Next(z))->SetB();
			    }
	    }
}
//
static void FaceFauxCrease(MeshType &m,float AngleRad)
{
  assert(HasPerFaceFlags(m));
  assert(HasFFAdjacency(m));

  typename MeshType::FaceIterator f;

  //initially everything is faux (e.g all internal)
  FaceSetF(m);
  for(f=m.face.begin();f!=m.face.end();++f)
  {
    if(!(*f).IsD())
    {
      for(int z=0;z<(*f).VN();++z)
      {
        if( face::IsBorder(*f,z) )  (*f).ClearF(z);
        else
        {
          if(Angle((*f).N(), (*f).FFp(z)->N()) > AngleRad)
            (*f).ClearF(z);
        }
      }
    }
  }
}

}; // end class

}	// End namespace tri
}	// End namespace vcg


#endif
