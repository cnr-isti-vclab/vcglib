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
Revision 1.3  2004/03/12 15:22:19  cignoni
Written some documentation and added to the trimes doxygen module

Revision 1.2  2004/03/05 21:49:21  cignoni
First working version for face face

Revision 1.1  2004/03/04 00:53:24  cignoni
Initial commit


****************************************************************************/
#ifndef __VCG_TRI_UPDATE_TOPOLOGY
#define __VCG_TRI_UPDATE_TOPOLOGY
#include <algorithm>
namespace vcg {
namespace tri {
/** \addtogroup trimesh */
/*@{*/

template <class UpdateMeshType>
class UpdateTopology
{

public:
typedef UpdateMeshType MeshType; 
typedef typename MeshType::VertexType     VertexType;
typedef typename MeshType::VertexPointer  VertexPointer;
typedef typename MeshType::VertexIterator VertexIterator;
typedef typename MeshType::FaceType       FaceType;
typedef typename MeshType::FacePointer    FacePointer;
typedef typename MeshType::FaceIterator   FaceIterator;



/// Auxiliairy data structure for computing face face adjacency information. 
// It identifies and edge storing two vertex pointer and a face pointer where it belong.
class PEdge
{
public:
	
	VertexPointer  v[2];		// the two Vertex pointer are ordered!
	FacePointer    f;		  // the face where this edge belong
	int      z;				      // index in [0..2] of the edge of the face

  PEdge() {}

void Set( FacePointer  pf, const int nz )
{
	assert(pf!=0);
	assert(nz>=0);
	assert(nz<3);
	
	v[0] = pf->V(nz);
	v[1] = pf->V((nz+1)%3);
	assert(v[0] != v[1]);

	if( v[0] > v[1] ) swap(v[0],v[1]);
	f    = pf;
	z    = nz;
}

inline bool operator <  ( const PEdge & pe ) const
{
	if( v[0]<pe.v[0] ) return true;
	else if( v[0]>pe.v[0] ) return false;
	else return v[1] < pe.v[1];
}

inline bool operator <=  ( const PEdge & pe ) const
{
	if( v[0]<pe.v[0] ) return true;
	else if( v[0]>pe.v[0] ) return false;
	else return v[1] <= pe.v[1];
}

inline bool operator >  ( const PEdge & pe ) const
{
	if( v[0]>pe.v[0] ) return true;
	else if( v[0]<pe.v[0] ) return false;
	else return v[1] > pe.v[1];
}

inline bool operator >=  ( const PEdge & pe ) const
{
	if( v[0]>pe.v[0] ) return true;
	else if( v[0]<pe.v[0] ) return false;
	else return v[1] >= pe.v[1];
}

inline bool operator == ( const PEdge & pe ) const
{
	return v[0]==pe.v[0] && v[1]==pe.v[1];
}

inline bool operator != ( const PEdge & pe ) const
{
	return v[0]!=pe.v[0] || v[1]!=pe.v[1];
}
};



static void FaceFace(MeshType &m)
{
  if(!m.HasFFTopology()) return;		

  vector<PEdge> e;
	FaceIterator pf;
	vector<PEdge>::iterator p;

	if( m.fn == 0 ) return;

	e.resize(m.fn*3);								// Alloco il vettore ausiliario
	p = e.begin();
	for(pf=m.face.begin();pf!=m.face.end();++pf)			// Lo riempio con i dati delle facce
		if( ! (*pf).IsD() )
			for(int j=0;j<3;++j)
			{
				(*p).Set(&(*pf),j);
				++p;
			}
	assert(p==e.end());
	sort(e.begin(), e.end());							// Lo ordino per vertici

	int ne = 0;											// Numero di edge reali

	vector<PEdge>::iterator pe,ps;
	for(ps = e.begin(),pe=e.begin();pe<=e.end();++pe)	// Scansione vettore ausiliario
	{
		if( pe==e.end() || *pe != *ps )					// Trovo blocco di edge uguali
		{
			vector<PEdge>::iterator q,q_next;
			for (q=ps;q<pe-1;++q)						// Scansione facce associate
			{
				assert((*q).z>=0);
				assert((*q).z< 3);
				q_next = q;
				++q_next;
				assert((*q_next).z>=0);
				assert((*q_next).z< 3);
				(*q).f->F(q->z) = (*q_next).f;				// Collegamento in lista delle facce
				(*q).f->Z(q->z) = (*q_next).z;
			}
			assert((*q).z>=0);
			assert((*q).z< 3);
			(*q).f->F((*q).z) = ps->f;
			(*q).f->Z((*q).z) = ps->z;
			ps = pe;
			++ne;										// Aggiorno il numero di edge
		}
	}
}

static void VertexFace(MeshType &m)
{
  if(!m.HasVFTopology()) return;		

	VertexIterator vi;
	FaceIterator fi;

	for(vi=m.vert.begin();vi!=m.vert.end();++vi)
	{
		(*vi).Fp() = 0;
		(*vi).Zp() = 0;
	}

	for(fi=m.face.begin();fi!=m.face.end();++fi)
	if( ! (*fi).IsD() )
	{
		for(int j=0;j<3;++j)
		{
			(*fi).Fv(j) = (*fi).V(j)->Fp();
			(*fi).Zv(j) = (*fi).V(j)->Zp();
			(*fi).V(j)->Fp() = &(*fi);
			(*fi).V(j)->Zp() = j;
		}
	}
}
	

}; // end class

/*@}*/
}	// End namespace
}	// End namespace


#endif
