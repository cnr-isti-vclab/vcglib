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
Revision 1.8  2004/06/02 16:42:44  ganovelli
typename for gcc compilation

Revision 1.7  2004/06/02 16:28:22  ganovelli
minor changes (swap =>> math::Swap)

Revision 1.6  2004/05/10 15:23:43  cignoni
Changed a FV -> VF in VertexFace topology computation

Revision 1.5  2004/05/06 15:24:38  pietroni
changed names to topology functions

Revision 1.4  2004/03/31 14:44:43  cignoni
Added Vertex-Face Topology

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
/** Generation of per-vertex and per-face topological information.
**/

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

	if( v[0] > v[1] ) math::Swap(v[0],v[1]);
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


/** Update the Face-Face topological relation by allowing to retrieve for each face what other faces shares their edges.
*/
static void FaceFace(MeshType &m)
{
  if(!m.HasFFTopology()) return;		

  std::vector<PEdge> e;
	FaceIterator pf;
	typename std::vector<PEdge>::iterator p;

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

	typename std::vector<PEdge>::iterator pe,ps;
	for(ps = e.begin(),pe=e.begin();pe<=e.end();++pe)	// Scansione vettore ausiliario
	{
		if( pe==e.end() || *pe != *ps )					// Trovo blocco di edge uguali
		{
			typename std::vector<PEdge>::iterator q,q_next;
			for (q=ps;q<pe-1;++q)						// Scansione facce associate
			{
				assert((*q).z>=0);
				assert((*q).z< 3);
				q_next = q;
				++q_next;
				assert((*q_next).z>=0);
				assert((*q_next).z< 3);
				(*q).f->FFp(q->z) = (*q_next).f;				// Collegamento in lista delle facce
				(*q).f->FFi(q->z) = (*q_next).z;
			}
			assert((*q).z>=0);
			assert((*q).z< 3);
			(*q).f->FFp((*q).z) = ps->f;
			(*q).f->FFi((*q).z) = ps->z;
			ps = pe;
			++ne;										// Aggiorno il numero di edge
		}
	}
}

/** Update the Vertex-Face topological relation by allowing to retrieve for each vertex the list of faces sharing this vertex..
*/
static void VertexFace(MeshType &m)
{
  if(!m.HasVFTopology()) return;		

	VertexIterator vi;
	FaceIterator fi;

	for(vi=m.vert.begin();vi!=m.vert.end();++vi)
	{
		(*vi).VFb() = 0;
		(*vi).VFi() = 0;
	}

	for(fi=m.face.begin();fi!=m.face.end();++fi)
	if( ! (*fi).IsD() )
	{
		for(int j=0;j<3;++j)
		{
			(*fi).VFp(j) = (*fi).V(j)->VFb();
			(*fi).VFi(j) = (*fi).V(j)->VFi();
			(*fi).V(j)->VFb() = &(*fi);
			(*fi).V(j)->VFi() = j;
		}
	}
}
	

}; // end class

/*@}*/
}	// End namespace
}	// End namespace


#endif
