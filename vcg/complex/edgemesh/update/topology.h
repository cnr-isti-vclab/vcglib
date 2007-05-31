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
Revision 1.3  2004/10/28 00:47:51  cignoni
Better Doxygen documentation

Revision 1.2  2004/05/10 14:42:17  ganovelli
nimor changes


****************************************************************************/
#ifndef __VCG_EDGE_UPDATE_TOPOLOGY
#define __VCG_EDGE_UPDATE_TOPOLOGY
#include <algorithm>
namespace vcg {
namespace edge {
/** \addtogroup edgemesh */
/*@{*/

template <class UpdateMeshType>
class UpdateTopology
{

public:
typedef UpdateMeshType MeshType; 
typedef typename MeshType::VertexType     VertexType;
typedef typename MeshType::VertexPointer  VertexPointer;
typedef typename MeshType::VertexIterator VertexIterator;
typedef typename MeshType::EdgeType       EdgeType;
typedef typename MeshType::EdgePointer    EdgePointer;
typedef typename MeshType::EdgeIterator   EdgeIterator;



/// Auxiliairy data structure for computing face face adjacency information. 
// It identifies and edge storing two vertex pointer and a face pointer where it belong.
class PVertex
{
public:
	
	VertexPointer  v;		// the two Vertex pointer are ordered!
	EdgePointer    e;		  // the edge where this vertex belong
	int      z;				      // index in [0..2] of the edge of the face

  PVertex() {}

void Set( EdgePointer  pe, const int nz )
{
	assert(pe!=0);
	assert(nz>=0);
	assert(nz<2);
	
	v= pe->V(nz);
	e    = pe;
	z    = nz;
}

inline bool operator <  ( const PVertex & pe ) const
{
	return ( v<pe.v );
}

inline bool operator <=  ( const PVertex & pe ) const
{
	return ( v<=pe.v );
}

inline bool operator >  ( const PVertex & pe ) const
{
	return ( v>pe.v );
}

inline bool operator >=  ( const PVertex & pe ) const
{
	return( v>pe.v );
}

inline bool operator == ( const PVertex & pe ) const
{
	return (v==pe.v);
}

inline bool operator != ( const PVertex & pe ) const
{
	return (v!=pe.v || v!=pe.v);
}
};



static void EdgeEdge(MeshType &m)
{
  if(!m.HasEETopology()) return;		

  vector<PVertex> v;
	EdgeIterator pf;
	typename vector<PVertex>::iterator p;

	if( m.en == 0 ) return;

	v.resize(m.en*2);								// Alloco il vettore ausiliario
	p = v.begin();
	for(pf=m.edges.begin();pf!=m.edges.end();++pf)			// Lo riempio con i dati delle facce
		if( ! (*pf).IsD() )
			for(int j=0;j<2;++j)
			{
				(*p).Set(&(*pf),j);
				++p;
			}
	assert(p==v.end());
	sort(v.begin(), v.end());							// Lo ordino per vertici

	int ne = 0;											// Numero di edge reali

	typename vector<PVertex>::iterator pe,ps;
	for(ps = v.begin(),pe=v.begin();pe<=v.end();++pe)	// Scansione vettore ausiliario
	{
		if( pe==v.end() || *pe != *ps )					// Trovo blocco di edge uguali
		{
			typename vector<PVertex>::iterator q,q_next;
			for (q=ps;q<pe-1;++q)						// Scansione edge associati
			{
				assert((*q).z>=0);
				assert((*q).z< 2);
				q_next = q;
				++q_next;
				assert((*q_next).z>=0);
				assert((*q_next).z< 2);
				(*q).e->EEp(q->z) = (*q_next).e;				// Collegamento in lista delle facce
				(*q).e->EEi(q->z) = (*q_next).z;
			}
			assert((*q).z>=0);
			assert((*q).z< 3);
			(*q).e->EEp((*q).z) = ps->e;
			(*q).e->EEi((*q).z) = ps->z;
			ps = pe;
			++ne;										// Aggiorno il numero di edge
		}
	}
}

static void VertexEdge(MeshType &m)
{
  if(!m.HasVETopology()) return;		

	VertexIterator vi;
	EdgeIterator ei;

	for(vi=m.vert.begin();vi!=m.vert.end();++vi)
	{
		(*vi).Ep() = 0;
		(*vi).Ei() = 0;
	}

	for(ei=m.edges.begin();ei!=m.edges.end();++ei)
	if( ! (*ei).IsD() )
	{
		for(int j=0;j<2;++j)
		{
			(*ei).Ev(j) = (*ei).V(j)->Ep();
			(*ei).Zv(j) = (*ei).V(j)->Ei();
			(*ei).V(j)->Ep() = &(*ei);
			(*ei).V(j)->Ei() = j;
		}
	}
}
	

}; // end class

/*@}*/
}	// End namespace
}	// End namespace


#endif
