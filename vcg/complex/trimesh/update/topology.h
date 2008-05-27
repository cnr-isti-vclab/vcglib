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
Revision 1.20  2008/04/04 10:27:34  cignoni
minor changes to the topology correctness checks

Revision 1.19  2007/05/29 00:07:06  ponchio
VFi++ -> ++VFi

Revision 1.18  2006/02/27 19:26:14  spinelli
minor bug in Face-Face topology loop fixed

Revision 1.17  2006/02/27 11:56:48  spinelli
minor bug in Face-Face topology loop fixed

Revision 1.16  2005/11/10 15:36:42  cignoni
Added clarifying comment in an assert

Revision 1.15  2004/10/20 07:33:10  cignoni
removed FaceBorderFlags (already present in update/flags.h)

Revision 1.14  2004/10/18 17:10:22  ganovelli
added  ::FaceBorderFLags

Revision 1.13  2004/10/01 15:58:00  ponchio
Added include <vector>

Revision 1.12  2004/09/09 13:02:12  ponchio
Linux compatible path in #include

Revision 1.11  2004/08/07 16:18:20  pietroni
addet testFFTopology and testVFTopology functions used to test the rispective topology....

Revision 1.10  2004/07/15 11:35:08  ganovelli
Vfb to VFp

Revision 1.9  2004/07/15 00:13:39  cignoni
Better doxigen documentation

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
#include <vector>
#include <vcg/simplex/face/pos.h>
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
	
	v[0] = pf->V0(nz);
	v[1] = pf->V1(nz);
	assert(v[0] != v[1]); // The face pointed by 'f' is Degenerate (two coincident vertexes)

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

inline bool operator == ( const PEdge & pe ) const
{
	return v[0]==pe.v[0] && v[1]==pe.v[1];
}

};

static void FillEdgeVector(MeshType &m, std::vector<PEdge> &e)
{
		FaceIterator pf;
	typename std::vector<PEdge>::iterator p;
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
}

/** Update the Face-Face topological relation by allowing to retrieve for each face what other faces shares their edges.
*/
static void FaceFace(MeshType &m)
{
  if(!m.HasFFTopology()) return;		
	if( m.fn == 0 ) return;

  std::vector<PEdge> e;
  FillEdgeVector(m,e);
	sort(e.begin(), e.end());							// Lo ordino per vertici

	int ne = 0;											// Numero di edge reali

	typename std::vector<PEdge>::iterator pe,ps;
	ps = e.begin();pe=e.begin();
	//for(ps = e.begin(),pe=e.begin();pe<=e.end();++pe)	// Scansione vettore ausiliario
	do
	{
		if( pe==e.end() || !(*pe == *ps) )					// Trovo blocco di edge uguali
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
		if(pe==e.end()) break;
		++pe;
	} while(true);
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
		(*vi).VFp() = 0;
		(*vi).VFi() = 0;
	}

	for(fi=m.face.begin();fi!=m.face.end();++fi)
	if( ! (*fi).IsD() )
	{
		for(int j=0;j<3;++j)
		{
			(*fi).VFp(j) = (*fi).V(j)->VFp();
			(*fi).VFi(j) = (*fi).V(j)->VFi();
			(*fi).V(j)->VFp() = &(*fi);
			(*fi).V(j)->VFi() = j;
		}
	}
}



/// Auxiliairy data structure for computing face face adjacency information. 
// It identifies and edge storing two vertex pointer and a face pointer where it belong.
class PEdgeTex
{
public:
	
	typename FaceType::TexCoordType  v[2];		// the two Vertex pointer are ordered!
	FacePointer    f;		  // the face where this edge belong
	int      z;				      // index in [0..2] of the edge of the face

  PEdgeTex() {}

void Set( FacePointer  pf, const int nz )
{
	assert(pf!=0);
	assert(nz>=0);
	assert(nz<3);
	
	v[0] = pf->WT(nz);
	v[1] = pf->WT((nz+1)%3);
	assert(v[0] != v[1]); // The face pointed by 'f' is Degenerate (two coincident vertexes)

	if( v[1] < v[0] ) swap(v[0],v[1]);
	f    = pf;
	z    = nz;
}

inline bool operator <  ( const PEdgeTex & pe ) const
{
	if( v[0]<pe.v[0] ) return true;
	else if( pe.v[0]<v[0] ) return false;
	else return v[1] < pe.v[1];
}
inline bool operator == ( const PEdgeTex & pe ) const
{
	return (v[0]==pe.v[0]) && (v[1]==pe.v[1]);
}
inline bool operator != ( const PEdgeTex & pe ) const
{
	return (v[0]!=pe.v[0]) || (v[1]!=pe.v[1]);
}

};


/** Update the Face-Face topological relation by allowing to retrieve for each face what other faces shares their edges.
*/
static void FaceFaceFromTexCoord(MeshType &m)
{
//  assert(HasFFTopology(m));
	assert(HasPerWedgeTexCoord(m));
	
  std::vector<PEdgeTex> e;
	FaceIterator pf;
	typename std::vector<PEdgeTex>::iterator p;

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

	typename std::vector<PEdgeTex>::iterator pe,ps;
	ps = e.begin();pe=e.begin();
	//for(ps = e.begin(),pe=e.begin();pe<=e.end();++pe)	// Scansione vettore ausiliario
	do
	{
		if( pe==e.end() || (*pe) != (*ps) )					// Trovo blocco di edge uguali
		{
			typename std::vector<PEdgeTex>::iterator q,q_next;
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
		if(pe==e.end()) break;
		++pe;
	} while(true);
}





///test correctness of VFtopology
static void TestVertexFace(MeshType &m)
{
	if(!m.HasVFTopology()) return;		

	VertexIterator vi;
	vcg::face::VFIterator<FaceType> VFi;

	for(vi=m.vert.begin();vi!=m.vert.end();++vi)
	{
		if (!vi->IsD())
		if(vi->VFp()!=0) // unreferenced vertices MUST have VF == 0;
		{
			assert(vi->VFp() >= &*m.face.begin());
			assert(vi->VFp() <= &m.face.back());
			VFi.f=vi->VFp();
			VFi.z=vi->VFi();
			while (!VFi.End())
			{
				assert(!VFi.F()->IsD());
				assert((VFi.F()->V(VFi.I()))==&(*vi));
				++VFi;
			}
		}
	}
}

///test correctness of FFtopology
static void TestFaceFace(MeshType &m)
{
	if(!m.HasFFTopology()) return;		

	FaceIterator Fi;
	
	for(Fi=m.face.begin();Fi!=m.face.end();++Fi)
	{
		if (!Fi->IsD())
		{
			for (int i=0;i<3;i++)
			{
				FaceType *f=Fi->FFp(i);
				int e=Fi->FFi(i);
				//invariant property of fftopology
				assert(f->FFp(e)=&(*Fi));
				// Test that the two faces shares the same edge
				VertexPointer v0= Fi->V0(i);
				VertexPointer v1= Fi->V1(i);
				assert( (f->V0(e)==v0) || (f->V1(e)==v0) );
				assert( (f->V0(e)==v1) || (f->V1(e)==v1) );
				
// Old unreadable test
//				assert(((f->V(e) == Fi->V(i))&&(f->V((e+1)%3)==Fi->V((i+1)%3)))||
//					   ((f->V(e)==Fi->V((i+1)%3))&&(f->V((e+1)%3)==Fi->V(i))));
			}
			
		}
	}
}

}; // end class

/*@}*/
}	// End namespace
}	// End namespace


#endif
