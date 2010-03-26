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
#include <vcg/complex/trimesh/base.h>
namespace vcg {
namespace tri {
/// \ingroup trimesh 

/// \headerfile topology.h vcg/complex/trimesh/update/topology.h

/// \brief Generation of per-vertex and per-face topological information.

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


/// \headerfile topology.h vcg/complex/trimesh/update/topology.h

/// \brief Auxiliairy data structure for computing face face adjacency information. 
/** 
It identifies and edge storing two vertex pointer and a face pointer where it belong.
*/

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
	assert(nz<pf->VN());
	
	v[0] = pf->V(nz);
	v[1] = pf->V(pf->Next(nz));
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

// Fill a vector with all the edges of the mesh.
// each edge is stored in the vector the number of times that it appears in the mesh, with the referring face. 
// optionally it can skip the faux edges (to retrieve only the real edges of a triangulated polygonal mesh)

static void FillEdgeVector(MeshType &m, std::vector<PEdge> &e, bool includeFauxEdge=true)
{
	FaceIterator pf;
	typename std::vector<PEdge>::iterator p;
	
	// Alloco il vettore ausiliario
	//e.resize(m.fn*3);			
	FaceIterator fi;
	int n_edges = 0;
	for(fi = m.face.begin(); fi != m.face.end(); ++fi) if(! (*fi).IsD()) n_edges+=(*fi).VN();
	e.resize(n_edges);
	
	p = e.begin();
	for(pf=m.face.begin();pf!=m.face.end();++pf)		
		if( ! (*pf).IsD() )
			for(int j=0;j<(*pf).VN();++j)
			if(includeFauxEdge || !(*pf).IsF(j))
				{
					(*p).Set(&(*pf),j);
					++p;
				}
			
	if(includeFauxEdge) assert(p==e.end());
	else e.resize(p-e.begin());
}

static void FillUniqueEdgeVector(MeshType &m, std::vector<PEdge> &Edges, bool includeFauxEdge=true)
{
	FillEdgeVector(m,Edges,includeFauxEdge);
	sort(Edges.begin(), Edges.end());		// Lo ordino per vertici

	typename std::vector< PEdge>::iterator newEnd = std::unique(Edges.begin(), Edges.end());
	typename std::vector<PEdge>::iterator   ei;

	Edges.resize(newEnd-Edges.begin());
}

/// \brief Update the Face-Face topological relation by allowing to retrieve for each face what other faces shares their edges.
static void FaceFace(MeshType &m)
{
	assert(HasFFAdjacency(m));
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
				//assert((*q).z< 3);
				q_next = q;
				++q_next;
				assert((*q_next).z>=0);
				assert((*q_next).z< (*q_next).f->VN());
				(*q).f->FFp(q->z) = (*q_next).f;				// Collegamento in lista delle facce
				(*q).f->FFi(q->z) = (*q_next).z;
			}
			assert((*q).z>=0);
			assert((*q).z< (*q).f->VN());
			(*q).f->FFp((*q).z) = ps->f;
			(*q).f->FFi((*q).z) = ps->z;
			ps = pe;
			++ne;										// Aggiorno il numero di edge
		}
		if(pe==e.end()) break;
		++pe;
	} while(true);
}


/// \brief Update the Vertex-Face topological relation.
/** 
The function allows to retrieve for each vertex the list of faces sharing this vertex.
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
		for(int j=0;j<(*fi).VN();++j)
		{
			(*fi).VFp(j) = (*fi).V(j)->VFp();
			(*fi).VFi(j) = (*fi).V(j)->VFi();
			(*fi).V(j)->VFp() = &(*fi);
			(*fi).V(j)->VFi() = j;
		}
	}
}


/// \headerfile topology.h vcg/complex/trimesh/update/topology.h

/// \brief Auxiliairy data structure for computing face face adjacency information. 
/** 
It identifies and edge storing two vertex pointer and a face pointer where it belong.
*/

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
	v[1] = pf->WT(pf->Next(nz));
	assert(v[0] != v[1]); // The face pointed by 'f' is Degenerate (two coincident vertexes)

    if( v[1] < v[0] ) std::swap(v[0],v[1]);
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


/// \brief Update the Face-Face topological relation

/** 
The function allows to retrieve for each face what other faces shares their edges.
*/

static void FaceFaceFromTexCoord(MeshType &m)
{
//  assert(HasFFTopology(m));
	assert(HasPerWedgeTexCoord(m));
	
  std::vector<PEdgeTex> e;
	FaceIterator pf;
	typename std::vector<PEdgeTex>::iterator p;

	if( m.fn == 0 ) return;

//	e.resize(m.fn*3);								// Alloco il vettore ausiliario
	FaceIterator fi;
	int n_edges = 0;
	for(fi = m.face.begin(); fi != m.face.end(); ++fi) if(! (*fi).IsD()) n_edges+=(*fi).VN();
	e.resize(n_edges);

	p = e.begin();
	for(pf=m.face.begin();pf!=m.face.end();++pf)			// Lo riempio con i dati delle facce
		if( ! (*pf).IsD() )
			for(int j=0;j<(*pf).VN();++j)
			{
				if( (*pf).WT(j) != (*pf).WT((*pf).Next(j)))
					 {
						(*p).Set(&(*pf),j);
						++p;
					 }
			}
	
	e.resize(p-e.begin());   // remove from the end of the edge vector the unitiailized ones
	assert(p==e.end());
	sort(e.begin(), e.end());		

	int ne = 0;											// number of real edges
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
				assert((*q_next).z< (*q_next).f->VN());
				(*q).f->FFp(q->z) = (*q_next).f;				// Collegamento in lista delle facce
				(*q).f->FFi(q->z) = (*q_next).z;
			}
			assert((*q).z>=0);
			assert((*q).z< (*q).f->VN());
			(*q).f->FFp((*q).z) = ps->f;
			(*q).f->FFi((*q).z) = ps->z;
			ps = pe;
			++ne;										// Aggiorno il numero di edge
		}
		if(pe==e.end()) break;
		++pe;
	} while(true);
}





/// \brief Test correctness of VFtopology
static void TestVertexFace(MeshType &m)
{
	SimpleTempData<typename MeshType::VertContainer, int > numVertex(m.vert,0);

	if(!m.HasVFTopology()) return;		
	
	FaceIterator fi;
	for(fi=m.face.begin();fi!=m.face.end();++fi)
	{
		if (!(*fi).IsD())
		{
			numVertex[(*fi).V0(0)]++;
			numVertex[(*fi).V1(0)]++;
			numVertex[(*fi).V2(0)]++;
		}
	}

	VertexIterator vi;
	vcg::face::VFIterator<FaceType> VFi;

	for(vi=m.vert.begin();vi!=m.vert.end();++vi)
	{
		if (!vi->IsD())
		if(vi->VFp()!=0) // unreferenced vertices MUST have VF == 0;
		{
			int num=0;
			assert(vi->VFp() >= &*m.face.begin());
			assert(vi->VFp() <= &m.face.back());
			VFi.f=vi->VFp();
			VFi.z=vi->VFi();
			while (!VFi.End())
			{
				num++;
				assert(!VFi.F()->IsD());
				assert((VFi.F()->V(VFi.I()))==&(*vi));
				++VFi;
			}
			int num1=numVertex[&(*vi)];
			assert(num==num1);
			/*assert(num>1);*/
		}
	}
}

/// \brief Test correctness of FFtopology (only for 2Manifold Meshes!)
static void TestFaceFace(MeshType &m)
{
	if(!m.HasFFTopology()) return;		

  for(FaceIterator fi=m.face.begin();fi!=m.face.end();++fi)
	{
    if (!fi->IsD())
		{
      for (int i=0;i<(*fi).VN();i++)
			{
        FaceType *ffpi=fi->FFp(i);
        int e=fi->FFi(i);
        //invariant property of FF topology for two manifold meshes
        assert(ffpi->FFp(e) == &(*fi));
        assert(ffpi->FFi(e) == i);

        // Test that the two faces shares the same edge
        // Vertices of the i-th edges of the first face
        VertexPointer v0i= fi->V0(i);
        VertexPointer v1i= fi->V1(i);
        // Vertices of the corresponding edge on the other face
        VertexPointer ffv0i= ffpi->V0(e);
        VertexPointer ffv1i= ffpi->V1(e);

        assert( (ffv0i==v0i) || (ffv0i==v1i) );
        assert( (ffv1i==v0i) || (ffv1i==v1i) );
			}
			
		}
	}
}

}; // end class

}	// End namespace
}	// End namespace


#endif
