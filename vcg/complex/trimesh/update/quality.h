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
Revision 1.2  2004/05/10 13:43:00  cignoni
Added use of VFIterator in VertexGeodesicFromBorder

Revision 1.1  2004/03/31 14:59:14  cignoni
First working version!

Revision 1.2  2004/03/29 14:26:57  cignoni
First working version!

****************************************************************************/

#ifndef __VCG_TRI_UPDATE_QUALITY
#define __VCG_TRI_UPDATE_QUALITY
#include <vcg/simplex/face/pos.h>

namespace vcg {
namespace tri {
/** \addtogroup trimesh */
/*@{*/
/// Generation of per-vertex and per-face Qualities according to various strategy, like geodesic distance from the border (UpdateQuality::VertexGeodesicFromBorder) or curvature ecc.

template <class UpdateMeshType>
class UpdateQuality
{
public:
  typedef UpdateMeshType MeshType; 
  typedef typename MeshType::VertexType     VertexType;
  typedef typename MeshType::VertexPointer  VertexPointer;
  typedef typename MeshType::VertexIterator VertexIterator;
  typedef typename MeshType::FaceType       FaceType;
  typedef typename MeshType::FacePointer    FacePointer;
  typedef typename MeshType::FaceIterator   FaceIterator;

class VQualityHeap
{
public:
	float q;
	VertexPointer p;
	inline VQualityHeap( VertexPointer np )
	{
		q = np->Q();
		p = np;
	}
		// Attenzione il minore e' maggiore
	inline bool operator <  ( const VQualityHeap & vq ) const { return q >  vq.q; }
	inline bool operator == ( const VQualityHeap & vq ) const { return q == vq.q; }
	inline bool operator >  ( const VQualityHeap & vq ) const { return q <  vq.q; }
	inline bool operator != ( const VQualityHeap & vq ) const { return q != vq.q; }
	inline bool operator <= ( const VQualityHeap & vq ) const { return q >= vq.q; }
	inline bool operator >= ( const VQualityHeap & vq ) const { return q <= vq.q; }
	inline bool is_valid() const { return q==p->Q(); }
};



// REQUIREMENT VF topology and Border FLags
// Calcola la qualita' come distanza geodesica dal bordo della mesh.
// Robusta funziona anche per mesh non manifold.
// La qualita' memorizzata indica la distanza assoluta dal bordo della mesh.
// Nota prima del 13/11/03 in alcuni casi rari SPT andava in loop perche' poteva capitare
// che per approx numeriche ben strane pw->Q() > pv->Q()+d ma durante la memorizzazione 
// della nuova distanza essa rimanesse uguale a prima. Patchato rimettendo i vertici nello 
// heap solo se migliorano la distanza di un epsilon == 1/100000 della mesh diag.
/** Compute, for each vertex of the mesh the geodesic distance from the border of the mesh itself;
Requirements: VF topology, Per Vertex Quality;
it uses the classical dijkstra Shortest Path Tree algorithm. 
The geodesic distance is approximated by allowing to walk only along edges of the mesh.
*/
static void VertexGeodesicFromBorder(MeshType &m)	// R1
{
	//Requirements
	assert(m.HasVFTopology());
	assert(m.HasPerVertexQuality());

	vector<VQualityHeap > heap;
	VertexIterator v;
	FaceIterator f;
	int j;

	for(v=m.vert.begin();v!=m.vert.end();++v)
		(*v).Q() = -1;
	for(f=m.face.begin();f!=m.face.end();++f)			// Inserisco nell'heap i v di bordo
		if(!(*f).IsD())
			for(j=0;j<3;++j)
				if( (*f).IsB(j) )
				{
					for(int k=0;k<2;++k)
					{
						VertexPointer pv = (*f).V((j+k)%3);
						if( pv->Q()==-1 )
						{
							pv->Q() = 0;
							heap.push_back(VQualityHeap(pv));
						}
					}
				}
	
 const MeshType::ScalarType loc_eps=m.bbox.Diag()/MeshType::ScalarType(100000);
 while( heap.size()!=0 )							// Shortest path tree
	{
		VertexPointer pv;
		pop_heap(heap.begin(),heap.end());
		if( ! heap.back().is_valid() )
		{
			heap.pop_back();
			continue;
		}
		pv = heap.back().p;
		heap.pop_back();
	  
		for(face::VFIterator<FaceType> vfi(pv) ; !vfi.End(); ++vfi )
		{
			for(int k=0;k<2;++k)
			{
				VertexPointer pw;
				float d;
				if(k==0) pw = vfi.f->V1(vfi.z);
				else     pw = vfi.f->V2(vfi.z);
				d = Distance(pv->P(),pw->P());
				if( pw->Q()==-1 || pw->Q() > pv->Q()+d + loc_eps)
				{
					pw->Q() = pv->Q()+d;
					heap.push_back(VQualityHeap(pw));
					push_heap(heap.begin(),heap.end());
				}
			}
		}
	}

	for(v=m.vert.begin();v!=m.vert.end();++v)
		if(v->Q()==-1)
			v->Q() = 0;
}

}; //end class
} // end namespace
} // end namespace
#endif