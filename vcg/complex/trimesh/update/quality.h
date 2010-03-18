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
Revision 1.7  2005/03/17 16:16:08  cignoni
removed small gcc compiling issues

Revision 1.6  2005/03/15 11:48:50  cignoni
Added missing include assert and improved comments and requirements of geodesic quality

Revision 1.5  2004/07/15 00:13:39  cignoni
Better doxigen documentation

Revision 1.4  2004/07/06 06:29:53  cignoni
removed assumption of a using namespace std and added a missing include

Revision 1.3  2004/06/24 15:15:12  cignoni
Better Doxygen documentation

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
#include <algorithm>
#include <vector>
#include <assert.h>

namespace vcg {
namespace tri {
/// \ingroup trimesh 

/// \headerfile quality.h vcg/complex/trimesh/update/quality.h

/// \brief Generation of per-vertex and per-face qualities. 
/**
 It works according to various strategy, like geodesic distance from the border (UpdateQuality::VertexGeodesicFromBorder) or curvature ecc.
 This class is templated over the mesh and (like all other Update* classes) has only static members; Typical usage:
\code
MyMeshType m;
UpdateQuality<MyMeshType>::VertexGeodesicFromBorder(m);
\endcode
*/

template <class UpdateMeshType>
class UpdateQuality
{
public:
  typedef UpdateMeshType MeshType; 
  typedef typename MeshType::ScalarType     ScalarType;
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



// *** IMPORTANT REQUIREMENTS 
//            VF topology 
//            Border FLags 
//        tri::UpdateTopology<SMesh>::VertexFace(sm);
//        tri::UpdateFlags<SMesh>::FaceBorderFromVF(sm);   
//
// Calcola la qualita' come distanza geodesica dal bordo della mesh.
// Robusta funziona anche per mesh non manifold.
// La qualita' memorizzata indica la distanza assoluta dal bordo della mesh.
// Nota prima del 13/11/03 in alcuni casi rari SPT andava in loop perche' poteva capitare
// che per approx numeriche ben strane pw->Q() > pv->Q()+d ma durante la memorizzazione 
// della nuova distanza essa rimanesse uguale a prima. Patchato rimettendo i vertici nello 
// heap solo se migliorano la distanza di un epsilon == 1/100000 della mesh diag.

/// \brief Compute, for each vertex of the mesh the geodesic distance from the border of the mesh itself.

/** 
It uses the classical Dijkstra Shortest Path Tree algorithm. 
The geodesic distance is approximated by allowing to walk only along edges of the mesh.
 
\warning VF topology, Per Vertex Quality and border flags already computed (see UpdateFlags::FaceBorderFromVF and UpdateTopology::VertexFace);

*/
static void VertexGeodesicFromBorder(MeshType &m)	// R1
{
	//Requirements
	assert(m.HasVFTopology());
	assert(m.HasPerVertexQuality());

  std::vector< VQualityHeap > heap;
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
	
 const ScalarType loc_eps=m.bbox.Diag()/ScalarType(100000);
 while( heap.size()!=0 )							// Shortest path tree
	{
		VertexPointer pv;
    std::pop_heap(heap.begin(),heap.end());
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
          std::push_heap(heap.begin(),heap.end());
				}
			}
		}
	}

	for(v=m.vert.begin();v!=m.vert.end();++v)
		if(v->Q()==-1)
			v->Q() = 0;
}


/** Assign to each vertex of the mesh a constant quality value. Useful for initialization.
*/
static void VertexConstant(MeshType &m, float q)
{
	VertexIterator vi;
	for(vi=m.vert.begin();vi!=m.vert.end();++vi) if(!(*vi).IsD()) 
		(*vi).Q()=q;
}

/** Clamp each vertex of the mesh with a range of values.
*/
static void VertexClamp(MeshType &m, float qmin, float qmax)
{
  VertexIterator vi;
  for(vi=m.vert.begin();vi!=m.vert.end();++vi) if(!(*vi).IsD())
    (*vi).Q()=std::min(qmax, std::max(qmin,(*vi).Q()));
}

/** Assign to each face of the mesh a constant quality value. Useful for initialization.
*/
static void FaceConstant(MeshType &m, float q)
{
	FaceIterator fi;
	for(fi=m.face.begin();fi!=m.face.end();++fi)		
		(*fi).Q()=q;
}


static void VertexFromGaussianCurvature(MeshType &m)
{ 
	VertexIterator vi;
	for(vi=m.vert.begin();vi!=m.vert.end();++vi) if(!(*vi).IsD()) 
		(*vi).Q() = (*vi).Kg();
}

static void VertexFromMeanCurvature(MeshType &m)
{ 
	VertexIterator vi;
	for(vi=m.vert.begin();vi!=m.vert.end();++vi) if(!(*vi).IsD()) 
		(*vi).Q() = (*vi).Kh();
}

/*
 *  Absolute Curvature                     
 *
 *                  2|H|                if K >= 0
 *  |k1| + |k2| = <
 *                  2 * sqrt(|H|^2-K)   otherwise
 *
 * defs and formulas taken from 
 *     
 * Improved curvature estimation for watershed segmentation of 3-dimensional meshes
 * S Pulla, A Razdan, G Farin - Arizona State University, Tech. Rep, 2001
 * and from
 * Optimizing 3D triangulations using discrete curvature analysis  
 * N Dyn, K Hormann, SJ Kim, D Levin - Mathematical Methods for Curves and Surfaces: Oslo, 2000 
 */

static void VertexFromAbsoluteCurvature(MeshType &m)
{ 
	VertexIterator vi;
	for(vi=m.vert.begin();vi!=m.vert.end();++vi) if(!(*vi).IsD()) 
	{
		if((*vi).Kg() >= 0) 
					(*vi).Q() = math::Abs( 2*(*vi).Kh() );
		else
		      (*vi).Q() = 2*math::Sqrt(math::Abs( (*vi).Kh()*(*vi).Kh() - (*vi).Kg())); 
	}
}

/*
 * RMS Curvature =   sqrt(4H^2-2K)
 * def and formula taken from 
 *     
 * Improved curvature estimation for watershed segmentation of 3-dimensional meshes
 * S Pulla, A Razdan, G Farin - Arizona State University, Tech. Rep, 2001
 */ 
static void VertexFromRMSCurvature(MeshType &m)
{ 
	VertexIterator vi;
	for(vi=m.vert.begin();vi!=m.vert.end();++vi) if(!(*vi).IsD()) 
		(*vi).Q() = math::Sqrt(math::Abs( 4*(*vi).Kh()*(*vi).Kh() - 2*(*vi).Kg())); 
}

}; //end class
} // end namespace
} // end namespace
#endif
