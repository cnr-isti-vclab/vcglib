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
Revision 1.12  2007/05/31 09:39:55  cignoni
Small gcc compiling  issues

Revision 1.11  2006/07/06 12:30:32  ganovelli
misleading comment removed

Revision 1.10  2005/12/14 17:14:13  pietroni
added assert on deleted flag condition

Revision 1.9  2005/02/08 14:38:05  turini
Warnings Correction

Revision 1.8  2004/05/17 08:26:28  turini
Changed : Parameters Order As In  vcg::tetra::SubSet.

Revision 1.7  2004/05/17 07:58:16  turini
Minor Changes To Compile Even Without using namespace std.

Revision 1.6  2004/05/14 11:43:17  turini
Changed mesh ClearFlag call.

Revision 1.5  2004/05/13 09:59:20  turini
Added typedef typename in InsertedV

Revision 1.4  2004/05/07 10:06:46  turini
include Plane3 removed.

Revision 1.3  2004/05/07 09:35:09  turini
Added History Info

****************************************************************************/


#ifndef __VCGLIB_TRISUBSET
#define __VCGLIB_TRISUBSET

#include <vcg/complex/trimesh/update/flag.h>

namespace vcg {
namespace tri {


template <class I_MESH_TYPE>
struct InsertedV
{
  typedef I_MESH_TYPE IMeshType; 
  typedef typename IMeshType::VertexPointer  VertexPointer;
  typedef typename IMeshType::FacePointer  FacePointer;

  InsertedV(VertexPointer _v, FacePointer _f, int _z)
	: v(_v), f(_f), z(_z)
  {}
  
  VertexPointer v;
  FacePointer f;
  int z;
  
  bool operator < (const InsertedV & o) const
  {
    return (v<o.v);
  }
  
  bool operator ==(const InsertedV & o)
  {
    return (v==o.v);
  }
  
  bool operator !=(const InsertedV & o)
  {
    return (v!=o.v);
  }
};


//  This function build a nesh from a subset of faces of another.
//	  @param : subSet, stl vector of face poitners.
//             m,  output mesh mesh.


template <class S_MESH_TYPE, class STL_CONT>
void SubSet(S_MESH_TYPE & m, STL_CONT & subSet)
{
  std::vector< InsertedV<S_MESH_TYPE> > newVertices;
  typename STL_CONT::const_iterator pfi;
  typename S_MESH_TYPE::VertexIterator vi;
  typename S_MESH_TYPE::FaceIterator fi;
  typedef typename S_MESH_TYPE::VertexType S_VertexType;
  std::vector<typename S_MESH_TYPE::VertexPointer> redirect;

  
	fi = vcg::tri::Allocator<S_MESH_TYPE>::AddFaces(m,subSet.size());
  for(pfi=subSet.begin(); pfi!=subSet.end(); ++pfi)
  {		
		assert(!(*pfi)->IsD());
		(*fi).ImportData(**pfi);
		for(int ii = 0 ; ii < (*fi).VN(); ++ii)
			(*fi).V(ii) = (S_VertexType*)(void*)(*pfi)->V(ii);
		++fi;
  }
  

  for(fi=m.face.begin(); fi!=m.face.end(); ++fi)
	 for(int ii = 0 ; ii < (*fi).VN(); ++ii)
		newVertices.push_back(InsertedV<S_MESH_TYPE>((*fi).V(ii), &(*fi),ii));
  
  sort(newVertices.begin(), newVertices.end());
  
  typename std::vector< InsertedV<S_MESH_TYPE> >::iterator curr, next;
  int pos=0;
  curr=next=newVertices.begin();
  while(next!=newVertices.end())
  {
    if((*curr)!=(*next))
	  pos++;
	(*next).f->V((*next).z)=(typename S_MESH_TYPE::VertexPointer)pos;
	curr=next;
	next++;
  }
  
  typename std::vector< InsertedV<S_MESH_TYPE> >::iterator newE=unique(newVertices.begin(), newVertices.end());
	
	vi = vcg::tri::Allocator<S_MESH_TYPE>::AddVertices(m,newE-newVertices.begin());
	for(curr=newVertices.begin(); curr!=newE; ++curr,++vi)
		(*vi).ImportData(*((*curr).v));

	for(vi=m.vert.begin(); vi!=m.vert.end(); ++vi)
		redirect.push_back(&(*vi));

	for(fi=m.face.begin(); fi!=m.face.end(); ++fi)
		{		
		  for(int ii = 0 ; ii < (*fi).VN(); ++ii)
			(*fi).V(ii)=redirect[(size_t)(*fi).V(ii)];
		}
	m.vn=(int)m.vert.size();
	m.fn=(int)m.face.size();
}


} // End Namespace TriMesh
} // End Namespace vcg


#endif


