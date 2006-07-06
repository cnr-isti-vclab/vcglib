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
  
  const bool operator <(const InsertedV & o)
  {
    return (v<o.v);
  }
  
  const bool operator ==(const InsertedV & o)
  {
    return (v==o.v);
  }
  
  const bool operator !=(const InsertedV & o)
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
  STL_CONT::iterator pfi;
  S_MESH_TYPE::VertexIterator vi;
  std::vector<S_MESH_TYPE::VertexPointer> redirect;
  
  for(pfi=subSet.begin(); pfi!=subSet.end(); ++pfi)
  {		
		assert(!(*pfi)->IsD());
		m.face.push_back(*(*pfi));
  }
  
  S_MESH_TYPE::FaceIterator fi;
  for(fi=m.face.begin(); fi!=m.face.end(); ++fi)
  {
    newVertices.push_back(InsertedV<S_MESH_TYPE>((*fi).V(0), &(*fi),0));
	newVertices.push_back(InsertedV<S_MESH_TYPE>((*fi).V(1), &(*fi),1));
	newVertices.push_back(InsertedV<S_MESH_TYPE>((*fi).V(2), &(*fi),2));
  }
  
  sort(newVertices.begin(), newVertices.end());
  
  std::vector< InsertedV<S_MESH_TYPE> >::iterator curr, next;
  int pos=0;
  curr=next=newVertices.begin();
  while(next!=newVertices.end())
  {
    if((*curr)!=(*next))
	  pos++;
	(*next).f->V((*next).z)=(S_MESH_TYPE::VertexPointer)pos;
	curr=next;
	next++;
  }
  
  std::vector< InsertedV<S_MESH_TYPE> >::iterator newE=unique(newVertices.begin(), newVertices.end());
  for(curr=newVertices.begin(); curr!=newE; ++curr)
    m.vert.push_back(*((*curr).v));
  
  for(vi=m.vert.begin(); vi!=m.vert.end(); ++vi)
    redirect.push_back(&(*vi));
  
  for(fi=m.face.begin(); fi!=m.face.end(); ++fi)
  {
    (*fi).V(0)=redirect[(int)(*fi).V(0)];
	(*fi).V(1)=redirect[(int)(*fi).V(1)];
	(*fi).V(2)=redirect[(int)(*fi).V(2)];
  }
  m.vn=(int)m.vert.size();
  m.fn=(int)m.face.size();
  vcg::tri::UpdateFlags<S_MESH_TYPE>::Clear(m);
}


} // End Namespace TriMesh
} // End Namespace vcg


#endif


