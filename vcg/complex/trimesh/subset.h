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



****************************************************************************/


#ifndef __VCGLIB_TRISUBSET
#define __VCGLIB_TRISUBSET


#include <vcg/Plane3.h>


namespace vcg {
namespace tri {


template < class STL_VERT_CONT, class STL_FACE_CONT > class Mesh;


template <class I_MESH_TYPE>
struct InsertedV
{
  InsertedV(I_MESH_TYPE::vertex_type *_v,
	        I_MESH_TYPE::face_pointer _f,	
			int _z):v(_v),f(_f),z(_z)
  {}
  
  I_MESH_TYPE::vertex_type *v;
  I_MESH_TYPE::face_pointer f;
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
//	  It assumes FF topology has been computed.
template <class S_MESH_TYPE, class STL_CONT>
void SubSet(STL_CONT & subSet, S_MESH_TYPE & m)
{
  vector< InsertedV<S_MESH_TYPE> > newVertices;
  STL_CONT::iterator pfi;
  S_MESH_TYPE::vertex_iterator vi;
  vector<S_MESH_TYPE::vertex_pointer> redirect;
  
  for(pfi=subSet.begin(); pfi!=subSet.end(); ++pfi)
    m.face.push_back(*(*pfi));
  
  S_MESH_TYPE::face_iterator fi;
  for(fi=m.face.begin(); fi!=m.face.end(); ++fi)
  {
    newVertices.push_back(InsertedV<S_MESH_TYPE>((*fi).V(0), &(*fi),0));
	newVertices.push_back(InsertedV<S_MESH_TYPE>((*fi).V(1), &(*fi),1));
	newVertices.push_back(InsertedV<S_MESH_TYPE>((*fi).V(2), &(*fi),2));
  }
  
  sort(newVertices.begin(), newVertices.end());
  
  vector< InsertedV<S_MESH_TYPE> >::iterator curr, next;
  int pos=0;
  curr=next=newVertices.begin();
  while(next!=newVertices.end())
  {
    if((*curr)!=(*next))
	  pos++;
	(*next).f->V((*next).z)=(S_MESH_TYPE::vertex_pointer)pos;
	curr=next;
	next++;
  }
  
  vector< InsertedV<S_MESH_TYPE> >::iterator newE=unique(newVertices.begin(), newVertices.end());
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
  m.vn = m.vert.size();
  m.fn = m.face.size();
  m.ClearFlags();
}


} // End Namespace TriMesh
} // End Namespace vcg


#endif


