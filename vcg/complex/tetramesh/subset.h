

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


#ifndef __VCGLIB_TETRASUBSET
#define __VCGLIB_TETRASUBSET


namespace vcg {
namespace tetra {

/** \addtogroup tetramesh */
/*@{*/


template <class I_TETRAMESH_TYPE>
struct InsertedVT
{
  InsertedVT(I_TETRAMESH_TYPE::VertexType *_v, I_TETRAMESH_TYPE::TetraType *_t,	int _z)
    : v(_v), t(_t), z(_z)
  {}

  I_TETRAMESH_TYPE::VertexType *v;
  I_TETRAMESH_TYPE::TetraType *t;
  int z;
  
  const bool operator <(const InsertedVT & o)
  {
    return (v<o.v);
  }
  
  const bool operator ==(const InsertedVT & o)
  {
    return (v==o.v);
  }
  
  const bool operator !=(const InsertedVT & o)
  {
    return (v!=o.v);
  }
};


/** Create a copy of the mesh with tetrahedron that are into the templated container
@param ST_CONT (Template Parameter) Specifies the type of the container of tetrahedron.
@param subSet Container of tetrahedron.
@param m destination mesh.
*/
/*
template <class S_MESH_TYPE, class STL_CONT>
void SubSet(STL_CONT & subSet, S_MESH_TYPE & m)
*/
template <class S_TETRAMESH_TYPE, class STL_CONT >
void SubSet(TM_TYPE & m, STL_CONT & subSet, )
{
  std::vector<InsertedVT> newVertices;
  STL_CONT::iterator pfi;
  newVertices.clear();
  
  for(pfi=subSet.begin(); pfi!=subSet.end(); ++pfi) 
    m.tetra.push_back((*pfi));
  
  S_TETRAMESH_TYPE::TetraIterator fi;
  for(fi=m.tetra.begin(); fi!=m.tetra.end(); ++fi)
  {
    newVertices.push_back(InsertedVT((*fi).V(0), &(*fi), 0));
	newVertices.push_back(InsertedVT((*fi).V(1), &(*fi), 1));
	newVertices.push_back(InsertedVT((*fi).V(2), &(*fi), 2));
	newVertices.push_back(InsertedVT((*fi).V(3), &(*fi), 3));
  }
  
  std::sort(newVertices.begin(), newVertices.end());
  
  std::vector<InsertedVT>::iterator curr,next;
  int pos=0;
  curr=next=newVertices.begin();
  while(next!=newVertices.end())
  {
    if((*curr)!=(*next))
	  pos++;
	(*next).t->V((*next).z)=(VertexType*)pos;
	curr=next;
	next++;
  }
  
  std::vector<InsertedVT >::iterator newE=std::unique(newVertices.begin(), newVertices.end());
  
  for(curr=newVertices.begin(); curr!=newE; ++curr)
    m.vert.push_back(*((*curr).v));
  
  for(fi=m.tetra.begin(); fi!=m.tetra.end(); ++fi)
  {
    (*fi).V(0)=&(m.vert[(int)(*fi).V(0)]);
	(*fi).V(1)=&(m.vert[(int)(*fi).V(1)]);
	(*fi).V(2)=&(m.vert[(int)(*fi).V(2)]);
	(*fi).V(3)=&(m.vert[(int)(*fi).V(3)]);
  }
  m.vn=m.vert.size();
  m.tn=m.tetra.size();
}


/*@}*/
}	// End namespace
}	// End namespace


#endif

