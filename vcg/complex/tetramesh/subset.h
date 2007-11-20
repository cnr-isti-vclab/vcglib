

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
Revision 1.4  2005/02/08 14:36:20  turini
Warnings Correction

Revision 1.3  2004/05/17 08:22:45  turini
Minor Changes and Now Use STLContainer of Tetrahedron Pointers.

Revision 1.2  2004/05/14 15:51:47  turini
Adjusted VCG Style

Revision 1.1  2004/05/14 15:43:41  turini
Initial Commit



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
  typedef I_TETRAMESH_TYPE ITetraMeshType; 
  typedef typename ITetraMeshType::VertexPointer  VertexPointer;
  typedef typename ITetraMeshType::TetraPointer  TetraPointer;

  InsertedVT(VertexPointer _v, TetraPointer _t,	int _z)
    : v(_v), t(_t), z(_z)
  {}

  VertexPointer v;
  TetraPointer t;
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
@param subSet Container of tetrahedron pointers !!!
@param m destination mesh.
*/
template <class S_TETRAMESH_TYPE, class STL_CONT >
void SubSet(S_TETRAMESH_TYPE & m, STL_CONT & subSet)
{
  std::vector< InsertedVT<S_TETRAMESH_TYPE> > newVertices;
  typename STL_CONT::iterator pfi;
  newVertices.clear();
  
  for(pfi=subSet.begin(); pfi!=subSet.end(); ++pfi) 
    m.tetra.push_back(*(*pfi));
  
  typename S_TETRAMESH_TYPE::TetraIterator fi;
  for(fi=m.tetra.begin(); fi!=m.tetra.end(); ++fi)
  {
    newVertices.push_back(InsertedVT<S_TETRAMESH_TYPE>((*fi).V(0), &(*fi), 0));
	newVertices.push_back(InsertedVT<S_TETRAMESH_TYPE>((*fi).V(1), &(*fi), 1));
	newVertices.push_back(InsertedVT<S_TETRAMESH_TYPE>((*fi).V(2), &(*fi), 2));
	newVertices.push_back(InsertedVT<S_TETRAMESH_TYPE>((*fi).V(3), &(*fi), 3));
  }
  
  std::sort(newVertices.begin(), newVertices.end());
  
  typename std::vector< InsertedVT<S_TETRAMESH_TYPE> >::iterator curr,next;
  int pos=0;
  curr=next=newVertices.begin();
  while(next!=newVertices.end())
  {
    if((*curr)!=(*next))
	  pos++;
	(*next).t->V((*next).z)=(typename S_TETRAMESH_TYPE::VertexPointer)pos;
	curr=next;
	next++;
  }
  
  typename std::vector< InsertedVT<S_TETRAMESH_TYPE> >::iterator newE=std::unique(newVertices.begin(), newVertices.end());
  
  for(curr=newVertices.begin(); curr!=newE; ++curr)
    m.vert.push_back(*((*curr).v));
  
  for(fi=m.tetra.begin(); fi!=m.tetra.end(); ++fi)
  {
    (*fi).V(0)=&(m.vert[(int)(*fi).V(0)]);
	(*fi).V(1)=&(m.vert[(int)(*fi).V(1)]);
	(*fi).V(2)=&(m.vert[(int)(*fi).V(2)]);
	(*fi).V(3)=&(m.vert[(int)(*fi).V(3)]);
  }
  m.vn=(int)m.vert.size();
  m.tn=(int)m.tetra.size();
}


/*@}*/
}	// End namespace
}	// End namespace


#endif

