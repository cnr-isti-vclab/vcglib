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
Revision 1.1  2007/07/31 12:31:34  ganovelli
added


****************************************************************************/
#ifndef __VCG_TETRA_TRI_CONVERTER
#define __VCG_TETRA_TRI_CONVERTER
#include <map>
#include <vector>
#include <vcg/space/tetra3.h>
#include <vcg/complex/tetramesh/allocate.h>
namespace vcg {

/** Class Boundary.
This is class for exporting the boundary of a d simplicial complex as a d-1 simplicial complex
*/
class Boundary{
public:

///this function build a triangle mesh using the same pointers to the tetrahedral mesh vertex
template <class TetraContainer, class TriangleMeshType>
static void OfTetramesh(TetraContainer &tetra,TriangleMeshType &trim)
{
typedef typename TetraContainer::iterator TetraIterator;
typedef typename TetraContainer::value_type TetraVertexType;
typedef typename TriangleMeshType::FaceType FaceType;
typedef typename TriangleMeshType::VertexType TriangleVertexType;

TetraIterator ti;
TetraVertexType *v0;
TetraVertexType *v2;

trim.Clear();
for (ti=tetra.begin();ti<tetra.end();ti++)
{
if (!(ti->IsD()))
	{
		if ((ti->IsBorderF(0))||(ti->IsBorderF(1))||(ti->IsBorderF(2))||(ti->IsBorderF(3)))
			for (int i=0;i<4;i++)
				if (ti->IsBorderF(i))
					{
						FaceType f=FaceType();
						f.ClearFlags();
						f.V(0)=(TriangleVertexType*)ti->V(Tetra::VofF(i,0));
						f.V(1)=(TriangleVertexType*)ti->V(Tetra::VofF(i,1));
						f.V(2)=(TriangleVertexType*)ti->V(Tetra::VofF(i,2));
						trim.face.push_back(f);
					}
		}
	}
}

template <class TriVertexType >
struct InsertedV{

typedef typename TriVertexType::FaceType FaceType;
InsertedV( TriVertexType *_v, FaceType* _f,int _z):v(_v),f(_f),z(_z){}

TriVertexType *v;
FaceType* f;
int z;

const bool operator <(const InsertedV & o){
return (v<o.v);
}
const bool operator ==(const InsertedV & o){
return (v==o.v);
}
const bool operator !=(const InsertedV & o){
return (v!=o.v);
}
};


/// this function build a triangle mesh using new pointers to the tetrahedral mesh vertex
template <class TetraContainer, class TriangleMeshType>
static void OfTetrameshCopy(TetraContainer &tetra,TriangleMeshType &trim)
{
	typedef typename TetraContainer::iterator TetraIterator;
	typedef typename TetraContainer::value_type::VertexType TetraVertexType;
	typedef typename TriangleMeshType::FaceType FaceType;
	typedef typename TriangleMeshType::FaceIterator FaceIterator;
	typedef typename TriangleMeshType::VertexIterator TriVertexIterator;
	typedef typename TriangleMeshType::VertexType TriVertexType;

	vector<InsertedV<TriVertexType> > newVertices;
	typename vector<InsertedV<TriVertexType> >::iterator curr,next;
	TriVertexIterator vi;
	vector<TriVertexType*> redirect;

	OfTetramesh(tetra,trim);

	FaceIterator fi;

	for(fi = trim.face.begin(); fi != trim.face.end(); ++fi){
		newVertices.push_back(InsertedV<TriVertexType>( (*fi).V(0),&(*fi),0));
		newVertices.push_back(InsertedV<TriVertexType>( (*fi).V(1),&(*fi),1));
		newVertices.push_back(InsertedV<TriVertexType>( (*fi).V(2),&(*fi),2));
	}

sort(newVertices.begin(),newVertices.end());


int pos = 0;
curr = next = newVertices.begin();
while( next != newVertices.end()){
	if((*curr)!=(*next))
		pos++;
	(*next).f->V( (*next).z) = (TriVertexType*)pos;
	curr = next;
	next++;
}

typename vector<InsertedV<TriVertexType> >::iterator newE = unique(newVertices.begin(),newVertices.end());
for(curr = newVertices.begin();curr!= newE;++curr)
	trim.vert.push_back(*((*curr).v));

for(vi = trim.vert.begin(); vi != trim.vert.end(); ++vi)
	redirect.push_back(&(*vi));

for(fi = trim.face.begin(); fi != trim.face.end(); ++fi){
	(*fi).V(0) = redirect[(int)(*fi).V(0)];
	(*fi).V(1) = redirect[(int)(*fi).V(1)];
	(*fi).V(2) = redirect[(int)(*fi).V(2)];
	}
trim.vn = trim.vert.size();
trim.fn = trim.face.size();
}

};// End class
} // End namespace


#endif
