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

Revision 1.1  2004/22/04 14:32  pietroni
Initial commit


****************************************************************************/
#ifndef __VCG_TETRA_TRI_CONVERTER
#define __VCG_TETRA_TRI_CONVERTER
#include <map>
#include <vector>
#include<vcg/space/tetra3.h>
namespace vcg {
namespace tetra {


 /**  Class TriConverter.
 This is class for convetr tetrahedral mesh into triangle mesh
		@param STL_VERT_CONT (Template Parameter) Specifies the type of the vertices container any the vertex type.
				@param STL_TETRA_CONT (Template Parameter) Specifies the type of the tetrahedrons container any the tetrahedrons type.
 */
template  < class TETRA_MESH ,class TRI_MESH >
class TriConverter
{

public:

	/// The tetrahedral mesh type
	typedef TETRA_MESH TetraMeshType;
	/// The triangle mesh type
	typedef TRI_MESH TriangleMeshType;

  /// The tetrahedron type
  typedef typename TetraMeshType::TetraType TetraType;
  /// The triangle type
  typedef typename TriangleMeshType::FaceType FaceType;

	/// The vertex type of tetrahedreal Mesh
	typedef typename TetraMeshType::VertexType TetraVertexType;
	/// The vertex type of triangular Mesh
	typedef typename TriangleMeshType::VertexType TriVertexType;
  
  /// The type of vertex iterator on tetrahedral mesh
  typedef typename TetraMeshType::VertexIterator TetraVertexIterator;
  /// The type of vertex iterator on tetrahedral mesh
  typedef typename TriangleMeshType::VertexIterator TriVertexIterator;

	/// The type of tetra iterator
	typedef typename TetraMeshType::TetraIterator TetraIterator;
  /// The type of tetra iterator
	typedef typename TriangleMeshType::FaceIterator FaceIterator;
  
  /// The type of const tetra iterator
  typedef typename TetraMeshType::const_TetraIterator const_TetraIterator;
  /// The type of const face iterator
  typedef typename TriangleMeshType::ConstFaceIterator ConstFaceIterator;

  /// The type of tetrahedrons container
  typedef typename TetraMeshType::TetraContainer TetraContainer;

   /// The type of const vertex pointer of tetrahedral mesh
  typedef typename TetraMeshType::const_VertexPointer const_VertexPointer;


public:
  
/***********************************************/
/** @Convert to triangle-mesh functions
**/
//@{

///this function build a triangle mesh using the same pointers to the tetrahedral mesh vertex
void Convert(TetraContainer &tetra,TriangleMeshType &trim)
{
  TetraIterator ti;
  
  TetraVertexType *v0;
  TetraVertexType *v1;
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
        v0=ti->V(Tetra::VofF(i,0));
        v1=ti->V(Tetra::VofF(i,1));
        v2=ti->V(Tetra::VofF(i,2));
        FaceType f=FaceType();
        f.ClearFlags();
        f.V(0)=v0;
        f.V(1)=v1;
        f.V(2)=v2;
        trim.face.push_back(f);
       }
    }
  }

}

struct InsertedV{
	InsertedV(	TriVertexType *_v,
				FaceType* _f,	
				int _z):v(_v),f(_f),z(_z){}

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


///this function build a triangle mesh using new pointers to the tetrahedral mesh vertex

void ConvertCopy(TetraContainer &tetra,TriangleMeshType &trim)
{
	vector<InsertedV > newVertices;
	typename vector<InsertedV>::iterator curr,next;
	TriVertexIterator vi;
	vector<TriVertexType*> redirect;

  Convert(tetra,trim);

	FaceIterator fi;

	 for(fi = trim.face.begin(); fi != trim.face.end(); ++fi){
		newVertices.push_back(InsertedV( (*fi).V(0),&(*fi),0));
		newVertices.push_back(InsertedV( (*fi).V(1),&(*fi),1));
		newVertices.push_back(InsertedV( (*fi).V(2),&(*fi),2));
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

	typename vector<InsertedV>::iterator newE = unique(newVertices.begin(),newVertices.end());
	for(curr = newVertices.begin();curr!= newE;++curr)
		trim.vert.push_back(*((*curr).v));

	for(vi = trim.vert.begin(); vi != trim.vert.end(); ++vi)
		redirect.push_back(&(*vi));

	 for(fi = trim.face.begin(); fi != trim.face.end(); ++fi){
		(*fi).V(0) =	redirect[(int)(*fi).V(0)];
		(*fi).V(1) =  redirect[(int)(*fi).V(1)];
		(*fi).V(2) =  redirect[(int)(*fi).V(2)];
	}
	trim.vn = trim.vert.size();
	trim.fn = trim.face.size();
}

	};// End class
}	// End namespace
}	// End namespace


#endif
