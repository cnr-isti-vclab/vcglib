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
Revision 1.1  2005/03/09 13:22:55  ganovelli
creation


****************************************************************************/
#ifndef __VCG_POINT_UPDATE_BOUNDING
#define __VCG_POINT_UPDATE_BOUNDING

#include <vcg/space/box3.h>

namespace vcg {
namespace vertex {

/** \addtogroup vertexmesh */
/*@{*/

/// Management, updating and computation of per-vertex and per-face normals.
/// This class is used to compute or update the normals that can be stored in the vertex or face component of a mesh.
template <class VERTEX_CONTAINER>
class UpdateBoundingBase
{

public:
typedef typename VERTEX_CONTAINER::value_type   VertexType;
typedef typename VERTEX_CONTAINER::value_type * VertexPointer;
typedef typename VERTEX_CONTAINER::iterator			VertexIterator;
typedef typename VERTEX_CONTAINER::value_type::ScalarType			ScalarType;

/// Calculates the vertex normal (if stored in the current face type)
static Box3<ScalarType> Box(VERTEX_CONTAINER &vert)
{
  Box3<ScalarType> res;res.SetNull();
	VertexIterator vi;
	for(vi= vert.begin();vi!= vert.end();++vi)
			if( !(*vi).IsD() )	res.Add((*vi).P());
	return res;
}

}; // end class

template <class VMType>
class UpdateBounding: public UpdateBoundingBase<typename VMType::VertexContainer> {
	public:
	typedef typename VMType::VertexContainer VertexContainer;

	static void Box(VMType &vm){
		vm.bbox = UpdateBoundingBase<typename VMType::VertexContainer>::Box(vm.vert);
		}

};

}	// End namespace
}	// End namespace


#endif
