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
Revision 1.1  2005/07/06 08:02:27  cignoni
Initial commit


****************************************************************************/
#ifndef __VCG_VERTEX_UPDATE_POSITION
#define __VCG_VERTEX_UPDATE_POSITION

namespace vcg {
namespace vrt {

/** \addtogroup vertexmesh */
/*@{*/

/// This class is used to update vertex position according to a transformation matrix.
template <class VERTEX_CONTAINER>
class UpdatePositionBase
{

public:
	typedef typename VERTEX_CONTAINER::value_type::ScalarType     ScalarType;
	typedef typename VERTEX_CONTAINER::iterator VertexIterator;

	/// Multiply 
	static void Matrix(VERTEX_CONTAINER &vert, const Matrix44<ScalarType> &mm)
	{
		VertexIterator vi;
		for(vi= vert.begin();vi!= vert.end();++vi)
						if(!(*vi).IsD()) (*vi).P()=mm*(*vi).cP();
	}


}; // end class

template <class VertexMeshType>
class UpdatePosition: public UpdatePositionBase<typename VertexMeshType::VertexContainer> {
public:
	typedef typename VertexMeshType::VertexContainer VertexContainer;

	static void Matrix(VertexMeshType &vm, const Matrix44<ScalarType> &mm)
	{ UpdatePositionBase<typename VertexMeshType::VertexContainer>::Matrix(vm.vert,mm);
	}

};


}	// End namespace
}	// End namespace


#endif
