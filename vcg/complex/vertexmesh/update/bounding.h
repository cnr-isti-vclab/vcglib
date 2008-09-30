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
Revision 1.2  2006/06/29 13:02:38  ganovelli
agiunta UpdateBoundingBase, superclasse di UpdateBounding, templated sul container di vertici.

Revision 1.1  2005/03/09 13:22:55  ganovelli
creation


****************************************************************************/
#ifndef __VCG_POINT_UPDATE_BOUNDING
#define __VCG_POINT_UPDATE_BOUNDING

#include <vcg/space/box3.h>

namespace vcg {
namespace vrt {

/** \addtogroup vertexmesh */
/*@{*/

template <class VERTEX_CONTAINER>
class UpdateBoundingBase
{

public:
	typedef typename VERTEX_CONTAINER::value_type   VertexType;
	typedef typename VERTEX_CONTAINER::value_type*  VertexPointer;
	typedef typename VERTEX_CONTAINER::iterator		VertexIterator;
	typedef typename VERTEX_CONTAINER::value_type::ScalarType	ScalarType;

	static Box3<ScalarType> Box(VERTEX_CONTAINER &vert)
	{
		Box3<ScalarType> res; res.SetNull();
		VertexIterator vi;
		for(vi= vert.begin();vi!= vert.end();++vi)
			if( !(*vi).IsD() )	res.Add((*vi).P());
		return res;
	}

}; // end class UpdateBoundingBase

template <class VMType>
class UpdateBounding: public UpdateBoundingBase<typename VMType::VertexContainer> 
{
	public:
	static void Box(VMType &vm)
	{
		vm.bbox = UpdateBoundingBase<VMType::VertexContainer>::Box(vm.vert);
	}
};

}	// End namespace vert
}	// End namespace vcg


#endif
