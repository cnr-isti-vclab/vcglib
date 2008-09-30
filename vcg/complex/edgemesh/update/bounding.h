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
#ifndef __VCG_EDGE_UPDATE_BOUNDING
#define __VCG_EDGE_UPDATE_BOUNDING

namespace vcg {
namespace edg {

/** \addtogroup edgemesh */
/*@{*/

/// Management, updating and computation of per-vertex and per-face normals.
/// This class is used to compute or update the normals that can be stored in the vertex or face component of a mesh.
template <class ComputeMeshType>
class UpdateBounding
{

public:
typedef ComputeMeshType MeshType; 
typedef typename MeshType::VertexType     VertexType;
typedef typename MeshType::VertexPointer  VertexPointer;
typedef typename MeshType::VertexIterator VertexIterator;
typedef typename MeshType::EdgeType       EdgeType;
typedef typename MeshType::EdgePointer    EdgePointer;
typedef typename MeshType::EdgeIterator   EdgeIterator;

/// Calculates the vertex normal (if stored in the current face type)
static void Box(ComputeMeshType &m)
{
  m.bbox.SetNull();
	VertexIterator vi;
	for(vi=m.vert.begin();vi!=m.vert.end();++vi)
			if( !(*vi).IsD() )	m.bbox.Add((*vi).P());

}


}; // end class

}	// End namespace
}	// End namespace


#endif
