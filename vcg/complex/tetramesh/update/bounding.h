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
Revision 1.1  2004/05/04 11:15:13  pietroni
First working version!


****************************************************************************/
#ifndef __VCG_TETRA_UPDATE_BOUNDING
#define __VCG_TETRA_UPDATE_BOUNDING

namespace vcg {
namespace tetra {

/** \addtogroup tetramesh */
/*@{*/

/// Management, updating and computation of bonding box on a tetrahedral mesh

template <class ComputeMeshType>
class UpdateBounding
{

public:
typedef ComputeMeshType MeshType; 
typedef typename MeshType::VertexType     VertexType;
typedef typename MeshType::VertexPointer  VertexPointer;
typedef typename MeshType::VertexIterator VertexIterator;
typedef typename MeshType::TetraType      TetraType;
typedef typename MeshType::TetraPointer   TetraPointer;
typedef typename MeshType::TetraIterator  TetraIterator;

/// Calculates the limits of bounding box of tetrahedral mesh
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
