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

// marco: removed types FaceType, FacePointer, FaceIterator to allow the use of this method from vertex meshes


/****************************************************************************
  History

$Log: not supported by cvs2svn $
Revision 1.2  2004/09/15 11:16:27  ganovelli
changed P() to cP()

Revision 1.1  2004/04/05 11:56:13  cignoni
First working version!

Revision 1.2  2004/03/12 15:22:19  cignoni
Written some documentation and added to the trimes doxygen module

Revision 1.1  2004/03/05 10:59:24  cignoni
Changed name from plural to singular (normals->normal)

Revision 1.1  2004/03/04 00:05:50  cignoni
First working version!

Revision 1.1  2004/02/19 13:11:06  cignoni
Initial commit


****************************************************************************/
#ifndef __VCG_TRI_UPDATE_BOUNDING
#define __VCG_TRI_UPDATE_BOUNDING

namespace vcg {
namespace tri {

/// \ingroup trimesh 

/// \headerfile bounding.h vcg/complex/trimesh/update/bounding.h

/// \brief Management, updating and computation of per-vertex and per-face normals.
/** 
This class is used to compute or update the normals that can be stored in the vertex or face component of a mesh.
*/

template <class ComputeMeshType>
class UpdateBounding
{

public:
typedef ComputeMeshType MeshType; 
typedef typename MeshType::VertexType     VertexType;
typedef typename MeshType::VertexPointer  VertexPointer;
typedef typename MeshType::VertexIterator VertexIterator;

/// \brief Calculates the bounding box of the \code <ComputeMeshType> \endcode m

static void Box(ComputeMeshType &m)
{
	m.bbox.SetNull();
	VertexIterator vi;
	for(vi = m.vert.begin(); vi != m.vert.end(); ++vi)
		if( !(*vi).IsD() )	m.bbox.Add((*vi).cP());

}


}; // end class

}	// End namespace
}	// End namespace


#endif
