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
Revision 1.1  2004/05/13 15:51:40  turini
Initial Commit



****************************************************************************/
#ifndef __VCG_TRI_UPDATE_CLEAR_FLAGS
#define __VCG_TRI_UPDATE_CLEAR_FLAGS

namespace vcg {
namespace tri {

/** \addtogroup trimesh */
/*@{*/

/// Updating of Mesh Vertexes and Faces Flags.
/// This class is used to clear the vertex and face flags of a mesh.
template <class ComputeMeshType>
class UpdateFlags
{

public:
typedef ComputeMeshType MeshType; 
typedef typename MeshType::VertexType     VertexType;
typedef typename MeshType::VertexPointer  VertexPointer;
typedef typename MeshType::VertexIterator VertexIterator;
typedef typename MeshType::FaceType       FaceType;
typedef typename MeshType::FacePointer    FacePointer;
typedef typename MeshType::FaceIterator   FaceIterator;

/// Rimmette a zero tutti i flags della mesh.
static void Clear()
{
	FaceIterator fi;
	VertexIterator vi;
	for(fi=face.begin(); fi!=face.end(); ++fi)
		(*fi).Flags() = 0;
	for(vi=vert.begin(); vi!=vert.end(); ++vi)
		(*vi).Flags() = 0;
}

}; // end class

}	// End namespace
}	// End namespace


#endif
