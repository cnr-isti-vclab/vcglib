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

$Log: position.h,v $
****************************************************************************/

#ifndef __VCG_TRI_UPDATE_TEXTURE
#define __VCG_TRI_UPDATE_TEXTURE

#include <vcg/space/plane.h>

namespace vcg {
namespace tri {

/** \addtogroup trimesh */
/*@{*/

/// This class is used to update vertex position according to a transformation matrix.
template <class ComputeMeshType>
class UpdateTexture
{

public:
typedef ComputeMeshType MeshType; 
typedef typename MeshType::ScalarType     ScalarType;
typedef typename MeshType::VertexType     VertexType;
typedef typename MeshType::VertexPointer  VertexPointer;
typedef typename MeshType::VertexIterator VertexIterator;
typedef typename MeshType::FaceType       FaceType;
typedef typename MeshType::FacePointer    FacePointer;
typedef typename MeshType::FaceIterator   FaceIterator;

static void WedgeTexFromPlanar(ComputeMeshType &m, Plane3<ScalarType> &pl)
{
	FaceIterator fi;
	for(fi=m.face.begin();vi!=m.face.end();++vi)
	        if(!(*fi).IsD()) 
							{
							
							}											
}

static void WedgeTexFromCamera(ComputeMeshType &m, Plane3<ScalarType> &pl)
{
	
}

}; // end class

}	// End namespace
}	// End namespace


#endif
