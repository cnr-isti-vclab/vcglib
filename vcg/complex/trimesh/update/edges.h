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
Revision 1.2  2004/05/12 18:52:35  ganovelli
removed call to ComputeRT and put its body here

Revision 1.1  2004/05/12 10:39:45  ganovelli
created

****************************************************************************/
#ifndef __VCG_TRI_UPDATE_EDGES
#define __VCG_TRI_UPDATE_EDGES

#include <vcg/space/plane3.h>

namespace vcg {
namespace tri {

/** \addtogroup trimesh */
/*@{*/

	/// Management, updating and computation of per-vertex and per-face normals.
	/// This class is used to compute or update the normals that can be stored in the vertex or face component of a mesh.
	template <class ComputeMeshType>
	class UpdateEdges
	{

	public:
		typedef ComputeMeshType MeshType; 
		typedef typename MeshType::VertexType     VertexType;
		typedef typename MeshType::VertexPointer  VertexPointer;
		typedef typename MeshType::VertexIterator VertexIterator;
		typedef typename MeshType::FaceType       FaceType;
		typedef typename MeshType::FacePointer    FacePointer;
		typedef typename MeshType::FaceIterator   FaceIterator;
		typedef typename MeshType::FaceType::ScalarType     ScalarType;

		/// Calculates the vertex normal (if stored in the current face type)
		static void Box(ComputeMeshType &m)
		{
			m.bbox.SetNull();
			VertexIterator vi;
			for(vi=m.vert.begin();vi!=m.vert.end();++vi)
				if( !(*vi).IsD() )	m.bbox.Add((*vi).P());

		}

		static void Set(FaceType &f)
		{
			f.Flags() = f.Flags() & (~(FaceType::NORMX|FaceType::NORMY|FaceType::NORMZ));
		
			// Primo calcolo degli edges
			f.edge[0] = f.V(1)->P(); f.edge[0] -= f.V(0)->P();
			f.edge[1] = f.V(2)->P(); f.edge[1] -= f.V(1)->P();
			f.edge[2] = f.V(0)->P(); f.edge[2] -= f.V(2)->P();
			// Calcolo di plane
			f.plane.SetDirection(f.edge[0]^f.edge[1]);
			f.plane.SetOffset(f.plane.Direction() * f.V(0)->P());
			f.plane.Normalize();
			// Calcolo migliore proiezione
			ScalarType nx = math::Abs(f.plane.Direction()[0]);
			ScalarType ny = math::Abs(f.plane.Direction()[1]);
			ScalarType nz = math::Abs(f.plane.Direction()[2]);
			ScalarType d;
			if(nx>ny && nx>nz) { f.Flags() |= FaceType::NORMX; d = 1/f.plane.Direction()[0]; }
			else if(ny>nz)     { f.Flags() |= FaceType::NORMY; d = 1/f.plane.Direction()[1]; }
			else               { f.Flags() |= FaceType::NORMZ; d = 1/f.plane.Direction()[2]; }

			// Scalatura spigoli
			f.edge[0] *= d;
			f.edge[1] *= d;
			f.edge[2] *= d;
		}

		static void Set(ComputeMeshType &m)
		{
			FaceIterator f;

			for(f = m.face.begin(); f!=m.face.end(); ++f)
				if(!(*f).IsD())
					Set(*f);
		}

	}; // end class

}	// End namespace
}	// End namespace


#endif
