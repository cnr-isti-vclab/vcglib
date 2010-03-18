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
Revision 1.4  2006/05/16 21:36:54  cignoni
Removed unused box function and rewrote initial comment.

Revision 1.3  2006/05/15 13:12:36  pietroni
Updating of edge values id divided into 2 functions ( the first one update only a face...) added also resetting of edges flags.. (see first line of Set function)

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

/// \ingroup trimesh 

	/// \headerfile edges.h vcg/complex/trimesh/update/edges.h

	/// \brief This class is used to compute or update the precomputed data used to efficiently compute point-face distances.
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
    typedef typename MeshType::FaceType::CoordType::ScalarType     ScalarType;

		static void Set(FaceType &f)
		{
			f.Flags() = f.Flags() & (~(FaceType::NORMX|FaceType::NORMY|FaceType::NORMZ));
		
			// Primo calcolo degli edges
			f.Edge(0) = f.V(1)->P(); f.Edge(0) -= f.V(0)->P();
			f.Edge(1) = f.V(2)->P(); f.Edge(1) -= f.V(1)->P();
			f.Edge(2) = f.V(0)->P(); f.Edge(2) -= f.V(2)->P();
			// Calcolo di plane
			f.Plane().SetDirection(f.Edge(0)^f.Edge(1));
			f.Plane().SetOffset(f.Plane().Direction().dot(f.V(0)->P()));
			f.Plane().Normalize();
			// Calcolo migliore proiezione
			ScalarType nx = math::Abs(f.Plane().Direction()[0]);
			ScalarType ny = math::Abs(f.Plane().Direction()[1]);
			ScalarType nz = math::Abs(f.Plane().Direction()[2]);
			ScalarType d;
			if(nx>ny && nx>nz) { f.Flags() |= FaceType::NORMX; d = 1/f.Plane().Direction()[0]; }
			else if(ny>nz)     { f.Flags() |= FaceType::NORMY; d = 1/f.Plane().Direction()[1]; }
			else               { f.Flags() |= FaceType::NORMZ; d = 1/f.Plane().Direction()[2]; }

			// Scalatura spigoli
			f.Edge(0)*=d;
			f.Edge(1)*=d;
			f.Edge(2)*=d;
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
