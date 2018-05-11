/****************************************************************************
* VCGLib                                                            o o     *
* Visual and Computer Graphics Library                            o     o   *
*                                                                _   O  _   *
* Copyright(C) 2004-2016                                           \/)\/    *
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

#ifndef _VCG_TETRA_TOPOLOGY
#define _VCG_TETRA_TOPOLOGY

#include <vcg/complex/algorithms/update/topology.h> 

namespace vcg {
namespace tetrahedron {
/** \addtogroup tetrahedron */
/*@{*/

/** Return a boolean that indicate if the j-th face of the tetra is a border.
    @param j Index of the face
    @return true if j is an face of border, false otherwise
*/
template <class TetraType>
inline bool IsBorder(TetraType const & t,  const int j )
{
  if(TetraType::HasTTAdjacency())
      return t.cTTp(j)==&t;
  assert(0);
  return true;
}

template <class TetraMesh, class TriMesh>
inline void TriMeshFromBorder(TetraMesh & tetramesh, TriMesh & trimesh)
{
    typedef typename TriMesh::VertexPointer VertexPointer;
    typedef typename TriMesh::FacePointer     FacePointer;

    RequireTTAdjacency(tetramesh);
    tri::UpdateTopology<TetraMesh>::TetraTetra(tetramesh);

    trimesh.Clear();

    std::vector<VertexPointer> verts;
    std::vector<FacePointer>   faces;

    typedef typename TetraMesh::TetraType TetraType;
    ForEachTetra(tetramesh, [&] (TetraType & t) {
        for (int i = 0; i < 4; ++i)
            if (IsBorder(t, i))
            {
                verts.push_back(t.V(Tetra::VofF(i, 0)));
                verts.push_back(t.V(Tetra::VofF(i, 1)));
                verts.push_back(t.V(Tetra::VofF(i, 2)));
            }
    });

    typedef typename TriMesh::VertexIterator VertexIterator;
    typedef typename TriMesh::FaceIterator     FaceIterator;
    
    VertexIterator vi = tri::Allocator<TriMesh>::AddVertices(trimesh, verts.size());
    FaceIterator   fi = tri::Allocator<TriMesh>::AddFaces(trimesh, verts.size() / 3);
    
    for (int i = 0; i < verts.size(); i += 3)
    {
        fi->Alloc(3);
        
        vi->P() = verts[i + 0]->P();
        fi->V(0) = &*vi;
        ++vi;

        vi->P() = verts[i + 1]->P();
        fi->V(1) = &*vi;
        ++vi;

        vi->P() = verts[i + 2]->P();
        fi->V(2) = &*vi;
        ++vi;

        ++fi;
    }
}

}
}

#endif
