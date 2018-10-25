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
#ifndef TETRAMESH_SUPPORT_H
#define TETRAMESH_SUPPORT_H

namespace vcg {
namespace tri {

template <class TetraMesh, class TriMesh>
inline void CreateTriMeshFromTTBorder(TetraMesh & tetramesh, TriMesh & trimesh)
{
    RequireTTAdjacency(tetramesh);
    tri::UpdateTopology<TetraMesh>::TetraTetra(tetramesh);
    trimesh.Clear();

    typedef typename TetraMesh::TetraType TetraType;
    ForEachTetra(tetramesh, [&] (TetraType & t) {
        for (int i = 0; i < 4; ++i)
            if (IsTTBorder(t, i))
            {
                tri::Allocator<TriMesh>::AddFace(trimesh, 
                                                 t.V(Tetra::VofF(i, 0)),
                                                 t.V(Tetra::VofF(i, 1)),
                                                 t.V(Tetra::VofF(i, 2)));
            }
    });

	vcg::tri::Clean<TriMesh>::RemoveDuplicateVertex(trimesh);
	vcg::tri::Allocator<TriMesh>::CompactEveryVector(trimesh);
}

} // end namespace tri
} // end namespace vcg

#endif // EXTRUDE_H
