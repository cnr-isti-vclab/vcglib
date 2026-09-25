/****************************************************************************
* VCGLib                                                            o o     *
* Visual and Computer Graphics Library                            o     o   *
*                                                                _   O  _   *
* Copyright(C) 2004-2026                                           \/)\/    *
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

#ifndef __VCG_MC_CLOSED_ISOSURFACE
#define __VCG_MC_CLOSED_ISOSURFACE

#include <algorithm>
#include <cmath>
#include <limits>
#include <vcg/complex/algorithms/clean.h>
#include <vcg/complex/algorithms/create/marching_cubes.h>
#include <vcg/complex/algorithms/create/mc_trivial_walker.h>

namespace vcg {
namespace tri {

/** \brief Marching cubes over \a volume, with the surface closed where it runs into the
  volume's boundary.

  Plain extraction leaves the surface open wherever the solid it bounds -- the region below
  \a threshold -- reaches the edge of the sampled grid. The boundary it leaves lies on the
  grid's outer sample planes, but it is not a set of loops that one plane at a time could
  fill: where the solid reaches an edge or a corner of the grid, the boundary runs from one
  face onto the next.

  So the grid is closed the way marching cubes closes anything: it is surrounded by a layer
  of samples outside the solid, and the surface closes itself, around the grid's edges and
  corners too. Every closing vertex lies on a grid edge running from a boundary sample out
  into that layer; moved back onto its sample, the closure lies flat on the grid's faces, as
  a marching-squares tessellation of the part of each face inside the solid. Vertices moved
  onto the same sample are merged, and the faces left without area along the grid's edges
  are removed.

  The result is the part of the solid inside the grid, closed and manifold, with its caps
  exactly on the outermost sample planes. Unlike plain extraction it uses every sample:
  TrivialWalker never reads the last layer of the volume it walks, which here is padding.

  Filling each face with a planar tessellation instead (CapPlanarBoundary) gives fewer
  triangles, but the planar tessellator bridges holes in time cubic in the outline: on a
  256^3 noise field one face, 178 loops and 9461 edges, took 278 s and was then rejected.
  This closure is linear in the number of boundary cells and cannot fail.

  Works on a padded copy of \a volume.
 */
template <class MeshType, class VolumeType>
void BuildClosedIsosurface(MeshType &m, VolumeType &volume, float threshold, CallBackPos *cb = 0)
{
  typedef typename MeshType::ScalarType ScalarType;
  typedef typename MeshType::CoordType CoordType;
  typedef typename VolumeType::ScalarType VolumeScalar;

  // One layer of padding before the data and two after, the second only because the walker
  // skips the last layer. A pad sample mirrors the nearest data sample across the threshold,
  // so it is outside the solid while staying at the scale of the field.
  const Point3i n = volume.ISize();
  const Point3i p = n + Point3i(3, 3, 3);
  typename VolumeType::Box3x box;
  box.min = volume.bbox.min - volume.voxel;
  box.max = box.min;
  for (int a = 0; a < 3; ++a) box.max[a] += volume.voxel[a] * VolumeScalar(p[a]);
  VolumeType padded;
  padded.Init(p, box);
  const VolumeScalar thr = VolumeScalar(threshold);
  const VolumeScalar outside = std::nextafter(thr, std::numeric_limits<VolumeScalar>::max());
  for (int k = 0; k < p[2]; ++k)
    for (int j = 0; j < p[1]; ++j)
      for (int i = 0; i < p[0]; ++i) {
        const int di = std::clamp(i - 1, 0, n[0] - 1);
        const int dj = std::clamp(j - 1, 0, n[1] - 1);
        const int dk = std::clamp(k - 1, 0, n[2] - 1);
        const VolumeScalar v = volume.Val(di, dj, dk);
        if (di == i - 1 && dj == j - 1 && dk == k - 1)
          padded.Val(i, j, k) = v;
        else
          padded.Val(i, j, k) = std::max(thr + std::abs(v - thr), outside);
      }

  typedef TrivialWalker<MeshType, VolumeType> Walker;
  typedef MarchingCubes<MeshType, Walker> Extractor;
  Walker walker;
  Extractor mc(m, walker);
  walker.template BuildMesh<Extractor>(m, padded, mc, threshold, cb);

  // The data's first and last sample planes, computed as the walker computes a vertex on
  // them, so a closing vertex moved onto one coincides exactly with a vertex already there.
  CoordType lo, hi;
  padded.IPfToPf(CoordType(ScalarType(1), ScalarType(1), ScalarType(1)), lo);
  padded.IPfToPf(CoordType(ScalarType(n[0]), ScalarType(n[1]), ScalarType(n[2])), hi);
  for (typename MeshType::VertexIterator vi = m.vert.begin(); vi != m.vert.end(); ++vi)
    for (int a = 0; a < 3; ++a)
      (*vi).P()[a] = std::min(std::max((*vi).P()[a], lo[a]), hi[a]);
  Clean<MeshType>::RemoveDuplicateVertex(m);
  // The closure that went around a grid edge has collapsed onto it.
  Clean<MeshType>::RemoveFaceOutOfRangeArea(m, ScalarType(0));
  Clean<MeshType>::RemoveUnreferencedVertex(m);
  Allocator<MeshType>::CompactEveryVector(m);
}

} // end namespace tri
} // end namespace vcg
#endif // __VCG_MC_CLOSED_ISOSURFACE
