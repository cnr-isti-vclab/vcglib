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
#ifndef VCG_VERTEX_INTERPOLATION_H
#define VCG_VERTEX_INTERPOLATION_H

#include <vcg/complex/complex.h>
#include <cmath>
#include <type_traits>

namespace vcg {
namespace tri {

/** \ingroup trimesh
 * \brief Attributes of a new vertex, interpolated from the vertices it lies between.
 *
 * Every algorithm that creates a vertex on an edge or inside a face (refinement,
 * clipping, cutting along a curve) has to give it the attributes of the surface there,
 * or it carries whatever the allocation left: a dark speck in a color, a spike in a
 * scalar field, a jump in the texture. This is the one place that does it, for the
 * per-vertex normal, color, quality and texture coordinate the mesh actually has.
 * The position is the caller's: it is often not the interpolated one (Butterfly, Loop).
 *
 * The normal is renormalized. Integer color channels are rounded, not truncated: with
 * weights that do not sum to exactly 1 in floating point, truncation turns a constant
 * color of 200 into 199.
 */
template <class MeshType>
class VertexInterpolator
{
public:
  typedef typename MeshType::VertexType VertexType;
  typedef typename MeshType::ScalarType ScalarType;

  /// \a nv gets the attributes at parameter \a t along the edge from \a v0 (t = 0) to \a v1 (t = 1).
  static void Lerp(const MeshType &m, VertexType &nv, const VertexType &v0, const VertexType &v1, ScalarType t)
  {
    const VertexType *v[2] = { &v0, &v1 };
    const ScalarType w[2] = { ScalarType(1) - t, t };
    Blend(m, nv, v, w, 2);
  }

  /// \a nv gets the attributes at barycentric coordinates \a ip in the triangle \a v0 \a v1 \a v2.
  static void Barycentric(const MeshType &m, VertexType &nv, const VertexType &v0, const VertexType &v1,
                          const VertexType &v2, const Point3<ScalarType> &ip)
  {
    const VertexType *v[3] = { &v0, &v1, &v2 };
    const ScalarType w[3] = { ip[0], ip[1], ip[2] };
    Blend(m, nv, v, w, 3);
  }

private:
  static void Blend(const MeshType &m, VertexType &nv, const VertexType *const v[], const ScalarType w[], int n)
  {
    if (tri::HasPerVertexNormal(m)) {
      typename VertexType::NormalType nrm = v[0]->cN() * w[0];
      for (int i = 1; i < n; ++i) nrm += v[i]->cN() * w[i];
      nv.N() = nrm.normalized();
    }
    if (tri::HasPerVertexColor(m)) {
      typedef typename VertexType::ColorType::ScalarType ChannelType;
      for (int c = 0; c < 4; ++c) {
        ScalarType sum = 0;
        for (int i = 0; i < n; ++i) sum += v[i]->cC()[c] * w[i];
        if constexpr (std::is_integral<ChannelType>::value)
          sum = std::round(std::min(ScalarType(std::numeric_limits<ChannelType>::max()), std::max(ScalarType(0), sum)));
        nv.C()[c] = ChannelType(sum);
      }
    }
    if (tri::HasPerVertexQuality(m)) {
      nv.Q() = v[0]->cQ() * w[0];
      for (int i = 1; i < n; ++i) nv.Q() += v[i]->cQ() * w[i];
    }
    if (tri::HasPerVertexTexCoord(m)) {
      nv.T().P() = v[0]->cT().P() * w[0];
      for (int i = 1; i < n; ++i) nv.T().P() += v[i]->cT().P() * w[i];
      nv.T().N() = v[0]->cT().N();
    }
  }
};

} // end namespace tri
} // end namespace vcg

#endif // VCG_VERTEX_INTERPOLATION_H
