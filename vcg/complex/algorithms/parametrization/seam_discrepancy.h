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

#ifndef VCG_COMPLEX_ALGORITHMS_PARAMETRIZATION_SEAM_DISCREPANCY_H
#define VCG_COMPLEX_ALGORITHMS_PARAMETRIZATION_SEAM_DISCREPANCY_H

#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/update/topology.h>
#include <vcg/space/point2.h>
#include <vcg/space/point3.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>
#include <vector>

namespace vcg {
namespace tri {

/// \brief How much the colors of a texture differ across its seams.
///
/// A texture seam is an edge that two faces share in 3D and place apart in UV: the surface is
/// continuous there, but its two sides are read from different places of the texture, so any
/// difference between those places is a visible line. This measures it, as the TexMetro tool
/// (Maggiordomo et al.) does: along the seam, the colors that the two sides give to the same
/// points are compared, and their mean distance is the discrepancy of the seam.
///
/// The measure is generic in what a texture is: the caller gives, as a template argument, an
/// object that says how large each texture is and what color a texel has, so that the library
/// needs no image code. It needs the faces to share the vertices of a seam: merge duplicated
/// vertices first (Clean::RemoveDuplicateVertex), or a seam looks like a border of the mesh.
template <class MeshType>
class SeamDiscrepancy
{
public:
    typedef typename MeshType::FaceType FaceType;
    typedef typename MeshType::FacePointer FacePointer;

    /// A seam, as the edge it is on in each of the two faces that share it.
    struct Seam
    {
        FacePointer f;
        int e;
        FacePointer g;
        int ge;
    };

    /// Which samples of the two sides are compared. The sample points of a side are taken along
    /// its edge, from the vertex the face lists first. Two consistently oriented faces list
    /// the shared edge in opposite directions, so PairByEdgeOrder, which compares the i-th
    /// sample of each, compares a point of one side with the mirrored point of the other, which
    /// is what TexMetro does; PairByPosition compares the same point on both sides.
    enum Pairing { PairByEdgeOrder, PairByPosition };

    struct Options
    {
        /// Samples along a seam for each texel of its length (the longer side's, rounded up).
        double samplesPerTexel = 5.0;
        Pairing pairing = PairByEdgeOrder;
    };

    /// Find the texture seams of m: the edges that exactly two faces share in 3D and that the
    /// per-wedge texture coordinates cut. Each is listed once. The face-face adjacency of m is
    /// rebuilt from its geometry and left that way, so a seam is not a border of it; a caller that
    /// then partitions the mesh into charts (UV_Utils::PerWedgeCharts) changes it, which does not
    /// matter to the seams listed, as they hold faces and not adjacency.
    static void FindSeams(MeshType &m, std::vector<Seam> &seams)
    {
        tri::RequireFFAdjacency(m);
        tri::RequirePerFaceWedgeTexCoord(m);
        UpdateTopology<MeshType>::FaceFace(m);
        seams.clear();
        for (typename MeshType::FaceIterator fi = m.face.begin(); fi != m.face.end(); ++fi)
        {
            if ((*fi).IsD()) continue;
            for (int i = 0; i < 3; i++)
            {
                if (face::IsBorder(*fi, i) || !face::IsManifold(*fi, i)) continue;
                FacePointer g = (*fi).FFp(i);
                if (g < &*fi) continue; // the other face lists it
                if (UpdateTopology<MeshType>::IsTexCoordSeam(*fi, i))
                {
                    Seam seam = { &*fi, i, g, (*fi).FFi(i) };
                    seams.push_back(seam);
                }
            }
        }
    }

    /// The discrepancy of a seam: the mean, over its samples, of the distance between the RGB
    /// colors that its two sides give to them, sampled bilinearly with the texture repeating.
    /// The texture coordinates are normalized, as on any mesh: a face reads the texture given by
    /// the index of its first wedge. `textures` must have
    ///     int Width(int texture) const;  int Height(int texture) const;
    ///     vcg::Point3d Texel(int texture, int x, int y) const;   // the color of a texel
    /// Returns NaN for a seam with no length, which has no sample.
    template <class Textures>
    static double Discrepancy(const Seam &seam, const Textures &textures, const Options &options = Options())
    {
        const FaceType &f = *seam.f;
        const FaceType &g = *seam.g;
        const int textureF = f.cWT(0).N();
        const int textureG = g.cWT(0).N();
        const Point2d scaleF(textures.Width(textureF), textures.Height(textureF));
        const Point2d scaleG(textures.Width(textureG), textures.Height(textureG));
        Point2d a1(f.cWT(seam.e).U(), f.cWT(seam.e).V());
        Point2d b1(f.cWT(f.Next(seam.e)).U(), f.cWT(f.Next(seam.e)).V());
        Point2d a2(g.cWT(seam.ge).U(), g.cWT(seam.ge).V());
        Point2d b2(g.cWT(g.Next(seam.ge)).U(), g.cWT(g.Next(seam.ge)).V());
        if (options.pairing == PairByPosition && g.cV(seam.ge) != f.cV(seam.e))
            std::swap(a2, b2); // g runs along the seam the other way: start where f starts
        const double length = std::max(Texels(b1 - a1, scaleF), Texels(b2 - a2, scaleG));
        const int samples = int(std::lround(options.samplesPerTexel * std::ceil(length)));
        if (samples <= 0)
            return std::numeric_limits<double>::quiet_NaN();
        double sum = 0;
        for (int i = 0; i < samples; ++i)
        {
            const double t = i / double(samples);
            const Point3d c1 = Sample(textures, textureF, Mix(a1, b1, t), scaleF);
            const Point3d c2 = Sample(textures, textureG, Mix(a2, b2, t), scaleG);
            sum += (c1 - c2).Norm();
        }
        return sum / samples;
    }

private:
    static Point2d Mix(const Point2d &a, const Point2d &b, double t) { return a * (1 - t) + b * t; }

    // The length in texels of a difference of normalized texture coordinates.
    static double Texels(const Point2d &delta, const Point2d &scale)
    {
        return std::hypot(delta.X() * scale.X(), delta.Y() * scale.Y());
    }

    // x modulo size, for any sign of x.
    static int Wrap(long long x, int size)
    {
        const int r = int(x % size);
        return r < 0 ? r + size : r;
    }

    // The color at a normalized texture coordinate, bilinear between the four texels around it,
    // texel centers at half-integers, row 0 at the top of the image and V growing upwards.
    template <class Textures>
    static Point3d Sample(const Textures &textures, int texture, const Point2d &uv, const Point2d &scale)
    {
        const double x = uv.X() * scale.X() - 0.5;
        const double y = scale.Y() - uv.Y() * scale.Y() - 0.5;
        const double fx = std::floor(x), fy = std::floor(y);
        const double wx = x - fx, wy = y - fy;
        const long long ix = (long long)fx, iy = (long long)fy;
        const int w = textures.Width(texture), h = textures.Height(texture);
        const int x0 = Wrap(ix, w), x1 = Wrap(ix + 1, w), y0 = Wrap(iy, h), y1 = Wrap(iy + 1, h);
        const Point3d top = textures.Texel(texture, x0, y0) * (1 - wx) + textures.Texel(texture, x1, y0) * wx;
        const Point3d bottom = textures.Texel(texture, x0, y1) * (1 - wx) + textures.Texel(texture, x1, y1) * wx;
        return top * (1 - wy) + bottom * wy;
    }
};

} // end Namespace tri
} // end Namespace vcg
#endif // VCG_COMPLEX_ALGORITHMS_PARAMETRIZATION_SEAM_DISCREPANCY_H
