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

#ifndef VCG_COMPLEX_ALGORITHMS_PARAMETRIZATION_TEXEL_COVERAGE_H
#define VCG_COMPLEX_ALGORITHMS_PARAMETRIZATION_TEXEL_COVERAGE_H

#include <vcg/complex/complex.h>
#include <wrap/callback.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>

namespace vcg {
namespace tri {

/// \brief How the triangles of a UV layout cover the texels of a texture.
///
/// The faces that use one texture are rasterized at one resolution, texel centres against
/// triangles, and the texels are counted: how many are covered, how many by more than one
/// triangle (their texture is overwritten, and what the lower ones held is lost), how many
/// sit on the border of the layout, and how many a bilinear lookup would blend with a texel of
/// another chart. These are the texel measures of the TexMetro tool (Maggiordomo et al.), which
/// computed them by rendering on a GPU; here the rasterization is done on the CPU, so it gives
/// the same numbers on every machine and needs no graphics context.
///
/// A texel is covered by a triangle when its centre is inside it. A centre exactly on an edge
/// shared by two triangles goes to one of them only (the top-left rule), so a layout that
/// tiles a region covers each texel once. Vertices are first rounded to a grid of 2^SubpixelBits
/// positions per texel, as a GPU rasterizer does: the guaranteed minimum of OpenGL is 4 bits and
/// the usual hardware value 8, and the rare counts (overwritten, lost) change by up to about one
/// percent between them, so the number is part of the measure and is an option.
template <class MeshType>
class TexelCoverage
{
public:
    /// Which neighbours of a covered texel are looked at for a clash with another chart. TexMetro
    /// looks at the four diagonal ones only, which its test `i != row && j != col` happens to select.
    enum ClashRule { ClashOnDiagonals, ClashOnAllNeighbours };

    struct Options
    {
        /// Bits of position inside a texel the vertices are rounded to, from 1 to 8.
        int subpixelBits = 8;
        ClashRule clashRule = ClashOnDiagonals;
        /// The coverage of a texture is computed in bands of rows that fit in this many bytes
        /// (eight per texel), so a texture too large to hold at once is still counted exactly.
        std::size_t maxBufferBytes = std::size_t(1) << 30;
    };

    struct Counters
    {
        int width = 0;
        int height = 0;
        /// Fragments: a texel covered by n triangles counts n.
        int64_t totalTexels = 0;
        /// The covered texels, and the uncovered ones that touch a covered texel (3x3): the texels
        /// a bilinear lookup near the layout may read.
        int64_t totalTexelsBilinear = 0;
        /// The texels covered by more than one triangle.
        int64_t overwrittenTexels = 0;
        /// The fragments that overwrite another: a texel covered by n triangles loses n - 1.
        int64_t lostTexels = 0;
        /// The covered texels with a neighbour of another chart.
        int64_t texelClashes = 0;
        /// The covered texels on the border of the texture, or next to an uncovered one.
        int64_t boundaryTexels = 0;
    };

    /// Count the texels of a width x height texture covered by the faces that use `texture`
    /// (the texture index of their first wedge). The UV of a face are scaled by uScale and vScale
    /// to texels of this resolution. faceChart gives the chart of each face (see
    /// UV_Utils::PerWedgeCharts): it tells apart a texel overwritten by its own chart from one
    /// shared with another. Faces with a negative chart, with no area, or a UV that is not finite
    /// or lies more than two million texels from the origin cover nothing. Where several faces
    /// cover a texel the last in face order is the one its neighbours are compared with.
    /// cb, if given, is told how far the bands are.
    static Counters Compute(const MeshType &m, const std::vector<int> &faceChart, int texture,
                            int width, int height, double uScale, double vScale,
                            const Options &opt = Options(), CallBackPos *cb = nullptr)
    {
        tri::RequirePerFaceWedgeTexCoord(m);
        Counters total;
        total.width = width;
        total.height = height;
        if (width < 1 || height < 1)
            return total;

        const int bits = std::min(8, std::max(1, opt.subpixelBits));
        const std::size_t rowBytes = std::size_t(width) * 2 * sizeof(uint32_t);
        // A band is a run of rows; the row above and the one below it are kept as well, as the
        // texels of its first and last rows have neighbours there. A band is at least one row.
        const std::size_t fittingRows = opt.maxBufferBytes / rowBytes;
        const int bandRows = int(std::min<std::size_t>(
            std::size_t(height), fittingRows > 3 ? fittingRows - 2 : 1));
        std::vector<uint32_t> count, id;
        for (int y0 = 0; y0 < height; y0 += bandRows)
        {
            if (cb)
                (*cb)(100 * y0 / height, "Rasterizing texels");
            const int y1 = std::min(height, y0 + bandRows);
            const int r0 = std::max(0, y0 - 1), r1 = std::min(height, y1 + 1);
            count.assign(std::size_t(r1 - r0) * width, 0);
            id.assign(count.size(), 0);
            Rasterize(m, faceChart, texture, width, uScale, vScale, bits, r0, r1, count, id);
            Count(count, id, width, height, r0, y0, y1, opt.clashRule, total);
        }
        return total;
    }

private:
    typedef typename MeshType::FaceType FaceType;

    // Fill the rows [r0, r1) of the coverage count and of the last chart to cover each texel.
    static void Rasterize(const MeshType &m, const std::vector<int> &faceChart, int texture,
                          int width, double uScale, double vScale, int bits, int r0, int r1,
                          std::vector<uint32_t> &count, std::vector<uint32_t> &id)
    {
        // The edge functions multiply two coordinate differences; this bound keeps the product
        // of two of them, in units of 1/256 texel, inside 64 bits.
        const double maxCoord = double(1 << 21);
        const double grid = double(int64_t(1) << bits);
        const int64_t half = int64_t(1) << (bits - 1);
        for (std::size_t f = 0; f < m.face.size(); ++f)
        {
            const FaceType &face = m.face[f];
            if (face.IsD() || faceChart[f] < 0 || face.cWT(0).N() != texture)
                continue;
            int64_t x[3], y[3];
            bool usable = true;
            for (int c = 0; c < 3 && usable; ++c)
            {
                const double px = double(face.cWT(c).U()) * uScale;
                const double py = double(face.cWT(c).V()) * vScale;
                usable = std::abs(px) <= maxCoord && std::abs(py) <= maxCoord; // false for NaN
                if (usable)
                {
                    x[c] = int64_t(std::llround(px * grid));
                    y[c] = int64_t(std::llround(py * grid));
                }
            }
            if (!usable)
                continue;
            const int64_t area = (x[1] - x[0]) * (y[2] - y[0]) - (y[1] - y[0]) * (x[2] - x[0]);
            if (area == 0)
                continue;
            if (area < 0) // both windings are covered
            {
                std::swap(x[1], x[2]);
                std::swap(y[1], y[2]);
            }
            const int minY = int(std::max<int64_t>(r0, std::min({ y[0], y[1], y[2] }) >> bits));
            const int maxY = int(std::min<int64_t>(r1 - 1, std::max({ y[0], y[1], y[2] }) >> bits));
            const int minX = int(std::max<int64_t>(0, std::min({ x[0], x[1], x[2] }) >> bits));
            const int maxX = int(std::min<int64_t>(width - 1, std::max({ x[0], x[1], x[2] }) >> bits));
            for (int py = minY; py <= maxY; ++py)
            {
                for (int px = minX; px <= maxX; ++px)
                {
                    const int64_t cx = int64_t(2 * px + 1) * half, cy = int64_t(2 * py + 1) * half;
                    bool inside = true;
                    for (int e = 0; e < 3 && inside; ++e)
                    {
                        const int a = e, b = (e + 1) % 3;
                        const int64_t dx = x[b] - x[a], dy = y[b] - y[a];
                        const int64_t side = dx * (cy - y[a]) - dy * (cx - x[a]);
                        // On the edge itself the texel belongs to the triangle for which it is a
                        // top or left edge: of the two that share it, exactly one.
                        const bool topLeft = dy > 0 || (dy == 0 && dx < 0);
                        inside = side > 0 || (side == 0 && topLeft);
                    }
                    if (!inside)
                        continue;
                    const std::size_t k = std::size_t(py - r0) * width + px;
                    ++count[k];
                    id[k] = uint32_t(faceChart[f]) + 1;
                }
            }
        }
    }

    // Add to total the counters of the rows [y0, y1), which lie in the buffers that start at row r0.
    static void Count(const std::vector<uint32_t> &count, const std::vector<uint32_t> &id,
                      int width, int height, int r0, int y0, int y1, ClashRule rule, Counters &total)
    {
        for (int y = y0; y < y1; ++y)
        {
            for (int x = 0; x < width; ++x)
            {
                const uint32_t n = count[std::size_t(y - r0) * width + x];
                bool anyClear = (y == 0 || y == height - 1 || x == 0 || x == width - 1);
                bool anySet = false;
                for (int j = std::max(y - 1, 0); j < std::min(y + 2, height); ++j)
                {
                    for (int i = std::max(x - 1, 0); i < std::min(x + 2, width); ++i)
                    {
                        if (count[std::size_t(j - r0) * width + i] == 0)
                            anyClear = true;
                        else
                            anySet = true;
                    }
                }
                if (n == 0)
                {
                    if (anySet)
                        ++total.totalTexelsBilinear;
                    continue;
                }
                total.totalTexels += n;
                ++total.totalTexelsBilinear;
                if (n > 1)
                {
                    ++total.overwrittenTexels;
                    total.lostTexels += n - 1;
                }
                if (anyClear)
                    ++total.boundaryTexels;
                const uint32_t chart = id[std::size_t(y - r0) * width + x];
                bool clash = false;
                for (int j = std::max(y - 1, 0); j < std::min(y + 2, height) && !clash; ++j)
                {
                    for (int i = std::max(x - 1, 0); i < std::min(x + 2, width) && !clash; ++i)
                    {
                        const bool centre = (i == x && j == y);
                        const bool diagonal = (i != x && j != y);
                        if (centre || (rule == ClashOnDiagonals && !diagonal))
                            continue;
                        const uint32_t other = id[std::size_t(j - r0) * width + i];
                        clash = other != 0 && other != chart;
                    }
                }
                if (clash)
                    ++total.texelClashes;
            }
        }
    }
};

} // end Namespace tri
} // end Namespace vcg
#endif // VCG_COMPLEX_ALGORITHMS_PARAMETRIZATION_TEXEL_COVERAGE_H
