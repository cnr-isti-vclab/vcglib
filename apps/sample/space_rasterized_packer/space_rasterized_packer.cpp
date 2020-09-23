/****************************************************************************
* VCGLib                                                            o o     *
* Visual and Computer Graphics Library                            o     o   *
*                                                                _   O  _   *
* Copyright(C) 2004-2019                                           \/)\/    *
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
#include <iostream>
#include <vector>
#include <cassert>
#include <string>

#include <vcg/space/point2.h>
#include <vcg/space/point3.h>
#include <vcg/complex/algorithms/outline_support.h>
#include <wrap/qt/outline2_rasterizer.h>
#include <vcg/space/rasterized_outline2_packer.h>

#include <wrap/qt/Outline2ToQImage.h>

// This sample shows how to pack a sequence of outlines at the given scale into
// a number of containers not known in advance using the PackBestEffort()
// function

typedef std::vector<vcg::Point2f> Outline2f;
typedef vcg::RasterizedOutline2Packer<float,QtOutline2Rasterizer> RasterizedPacker;
typedef RasterizedPacker::Parameters PackingParams;


void FillOutlineVec(std::vector<Outline2f>& outlines)
{
    vcg::tri::OutlineUtil<float>::BuildRandomOutlineVec(30, outlines);
    for (auto& outline : outlines)
        for (auto& p : outline)
            p *= 200.0;
}

int main(int argc, char **argv)
{
    // the vector of outlines to pack
    std::vector<Outline2f> outlines;
    FillOutlineVec(outlines);

    // containerIndices maps each outline index to the container in which is packed
    std::vector<int> containerIndices(outlines.size(), -1); // -1 means not packed to any container

    // packingTransforms maps each outline index to its transform */
    std::vector<vcg::Similarity2f> packingTransforms(outlines.size(), vcg::Similarity2f{});

    /* size of the packing area */
    const vcg::Point2i grid_size(400, 600);

    PackingParams params;
    params.costFunction  = PackingParams::LowestHorizon;
    params.doubleHorizon = true;
    params.innerHorizon  = true;   /* use inner horizons to pack charts in between gaps */
    params.permutations  = false;  /* do not use permutations (they are ignored by PackBestEffort() */
    params.rotationNum   = 16;     /* 16 rotations per chart */
    params.gutterWidth   = 2;      /* 2 pixels gutter */
    params.minmax        = false;  /* do not combine costs of the two horizons when evaluating placements */

    int totPacked = 0;
    int nc = 0; /* current container index */
    while (true) {

        std::cout << "Packing into container " << nc << std::endl;

        // build a vector with the outlines not yet packed, and a corresponding vector
        // of indices into the 'outlines' vector
        std::vector<int> outlineIndex_iter;
        std::vector<Outline2f> outlines_iter;
        for (unsigned i = 0; i < containerIndices.size(); ++i) {
            if (containerIndices[i] == -1) {
                outlineIndex_iter.push_back(i);
                outlines_iter.push_back(outlines[i]);
            }
        }

        std::vector<vcg::Similarity2f> transforms;
        std::vector<int> polyToContainer;
        int numPacked = RasterizedPacker::PackBestEffort(outlines_iter, {grid_size}, transforms, polyToContainer, params);

        totPacked += numPacked;

        if (numPacked == 0) {
            std::cerr << "Failed to pack any outline at the current iteration. Stopping." << std::endl;
            std::exit(-1);
        } else {
            for (unsigned i = 0; i < outlines_iter.size(); ++i) {
                if (polyToContainer[i] != -1) {
                    assert(polyToContainer[i] == 0); // We only use a single container in this example
                    int outline_i = outlineIndex_iter[i];
                    assert(containerIndices[outline_i] == -1);
                    containerIndices[outline_i] = nc;
                    packingTransforms[outline_i] = transforms[i];
                }
            }
        }

        if (totPacked == outlines.size())
            break;
        else
            nc++;

        if (nc > 10)
            std::cerr << "Warning: packing into more than 10 containers! (nc = " << nc << ")" << std::endl;
    }

    /* for each container, generate an image */
    for (int i = 0; i < nc + 1; ++i) {
        std::vector<Outline2f> outlines_iter;
        std::vector<vcg::Similarity2f> transforms_iter;
        for (unsigned k = 0; k < containerIndices.size(); ++k) {
            assert(containerIndices[k] != -1 && "Some outlines were not packed");
            if (containerIndices[k] == i) {
                outlines_iter.push_back(outlines[k]);
                transforms_iter.push_back(packingTransforms[k]);
            }
        }
        std::string filename = std::string("container_") + std::to_string(i) + std::string(".png");
        Outline2Dumper::Param pp;
        pp.width = grid_size.X();
        pp.height = grid_size.Y();
        pp.fill = true;
        pp.randomColor = true;

        Outline2Dumper::dumpOutline2VecPNG(filename.c_str(), outlines_iter, transforms_iter, pp);
    }
}


