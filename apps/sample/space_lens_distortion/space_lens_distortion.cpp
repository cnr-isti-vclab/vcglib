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

// A lens kept beside a pinhole camera: how a photo calibrated with a real lens (by OpenCV,
// COLMAP and the like) is related to the image the pinhole camera describes.

#include <cstdio>
#include <vcg/math/lens_distortion.h>

int main()
{
  // A strong barrel lens, with the parameters as a calibration tool exports them
  vcg::LensDistortiond lens;
  lens.k1 = -0.28;
  lens.k2 = 0.07;
  lens.p1 = 0.0004;
  lens.p2 = -0.0002;

  // The pinhole intrinsics of a 4000x3000 photo, in pixels, y down the image
  const double fx = 3200, fy = 3200, cx = 2000, cy = 1500;

  // A pixel of the undistorted image: where in the photo as shot is its colour?
  // Undistorting a whole image is this, for every pixel of the result.
  const vcg::Point2d u(3800, 2800);
  const vcg::Point2d n((u[0] - cx) / fx, (u[1] - cy) / fy);
  const vcg::Point2d d = lens.Distort(n);
  printf("undistorted pixel (%.1f, %.1f) comes from photo pixel (%.3f, %.3f)\n",
         u[0], u[1], d[0] * fx + cx, d[1] * fy + cy);

  // The other way, for a point picked on the photo as shot: iterative
  vcg::Point2d back;
  if (lens.Undistort(d, back))
    printf("and back to (%.3f, %.3f)\n", back[0] * fx + cx, back[1] * fy + cy);
  return 0;
}
