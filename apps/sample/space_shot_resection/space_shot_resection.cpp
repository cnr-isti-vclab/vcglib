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

// Camera resection: correcting a camera from a few points whose position is known and whose
// pixels in its image have been picked.

#include <cstdio>
#include <vector>
#include <vcg/math/shot_resection.h>

int main()
{
  // A camera for a 4000x3000 image, 3200 pixels of focal length, looking at the origin
  vcg::Shotd truth;
  truth.Intrinsics.ViewportPx = vcg::Point2i(4000, 3000);
  truth.Intrinsics.CenterPx = vcg::Point2d(2000, 1500);
  truth.Intrinsics.PixelSizeMm = vcg::Point2d(0.01, 0.01);
  truth.Intrinsics.FocalMm = 32;
  truth.LookAt(3, 2, 10, 0, 0, 0, 0, 1, 0);

  // Six points on the object, and where that camera sees them
  const std::vector<vcg::Point3d> points = {{-1, -1, 0}, {1, -1, 0.5}, {1, 1, 0}, {-1, 1, -0.5}, {0, 0, 1}, {0.5, -0.3, -1}};
  std::vector<vcg::Point2d> pixels;
  for (const vcg::Point3d &p : points)
    pixels.push_back(truth.Project(p));

  // A guess that is off in where it stands, where it looks, and its focal length
  vcg::Shotd shot = truth;
  shot.LookAt(3.4, 1.7, 9.5, 0.1, -0.1, 0, 0, 1, 0);
  shot.Intrinsics.FocalMm = 30;
  printf("before: %.3f pixels off\n", vcg::ShotResectiond::RmsError(shot, points, pixels));

  const double rms = vcg::ShotResectiond::Fit(shot, points, pixels, true);
  printf("after:  %.2e pixels off, focal %.4f mm\n", rms, shot.Intrinsics.FocalMm);
  return 0;
}
