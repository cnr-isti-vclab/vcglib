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

#ifndef __VCGLIB_LENS_DISTORTION
#define __VCGLIB_LENS_DISTORTION

#include <cmath>
#include <vcg/space/point2.h>

namespace vcg
{

/** \brief The distortion of a real lens, as the rational Brown-Conrady model of OpenCV.

  vcg::Camera is a pinhole camera and has no lens model: a lens is kept beside it, and is
  meant to be applied where an image is read, so that the image a program works on is the
  one the pinhole camera describes and every projection stays a pinhole one.

  The model maps **normalized** camera coordinates, a point (x, y, z) in camera space seen
  as (x/z, y/z), from where a pinhole camera would see it to where the lens puts it. The
  coordinates follow OpenCV: x to the right, y **down** the image, z forward. So a pixel of
  the photo as shot and a pixel of the undistorted image are related through the same
  intrinsics K of the pinhole camera:

      photo pixel = K Distort(n),   undistorted pixel = K n

  and the model does not depend on the image resolution. With fx, fy, cx, cy in pixels, y
  down, n = ((u - cx) / fx, (v - cy) / fy).

  The parameters are OpenCV's, in its order k1 k2 p1 p2 k3 k4 k5 k6, the order in which
  OpenCV, COLMAP (SIMPLE_RADIAL, RADIAL, OPENCV, FULL_OPENCV) and most photogrammetry tools
  export them, so they copy across unchanged:

      r2 = x^2 + y^2
      radial = (1 + k1 r2 + k2 r2^2 + k3 r2^3) / (1 + k4 r2 + k5 r2^2 + k6 r2^3)
      x' = x radial + 2 p1 x y + p2 (r2 + 2 x^2)
      y' = y radial + p1 (r2 + 2 y^2) + 2 p2 x y

  All zero is the identity.
 */
template <class S>
class LensDistortion
{
public:
  typedef S ScalarType;

  S k1 = 0, k2 = 0, p1 = 0, p2 = 0, k3 = 0, k4 = 0, k5 = 0, k6 = 0;

  bool IsIdentity() const
  {
    return k1 == 0 && k2 == 0 && p1 == 0 && p2 == 0 && k3 == 0 && k4 == 0 && k5 == 0 && k6 == 0;
  }

  bool operator==(const LensDistortion &o) const
  {
    return k1 == o.k1 && k2 == o.k2 && p1 == o.p1 && p2 == o.p2
        && k3 == o.k3 && k4 == o.k4 && k5 == o.k5 && k6 == o.k6;
  }
  bool operator!=(const LensDistortion &o) const { return !(*this == o); }

  /// Where the lens puts the normalized point n: closed form. Undistorting an image needs
  /// only this direction, for each pixel of the result the pixel of the photo it comes from.
  Point2<S> Distort(const Point2<S> &n) const
  {
    const S x = n[0], y = n[1];
    const S r2 = x * x + y * y;
    const S radial = Radial(r2);
    return Point2<S>(x * radial + 2 * p1 * x * y + p2 * (r2 + 2 * x * x),
                     y * radial + p1 * (r2 + 2 * y * y) + 2 * p2 * x * y);
  }

  /// Whether the lens is invertible around n: false where it folds back on itself, which a
  /// calibrated model can do beyond the region its photos constrained (a negative k3 bends
  /// the radial curve down). Past a fold Distort lands inside the photo again, so an image
  /// undistorted there would repeat the photo mirrored; such pixels are better left empty.
  bool IsInvertibleAt(const Point2<S> &n) const
  {
    S jxx, jxy, jyy;
    Jacobian(n, jxx, jxy, jyy);
    return jxx * jyy - jxy * jxy > 0;
  }

  /** The normalized point n whose Distort is d, by Newton's method from n = d: the inverse
    has no closed form. Returns false, with the last iterate in n, when it does not converge
    within maxIter steps to a correction below tolerance: for a point the lens never reaches,
    beyond the fold of a model that has one (see IsInvertibleAt).
    Needed only for points measured on the photo as shot.
   */
  bool Undistort(const Point2<S> &d, Point2<S> &n, int maxIter = 20, S tolerance = S(1e-12)) const
  {
    n = d;
    for (int i = 0; i < maxIter; ++i) {
      S jxx, jxy, jyy;
      Jacobian(n, jxx, jxy, jyy);
      const S det = jxx * jyy - jxy * jxy;
      if (det == 0)
        return false;
      const Point2<S> f = Distort(n) - d;
      const Point2<S> step((jyy * f[0] - jxy * f[1]) / det, (jxx * f[1] - jxy * f[0]) / det);
      n -= step;
      if (step.SquaredNorm() < tolerance * tolerance)
        return true;
    }
    return false;
  }

private:
  // The Jacobian of Distort at n; it is symmetric
  void Jacobian(const Point2<S> &n, S &jxx, S &jxy, S &jyy) const
  {
    const S x = n[0], y = n[1];
    const S r2 = x * x + y * y;
    // radial and its derivative in r2, by the quotient rule
    const S num = 1 + r2 * (k1 + r2 * (k2 + r2 * k3));
    const S den = 1 + r2 * (k4 + r2 * (k5 + r2 * k6));
    const S dNum = k1 + r2 * (2 * k2 + r2 * 3 * k3);
    const S dDen = k4 + r2 * (2 * k5 + r2 * 3 * k6);
    const S radial = num / den;
    const S dRadial = (dNum * den - num * dDen) / (den * den);
    jxx = radial + 2 * x * x * dRadial + 2 * p1 * y + 6 * p2 * x;
    jyy = radial + 2 * y * y * dRadial + 6 * p1 * y + 2 * p2 * x;
    jxy = 2 * x * y * dRadial + 2 * p1 * x + 2 * p2 * y;
  }

  S Radial(S r2) const
  {
    return (1 + r2 * (k1 + r2 * (k2 + r2 * k3))) / (1 + r2 * (k4 + r2 * (k5 + r2 * k6)));
  }
};

typedef LensDistortion<float>  LensDistortionf;
typedef LensDistortion<double> LensDistortiond;

} // end namespace vcg

#endif // __VCGLIB_LENS_DISTORTION
