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

#ifndef __VCGLIB_SHOT_RESECTION
#define __VCGLIB_SHOT_RESECTION

#include <cmath>
#include <stdexcept>
#include <vector>
#include <vcg/math/shot.h>
#include <wrap/newuoa/include/newuoa.h>

namespace vcg
{

/** \brief Camera resection: where a perspective Shot stands and where it looks, and optionally
  its focal length, from points of known position and the pixels where its image shows them.

  The fit starts from the shot it is given and refines it, minimizing the root mean square
  distance between Shot::Project of each point and its pixel, so the pixels are in Project's
  own convention: viewport pixels, y up. The shot has to be roughly right to start with, as
  for any local method: near enough for the points to be in front of it and on the right
  side of each other.

  Nothing else of the shot changes: the principal point, the pixel size and the viewport are
  the calibration's, and a fit that moved them to explain a few picked points would undo it.

  The parameters are a rotation and a translation of the camera **in its own frame**, about
  the centroid of the points, plus a change of FocalMm. Each is scaled so that a unit step
  moves the projected points by about one pixel, which makes a single trust radius, in
  pixels, mean the same for all of them; the scales are measured on the points themselves,
  so the fit is deterministic. It is minimized by NEWUOA, Powell's derivative-free method,
  in rounds: each round measures the scales again where the previous one ended, so that a
  start far from the solution, where the first scales no longer hold, still converges.

  Three points determine the six parameters of where the camera stands and looks, four with
  the focal length; more are averaged in the least squares sense. Throws
  std::invalid_argument for fewer, for an orthographic shot, or for points all in one place.
 */
template <class S>
class ShotResection
{
public:
  typedef Shot<S> ShotType;

  /// The root mean square distance, in pixels, between where shot projects the points and
  /// the pixels given for them.
  static S RmsError(const ShotType &shot, const std::vector<Point3<S>> &points,
                    const std::vector<Point2<S>> &pixels)
  {
    S sum = 0;
    for (size_t i = 0; i < points.size(); ++i)
      sum += (shot.Project(points[i]) - pixels[i]).SquaredNorm();
    return points.empty() ? S(0) : std::sqrt(sum / S(points.size()));
  }

  /// Refines shot so that it projects the points onto their pixels; returns the root mean
  /// square error left, in pixels. Rounds stop when one improves the error by less than
  /// tolerance pixels, or after maxRounds.
  static S Fit(ShotType &shot, const std::vector<Point3<S>> &points,
               const std::vector<Point2<S>> &pixels, bool fitFocal = false,
               S tolerance = S(1e-9), int maxRounds = 20)
  {
    if (points.size() != pixels.size())
      throw std::invalid_argument("ShotResection: as many pixels as points are needed");
    if (points.size() < (fitFocal ? 4u : 3u))
      throw std::invalid_argument(fitFocal ? "ShotResection: at least 4 points are needed to fit the focal length"
                                           : "ShotResection: at least 3 points are needed");
    if (shot.Intrinsics.cameraType != Camera<S>::PERSPECTIVE)
      throw std::invalid_argument("ShotResection: only a perspective shot can be fitted");

    Point3<S> centroid(0, 0, 0);
    for (const Point3<S> &p : points)
      centroid += p;
    centroid /= S(points.size());
    S spread = 0;
    for (const Point3<S> &p : points)
      spread += (p - centroid).SquaredNorm();
    spread = std::sqrt(spread / S(points.size()));
    if (spread == 0)
      throw std::invalid_argument("ShotResection: the points are all in one place");

    S rms = RmsError(shot, points, pixels);
    for (int round = 0; round < maxRounds; ++round) {
      Objective objective(shot, points, pixels, centroid, spread, fitFocal);
      S x[7] = {0, 0, 0, 0, 0, 0, 0};
      // The trust region starts as large as the error, in pixels, and ends well below the
      // tolerance asked of the rounds
      min_newuoa<S, Objective>(objective.Size(), x, objective, std::max(rms, S(1)), tolerance * S(1e-3), 5000);
      const ShotType fitted = objective.ToShot(x);
      const S fittedRms = RmsError(fitted, points, pixels);
      if (!(fittedRms < rms))
        break;
      shot = fitted;
      const S gain = rms - fittedRms;
      rms = fittedRms;
      if (gain < tolerance)
        break;
    }
    return rms;
  }

private:
  // The parameters about one shot: a rotation vector and a translation in the camera's frame,
  // about the centroid as the camera sees it, and a change of FocalMm, each scaled to move
  // the projected points by about one pixel per unit.
  class Objective
  {
  public:
    Objective(const ShotType &reference, const std::vector<Point3<S>> &points,
              const std::vector<Point2<S>> &pixels, const Point3<S> &centroid, S spread, bool fitFocal)
      : reference(reference), points(points), pixels(pixels), fitFocal(fitFocal)
    {
      // Not ConvertWorldToCameraCoordinates, which turns z to point forward
      centroidCam = reference.Extrinsics.Rot() * (centroid - reference.Extrinsics.Tra());
      // Steps small enough to measure a derivative, in the units of each parameter
      const S steps[7] = {S(1e-4), S(1e-4), S(1e-4), S(1e-4) * spread, S(1e-4) * spread, S(1e-4) * spread,
                          S(1e-4) * reference.Intrinsics.FocalMm};
      for (int i = 0; i < Size(); ++i) {
        S raw[7] = {0, 0, 0, 0, 0, 0, 0};
        raw[i] = steps[i];
        const S moved = Displacement(Apply(raw)) / steps[i];
        // A parameter that moves nothing stays where it is
        scale[i] = moved > 0 ? S(1) / moved : S(0);
      }
    }

    int Size() const { return fitFocal ? 7 : 6; }

    S operator()(int n, S *x) const
    {
      const ShotType shot = ToShot(x, n);
      S sum = 0;
      for (size_t i = 0; i < points.size(); ++i)
        sum += (shot.Project(points[i]) - pixels[i]).SquaredNorm();
      return sum / S(points.size());
    }

    ShotType ToShot(const S *x, int n = -1) const
    {
      S raw[7] = {0, 0, 0, 0, 0, 0, 0};
      for (int i = 0; i < (n < 0 ? Size() : n); ++i)
        raw[i] = x[i] * scale[i];
      return Apply(raw);
    }

  private:
    const ShotType &reference;
    const std::vector<Point3<S>> &points;
    const std::vector<Point2<S>> &pixels;
    bool fitFocal;
    Point3<S> centroidCam;
    S scale[7] = {1, 1, 1, 1, 1, 1, 1};

    // The reference shot moved by raw parameters. A camera point q becomes
    // dR (q - c) + c + dt, with c the centroid in camera coordinates: so the new rotation is
    // dR R and the viewpoint moves by R^T (c - dR^T (c + dt)).
    ShotType Apply(const S *raw) const
    {
      ShotType shot = reference;
      const Point3<S> axis(raw[0], raw[1], raw[2]);
      const S angle = axis.Norm();
      Matrix44<S> dR;
      dR.SetIdentity();
      if (angle > 0)
        dR.SetRotateRad(angle, axis / angle);
      const Matrix44<S> R = reference.Extrinsics.Rot();
      Matrix44<S> dRt = dR, Rt = R;
      Transpose(dRt);
      Transpose(Rt);
      const Point3<S> dt(raw[3], raw[4], raw[5]);
      shot.Extrinsics.SetRot(dR * R);
      shot.Extrinsics.SetTra(reference.Extrinsics.Tra() + Rt * (centroidCam - dRt * (centroidCam + dt)));
      shot.Intrinsics.FocalMm = reference.Intrinsics.FocalMm + raw[6];
      return shot;
    }

    // How far, root mean square in pixels, the points move from where the reference shows them
    S Displacement(const ShotType &shot) const
    {
      S sum = 0;
      for (const Point3<S> &p : points)
        sum += (shot.Project(p) - reference.Project(p)).SquaredNorm();
      return std::sqrt(sum / S(points.size()));
    }
  };
};

typedef ShotResection<float>  ShotResectionf;
typedef ShotResection<double> ShotResectiond;

} // end namespace vcg

#endif // __VCGLIB_SHOT_RESECTION
