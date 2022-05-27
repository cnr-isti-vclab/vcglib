/*****************************************************************************
 * VCGLib                                                            o o     *
 * Visual and Computer Graphics Library                            o     o   *
 *                                                                _   O  _   *
 * Copyright(C) 2004-2022                                           \/)\/    *
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

#ifndef __VCGLIB_SHOT
#define __VCGLIB_SHOT

#include <vcg/math/camera.h>
#include <vcg/math/similarity.h>
#include <vcg/space/point2.h>
#include <vcg/space/point3.h>

namespace vcg {

/**
 * @class Shot
 *
 * Shot is made of two elements:
 * - the Instrinsics paramaters, which are stored as a Camera type (see vcg/math/camera) and that
 *   determines how a point in the frame of the camera is projected in the 2D projection plane
 * - the Extrinsics parameters, which are stored in the class Shot (type ReferenceFrame)
 *   and that describe viewpoint and view direction.
 *
 * Some important notes about the usage of this class:
 * - The World coordinates system is assumed to be RIGHT-HANDED.
 * - The Shot reference frame is assumed to be RIGHT-HANDED.
 * - The associated Camera is assumed to point in the negative direction of the Z axis of the Shot
 *   coordinates system (reference frame). As a consequence, the Camera coordinates system is
 *   LEFT-HANDED.
 * - The Extrinsics parameters are kept as a rotation matrix "rot" and a translation vector "tra"
 *   The translation matrix "tra" corresponds to the viewpoint of the Shot while the rotation
 *   matrix "rot" corresponds to the axis of the reference frame by row, i.e.
 *    rot[0][0|1|2] == X axis
 *    rot[1][0|1|2] == Y axis
 *    rot[2][0|1|2] == Z axis
 *
 * It follows that the matrix made with the upper left 3x3 equal to rot and the 4th colum equal to
 * tra and (0,0,0,1) in the bottom row transform a point from world coordiantes to the reference
 * frame of the shot.
 */
template<class S, class RotationType = Matrix44<S>>
class Shot
{
public:
	typedef Camera<S> CameraType;
	typedef S         ScalarType;

	class ReferenceFrame;

	Camera<S>      Intrinsics; // the camera that made the shot
	ReferenceFrame Extrinsics; // the position and orientation of the camera

	Shot();

	Shot(const Camera<S>& i, const ReferenceFrame& e);

	Shot(const Camera<S>& c);

	template<class Q>
	static inline Shot Construct(const Shot<Q>& b);

	vcg::Point3<S> Axis(const int& i) const;

	const vcg::Point3<S> GetViewDir() const;
	const vcg::Point3<S> GetViewPoint() const;

	void SetViewPoint(const vcg::Point3<S>& viewpoint);

	float GetFovFromFocal() const;

	void LookAt(const vcg::Point3<S>& point, const vcg::Point3<S>& up);

	void LookAt(
		const S& eye_x,
		const S& eye_y,
		const S& eye_z,
		const S& at_x,
		const S& at_y,
		const S& at_z,
		const S& up_x,
		const S& up_y,
		const S& up_z);

	void LookTowards(const vcg::Point3<S>& z_dir, const vcg::Point3<S>& up);

	void ConvertFocalToMM(S ccdwidth);

	void RescalingWorld(S scalefactor, bool adjustIntrinsics);

	void ApplyRigidTransformation(const Matrix44<S>& M);

	void ApplySimilarity(Matrix44<S> M);

	void ApplySimilarity(const Similarity<S>& Sim);

	vcg::Point3<S> ConvertWorldToCameraCoordinates(const vcg::Point3<S>& p) const;

	vcg::Point3<S> ConvertCameraToWorldCoordinates(const vcg::Point3<S>& p) const;

	vcg::Point3<S> ConvertCameraToWorldCoordinates_Substitute(const vcg::Point3<S>& p) const;

	vcg::Point2<S> Project(const vcg::Point3<S>& p) const;

	vcg::Point3<S> UnProject(const vcg::Point2<S>& p, const S& d) const;

	vcg::Point3<S> UnProject_Substitute(const vcg::Point2<S>& p, const S& d) const;

	S Depth(const vcg::Point3<S>& p) const;

	Matrix44<S> GetExtrinsicsToWorldMatrix() const;

	Matrix44<S> GetWorldToExtrinsicsMatrix() const;

	void MultMatrix(vcg::Matrix44<S> m44);

	void MultSimilarity(const Similarity<S>& s);

	bool IsValid() const;

	bool operator==(const Shot<S, RotationType>& oth) const;
	bool operator!=(const Shot<S, RotationType>& oth) const;
};

template<class S, class RotationType>
class Shot<S, RotationType>::ReferenceFrame
{
	friend class Shot<S, RotationType>;

public:
	ReferenceFrame();

	void SetIdentity();
	void         SetTra(const Point3<S>& tr);
	void         SetRot(const RotationType& rt);
	Point3<S>    Tra() const;
	RotationType Rot() const;

	bool operator==(const Shot<S, RotationType>::ReferenceFrame& oth) const;
	bool operator!=(const Shot<S, RotationType>::ReferenceFrame& oth) const;

private:
	RotationType rot; // rotation
	Point3<S>    tra; // viewpoint
};

//--- utility definitions
typedef Shot<float>  Shotf;
typedef Shot<double> Shotd;

} // namespace vcg

#include "shot.ipp"

#endif // __VCGLIB_SHOT
