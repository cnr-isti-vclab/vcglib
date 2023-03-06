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

#include "shot.h"

namespace vcg {

template<class S, class RotationType>
Shot<S, RotationType>::Shot() : Intrinsics(), Extrinsics()
{
	Extrinsics.SetIdentity();
}

template<class S, class RotationType>
Shot<S, RotationType>::Shot(const Camera<S>& i, const ReferenceFrame& e) :
		Intrinsics(), Extrinsics()
{
	Intrinsics = i;
	Extrinsics = e;
}

template<class S, class RotationType>
Shot<S, RotationType>::Shot(const Camera<S>& c) : Intrinsics(), Extrinsics()
{
	Intrinsics = c;
	Extrinsics.SetIdentity();
}

template<class S, class RotationType>
template<class Q>
Shot<S, RotationType> Shot<S, RotationType>::Construct(const Shot<Q>& b)
{
	ReferenceFrame r;
	r.SetRot(Matrix44<S>::Construct(b.Extrinsics.Rot()));
	r.SetTra(Point3<S>::Construct(b.Extrinsics.Tra()));
	return Shot(Camera<S>::Construct(b.Intrinsics), r);
}

/**
 * @brief get the i-th axis of the coordinate system of the camera
 */
template<class S, class RotationType>
vcg::Point3<S> Shot<S, RotationType>::Axis(const int& i) const
{
	vcg::Matrix44<S> m;
	Extrinsics.rot.ToMatrix(m);
	vcg::Point3<S> aa = m.GetRow3(i);
	return aa;
}

/**
 * @brief Get the viewdir
 */
template<class S, class RotationType>
const vcg::Point3<S> Shot<S, RotationType>::GetViewDir() const
{
	return Extrinsics.Rot().GetRow3(2);
}

/**
 * @brief Get the viewpoint
 */
template<class S, class RotationType>
const vcg::Point3<S> Shot<S, RotationType>::GetViewPoint() const
{
	return Extrinsics.tra;
}

/**
 * @brief set the viewpoint
 */
template<class S, class RotationType>
void Shot<S, RotationType>::SetViewPoint(const vcg::Point3<S>& viewpoint)
{
	Extrinsics.SetTra(viewpoint);
}

/**
 * @brief get fov from focal
 */
template<class S, class RotationType>
float Shot<S, RotationType>::GetFovFromFocal() const
{
	double viewportYMm = Intrinsics.PixelSizeMm[1] * Intrinsics.ViewportPx[1];
	return 2 * (vcg::math::ToDeg(atanf(viewportYMm / (2 * Intrinsics.FocalMm))));
}

/**
 * @brief look at (point+up)
 */
template<class S, class RotationType>
void Shot<S, RotationType>::LookAt(const vcg::Point3<S>& z_dir, const vcg::Point3<S>& up)
{
	LookTowards(z_dir - GetViewPoint(), up);
}

/**
 * @brief look at (opengl-like)
 */
template<class S, class RotationType>
void Shot<S, RotationType>::LookAt(
	const S& eye_x,
	const S& eye_y,
	const S& eye_z,
	const S& at_x,
	const S& at_y,
	const S& at_z,
	const S& up_x,
	const S& up_y,
	const S& up_z)
{
	SetViewPoint(Point3<S>(eye_x, eye_y, eye_z));
	LookAt(Point3<S>(at_x, at_y, at_z), Point3<S>(up_x, up_y, up_z));
}

/**
 * @brief look towards (dir+up)
 */
template<class S, class RotationType>
void Shot<S, RotationType>::LookTowards(const vcg::Point3<S>& z_dir, const vcg::Point3<S>& up)
{
	vcg::Point3<S> x_dir = up ^ -z_dir;
	vcg::Point3<S> y_dir = -z_dir ^ x_dir;

	Matrix44<S> m;
	m.SetIdentity();
	*(vcg::Point3<S>*) &m[0][0] = x_dir / x_dir.Norm();
	*(vcg::Point3<S>*) &m[1][0] = y_dir / y_dir.Norm();
	*(vcg::Point3<S>*) &m[2][0] = -z_dir / z_dir.Norm();

	Extrinsics.rot.FromMatrix(m);
}

/**
 * @brief Sometimes the focal is given in pixels. In this case, this function can be used to convert
 * it in millimiters given the CCD width (in mm). This method should be moved in vcg::Camera().
 * Equivalent focal length is obtained by setting the ccd width to 35 mm.
 */
template<class S, class RotationType>
void Shot<S, RotationType>::ConvertFocalToMM(S ccdwidth)
{
	double ccd_width          = ccdwidth; // ccd is assumed conventionally to be 35mm
	double ccd_height         = (ccd_width * Intrinsics.ViewportPx[1]) / Intrinsics.ViewportPx[0];
	Intrinsics.PixelSizeMm[0] = (ccd_width / Intrinsics.ViewportPx[0]);
	Intrinsics.PixelSizeMm[1] = (ccd_height / Intrinsics.ViewportPx[1]);
	Intrinsics.FocalMm =
		(ccd_width * Intrinsics.FocalMm) / Intrinsics.ViewportPx[0]; // NOW FOCAL IS IN MM
}

/**
 * @brief Sometimes the 3D World coordinates are known up to a scale factor. This method adjust the
 * camera/shot parameters to account for the re-scaling of the World. If the intrisic parameters are
 * just reasonable values the cameras need only a re-positioning.
 */
template<class S, class RotationType>
void Shot<S, RotationType>::RescalingWorld(S scalefactor, bool adjustIntrinsics)
{
	// adjust INTRINSICS (if required)

	if (adjustIntrinsics) {
		Intrinsics.FocalMm = Intrinsics.FocalMm * scalefactor;
		double ccdwidth = static_cast<double>(Intrinsics.ViewportPx[0] * Intrinsics.PixelSizeMm[0]);
		double ccdheight =
			static_cast<double>(Intrinsics.ViewportPx[1] * Intrinsics.PixelSizeMm[1]);

		Intrinsics.PixelSizeMm[0] = (ccdwidth * scalefactor) / Intrinsics.ViewportPx[0];
		Intrinsics.PixelSizeMm[1] = (ccdheight * scalefactor) / Intrinsics.ViewportPx[1];
	}

	// adjust EXTRINSICS

	// rotation remains the same (!)
	// nothing to do..

	// the viewpoint should be modified according to the scale factor
	Extrinsics.tra *= scalefactor;
}

/**
 * @brief Given a pure roto-translation matrix (4-by-4) modify the reference frame accordingly.
 */
template<class S, class RotationType>
void Shot<S, RotationType>::ApplyRigidTransformation(const Matrix44<S>& M)
{
	Matrix44<S> rotM;
	Extrinsics.rot.ToMatrix(rotM);
	// roto-translate the viewpoint
	Extrinsics.tra     = M * Extrinsics.tra;
	Matrix44<S> newRot = rotM * M.transpose();
	newRot[3][0] = newRot[3][1] = newRot[3][2] = 0.0;

	Extrinsics.SetRot(newRot);
}

/**
 * @brief Given a similarity transformation modifies the reference frame accordingly.
 */
template<class S, class RotationType>
void Shot<S, RotationType>::ApplySimilarity(Matrix44<S> M)
{
	Matrix44<S> rotM;
	Extrinsics.rot.ToMatrix(rotM);

	// normalize
	M       = M * (1 / M.ElementAt(3, 3));
	M[3][3] = 1; // just for numeric precision

	// compute scale factor
	ScalarType scalefactor = 1.0 / pow(ScalarType(M.Determinant()), 1 / ScalarType(3.0));

	// roto-translate the viewpoint
	Extrinsics.tra = M * Extrinsics.tra;

	vcg::Matrix44<S> M2 = M;

	M2       = M2 * scalefactor; // remove the scaling
	M2[3][3] = 1.0;
	M2[0][3] = M2[1][3] = M2[2][3] = 0; // remove the translation

	rotM = rotM * M2.transpose();
	Extrinsics.SetRot(rotM);
}

/**
 * @brief Given a similarity transformation modifies the reference frame accordingly.
 */
template<class S, class RotationType>
void Shot<S, RotationType>::ApplySimilarity(const Similarity<S>& Sm)
{
	Matrix44<S> rotM;
	Extrinsics.rot.ToMatrix(rotM);

	// similarity decomposition
	vcg::Matrix44<S> R;
	Sm.rot.ToMatrix(R);
	vcg::Matrix44<S> T;
	T.SetIdentity();
	T.ElementAt(0, 3) = Sm.tra[0];
	T.ElementAt(1, 3) = Sm.tra[1];
	T.ElementAt(2, 3) = Sm.tra[2];
	vcg::Matrix44d S44;
	S44.SetIdentity();
	S44 *= Sm.sca;
	S44.ElementAt(3, 3) = 1.0;

	vcg::Matrix44<S> M = T * R * S44;

	// roto-translate the viewpoint
	Extrinsics.tra = M * Extrinsics.tra;

	vcg::Matrix44<S> M2 = M;

	M2 = M2 * (1.0 / Sm.sca);

	Extrinsics.rot = rotM * M2.transpose();

	Extrinsics.rot.ElementAt(3, 0) = 0;
	Extrinsics.rot.ElementAt(3, 1) = 0;
	Extrinsics.rot.ElementAt(3, 2) = 0;
	Extrinsics.rot.ElementAt(3, 3) = 1;
}

/**
 * @brief Convert a 3d point from world to camera coordinates (do not confuse with the Shot
 * reference frame)
 */
template<class S, class RotationType>
vcg::Point3<S> Shot<S, RotationType>::ConvertWorldToCameraCoordinates(const vcg::Point3<S>& p) const
{
	Matrix44<S> rotM;
	Extrinsics.rot.ToMatrix(rotM);
	vcg::Point3<S> cp = rotM * (p - GetViewPoint());
	cp[2]             = -cp[2];
	return cp;
}

/**
 * @brief Convert a 3d point from camera coordinates (do not confuse with the Shot reference frame)
 * to world coordinates
 */
template<class S, class RotationType>
vcg::Point3<S> Shot<S, RotationType>::ConvertCameraToWorldCoordinates(const vcg::Point3<S>& p) const
{
	Matrix44<S>    rotM;
	vcg::Point3<S> cp = p;
	cp[2]             = -cp[2];
	Extrinsics.rot.ToMatrix(rotM);
	cp = rotM.transpose() * cp + GetViewPoint();
	return cp;
}

/**
 * @brief Convert a 3d point from camera to world coordinates, uses inverse instead of trranspose
 * for non-exactly-rigid rotation matrices (such as calculated by tsai and garcia)
 */
template<class S, class RotationType>
vcg::Point3<S>
Shot<S, RotationType>::ConvertCameraToWorldCoordinates_Substitute(const vcg::Point3<S>& p) const
{
	Matrix44<S>    rotM;
	vcg::Point3<S> cp = p;
	cp[2]             = -cp[2];
	Extrinsics.rot.ToMatrix(rotM);
	cp = Inverse(rotM) * cp + GetViewPoint();
	return cp;
}

/**
 * @brief Project a 3d point from world coordinates to 2d camera viewport (the value returned is in
 * pixel)
 */
template<class S, class RotationType>
vcg::Point2<S> Shot<S, RotationType>::Project(const vcg::Point3<S>& p) const
{
	Point3<S> cp = ConvertWorldToCameraCoordinates(p);
	Point2<S> pp = Intrinsics.Project(cp);
	Point2<S> vp = Intrinsics.LocalToViewportPx(pp);
	return vp;
}

/**
 * @brief Inverse projection from 2d camera viewport (in pixels) to 3d world coordinates (it
 * requires the original depth of the point to unproject)
 */
template<class S, class RotationType>
vcg::Point3<S> Shot<S, RotationType>::UnProject(const vcg::Point2<S>& p, const S& d) const
{
	Point2<S> lp = Intrinsics.ViewportPxToLocal(p);
	Point3<S> cp = Intrinsics.UnProject(lp, d);
	Point3<S> wp = ConvertCameraToWorldCoordinates(cp);
	return wp;
}

/**
 * @brief Inverse projection from 2d camera viewport (in pixels) to 3d world coordinates (it
 * requires the original depth of the projected point) uses inverse instead of trranspose for
 * non-exactly-rigid rotation matrices (such as calculated by tsai and garcia)
 */
template<class S, class RotationType>
vcg::Point3<S>
Shot<S, RotationType>::UnProject_Substitute(const vcg::Point2<S>& p, const S& d) const
{
	Point2<S> lp = Intrinsics.ViewportPxToLocal(p);
	Point3<S> cp = Intrinsics.UnProject(lp, d);
	Point3<S> wp = ConvertCameraToWorldCoordinates_Substitute(cp);
	return wp;
}

/**
 * @brief Returns the distance of point p from camera plane (z depth), required for unprojection
 * operation
 */
template<class S, class RotationType>
S Shot<S, RotationType>::Depth(const vcg::Point3<S>& p) const
{
	return ConvertWorldToCameraCoordinates(p).Z();
}

/**
 * @brief Returns the (4-by-4) matrix M such that 3dpoint_in_world_coordinates = M *
 * 3dpoint_in_local_coordinates
 */
template<class S, class RotationType>
Matrix44<S> Shot<S, RotationType>::GetExtrinsicsToWorldMatrix() const
{
	Matrix44<S> rotM;
	Extrinsics.rot.ToMatrix(rotM);
	return Matrix44<S>().SetTranslate(Extrinsics.tra) * rotM.transpose();
}

/**
 * @brief Returns the (4-by-4) matrix M such that 3dpoint_in_local_coordinates = M *
 * 3dpoint_in_world_coordinates
 */
template<class S, class RotationType>
Matrix44<S> Shot<S, RotationType>::GetWorldToExtrinsicsMatrix() const
{
	Matrix44<S> rotM;
	Extrinsics.rot.ToMatrix(rotM);
	return rotM * Matrix44<S>().SetTranslate(-Extrinsics.tra);
}

/**
 * @brief multiply the current reference frame for the matrix passed
 * note: it is up to the caller to check the the matrix passed is a pure rototranslation
 */
template<class S, class RotationType>
void Shot<S, RotationType>::MultMatrix(vcg::Matrix44<S> m44)
{
	Extrinsics.tra = m44 * Extrinsics.tra;
	m44[0][3] = m44[1][3] = m44[2][3] = 0.0;                   // set no translation
	const S k                         = m44.GetRow3(0).Norm(); // compute scaling (assumed uniform)
	Extrinsics.rot                    = Extrinsics.rot * m44.transpose() * (1 / k);
}

/**
 * @brief Multiply the current reference frame for the similarity passed
 * note: it is up to the caller to check the the matrix passed is a pure rototranslation
 */
template<class S, class RotationType>
void Shot<S, RotationType>::MultSimilarity(const Similarity<S>& s)
{
	MultMatrix(s.Matrix());
}

template<class S, class RotationType>
bool Shot<S, RotationType>::IsValid() const
{
	return Intrinsics.PixelSizeMm[0] > 0 && Intrinsics.PixelSizeMm[1] > 0;
}

template<class S, class RotationType>
bool Shot<S, RotationType>::operator==(const Shot<S, RotationType> &oth) const
{
	return Intrinsics == oth.Intrinsics && Extrinsics == oth.Extrinsics;
}

template<class S, class RotationType>
bool Shot<S, RotationType>::operator!=(const Shot<S, RotationType> &oth) const
{
	return !(*this == oth);
}

template<class S, class RotationType>
Shot<S, RotationType>::ReferenceFrame::ReferenceFrame() : rot(), tra()
{
}

template<class S, class RotationType>
void Shot<S, RotationType>::ReferenceFrame::SetIdentity()
{
	rot.SetIdentity();
	tra = Point3<S>(0.0, 0.0, 0.0);
}

template<class S, class RotationType>
void vcg::Shot<S, RotationType>::ReferenceFrame::SetTra(const Point3<S> &tr)
{
	tra = tr;
}

template<class S, class RotationType>
void vcg::Shot<S, RotationType>::ReferenceFrame::SetRot(const RotationType &rt)
{
	rot = rt;
}

template<class S, class RotationType>
Point3<S> vcg::Shot<S, RotationType>::ReferenceFrame::Tra() const
{
	return tra;
}

template<class S, class RotationType>
RotationType vcg::Shot<S, RotationType>::ReferenceFrame::Rot() const
{
	return rot;
}

template<class S, class RotationType>
bool vcg::Shot<S, RotationType>::ReferenceFrame::operator==(
	const Shot<S, RotationType>::ReferenceFrame& oth) const
{
	return rot == oth.rot && tra == oth.tra;
}

template<class S, class RotationType>
bool vcg::Shot<S, RotationType>::ReferenceFrame::operator!=(
	const Shot<S, RotationType>::ReferenceFrame& oth) const
{
	return !(*this == oth);
}

} // namespace vcg
