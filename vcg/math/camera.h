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

/****************************************************************************
  History
$Log: not supported by cvs2svn $

Revision 1.28  2008/09/09 11:13:27  dellepiane
new functions to handle distortion: should not affect previous stuff. tested but still error
prone...

Revision 1.28  2006/12/21 00:13:27  cignoni
Corrected a syntax error detected only by gcc.
Corrected the order of initialization in the constructor to match the declaration order

Revision 1.27  2006/12/18 16:02:55  matteodelle
minor eroor correction on variable names

Revision 1.26  2006/12/18 09:46:38  callieri
camera+shot revamp: changed field names to something with more sense, cleaning of various functions,
correction of minor bugs/incongruences, removal of the infamous reference in shot.

Revision 1.25  2005/12/12 16:52:55  callieri
Added Unproject, from 2D local space + Zdepth to 3D camera space. Added ViewportToLocal, inverse of
LocalToViewport

Revision 1.24  2005/12/01 01:03:37  cignoni
Removed excess ';' from end of template functions, for gcc compiling

Revision 1.23  2005/10/12 16:43:32  ponchio
Added IsOrtho...

Revision 1.22  2005/07/11 13:12:34  cignoni
small gcc-related compiling issues (typenames,ending cr, initialization order)

Revision 1.21  2005/07/01 10:55:42  cignoni
Removed default values from the implementation of SetCavalieri and SetIsometric

Revision 1.20  2005/06/29 14:59:03  spinelli
aggiunto:
- l' enum dei tipi  PERSPECTIVE,  ORTHO, ISOMETRIC,  CAVALIERI
- inline void SetCavalieri(...)
- inline void SetIsometric(...)

- modificato
- void SetOrtho( .. )

Revision 1.19  2005/02/22 10:57:58  tommyfranken
Corrected declaration and some syntax errors in GetFrustum

Revision 1.18  2005/02/21 18:11:07  ganovelli
GetFrustum moved from gl/camera to math/camera.h

Revision 1.17  2005/02/15 14:55:52  tommyfranken
added principal point

Revision 1.16  2005/01/18 16:40:50  ricciodimare
*** empty log message ***

Revision 1.15  2005/01/18 15:14:22  ponchio
Far and end are reserved.

Revision 1.14  2005/01/14 15:28:33  ponchio
vcg/Point.h -> vcg/point.h   (again!)

Revision 1.13  2005/01/05 13:25:29  ganovelli
aggiunte conversione di coordinate

Revision 1.12  2004/12/16 11:22:30  ricciodimare
*** empty log message ***

Revision 1.11  2004/12/16 11:21:03  ricciodimare
*** empty log message ***

Revision 1.10  2004/12/15 18:45:50  tommyfranken
*** empty log message ***

<<<<<<< camera.h
=======
Revision 1.8  2004/11/23 10:15:38  cignoni
removed comment in comment gcc warning

Revision 1.7  2004/11/03 09:40:53  ganovelli
Point?.h to point?.h

Revision 1.6  2004/11/03 09:32:50  ganovelli
SetPerspective and SetFrustum added (same parameters as in opengl)

>>>>>>> 1.8
Revision 1.4  2004/10/07 14:39:57  fasano
Remove glew.h include

Revision 1.3  2004/10/07 14:22:38  ganovelli
y axis reverse in projecting (!)

Revision 1.2  2004/10/05 19:04:25  ganovelli
version 5-10-2004 in progress

Revision 1.1  2004/09/15 22:58:05  ganovelli
re-creation

Revision 1.2  2004/09/06 21:41:30  ganovelli
*** empty log message ***

Revision 1.1  2004/09/03 13:01:51  ganovelli
creation

****************************************************************************/

#ifndef __VCGLIB_CAMERA
#define __VCGLIB_CAMERA

#include <vcg/math/similarity.h>
#include <vcg/space/point2.h>
#include <vcg/space/point3.h>

namespace vcg {

template<class S>
class Camera
{
public:
	typedef S ScalarType;

	enum CameraType { PERSPECTIVE = 0, ORTHO = 1, ISOMETRIC = 2, CAVALIERI = 3 };

	Camera();

	template<class Q>
	static Camera Construct(const Camera<Q>& t);

	void SetOrtho(S l, S r, S b, S t, vcg::Point2<int> viewport);
	bool IsOrtho() const;

	//--- Set-up methods

	void SetPerspective(S AngleDeg, S AspectRatio, S Focal, vcg::Point2<int> Viewport);
	void SetCavalieri(S sx, S dx, S bt, S tp, S Focal, vcg::Point2<int> Viewport);
	void SetIsometric(S sx, S dx, S bt, S tp, S Focal, vcg::Point2<int> Viewport);
	void SetFrustum(S dx, S sx, S bt, S tp, S Focal, vcg::Point2<int> Viewport);

	vcg::Matrix44<S> GetMatrix(S nearVal, S farVal);
	void             GetFrustum(S& sx, S& dx, S& bt, S& tp, S& nr) const;

	//--- Space transformation methods

	vcg::Point2<S> Project(const vcg::Point3<S>& p) const;
	vcg::Point3<S> UnProject(const vcg::Point2<S>& p, const S& d) const;
	vcg::Point2<S> LocalToViewportPx(const vcg::Point2<S>& p) const;
	vcg::Point2<S> ViewportPxToLocal(const vcg::Point2<S>& p) const;
	vcg::Point2<S> ViewportPxTo_neg1_1(const vcg::Point2<S>& p) const;
	vcg::Point2<S> Neg1_1ToViewportPx(const vcg::Point2<S>& p) const;
	vcg::Point2<S> LocalTo_0_1(const vcg::Point2<S>& p) const;
	vcg::Point2<S> LocalTo_neg1_1(const vcg::Point2<S>& p) const;

	vcg::Point2<S> UndistortedToDistorted(vcg::Point2<S> u) const;
	vcg::Point2<S> DistortedToUndistorted(vcg::Point2<S> d) const;

	bool operator==(const Camera<S>& oth) const;
	bool operator!=(const Camera<S>& oth) const;

	//------ camera intrinsics
	ScalarType FocalMm;      /// Focal Distance: the distance between focal center and image plane.
							 /// Expressed in mm
	Point2<int> ViewportPx;  /// Dimension of the Image Plane (in pixels)
	Point2<S>   PixelSizeMm; /// Dimension in mm of a single pixel
	Point2<S>   CenterPx;    /// Position of the projection of the focal center on the image plane.
							 /// Expressed in pixels
	Point2<S>
		DistorCenterPx; /// Position of the radial distortion center on the image plane in pixels
	std::array<S, 4> k; /// 1st & 2nd order radial lens distortion coefficient (only the first 2 terms are used)

	CameraType cameraType; /// Type of camera: PERSPECTIVE,ORTHO,ISOMETRIC,CAVALIERI
};

} // namespace vcg

#include "camera.ipp"

#endif // __VCGLIB_CAMERA
