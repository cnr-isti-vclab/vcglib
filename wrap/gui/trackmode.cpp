/****************************************************************************
* VCGLib                                                            o o     *
* Visual and Computer Graphics Library                            o     o   *
*                                                                _   O  _   *
* Copyright(C) 2004                                                \/)\/    *
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
Revision 1.2  2004/03/25 14:55:25  ponchio
Adding copyright.


****************************************************************************/

#include "trackmode.h"
#include <vcg/space/point3.h>
#include <vcg/math/similarity.h>

using namespace vcg;

Similarityf SphereMode::Apply(const Point3f &p, const Similarityf & /* m */) {
  float u = p[0];
	float w = p[1];
  float thr = 1/math::Sqrt(2.0f);       //in the plane x-y distance from origin, above this use hyperboloid
	
  float dist = math::Sqrt(u * u + w * w);
	Point3f result;
	if(dist < thr) 	{                    // First case: The ray is nearer to the sphere than r/sqrt(2) 
    float z = math::Sqrt(1 - u * u - w* w);
		result = Point3f(u, w, z);
	} else {                             // Second case: The ray should hit the 1/d hyperboloid
		float a = thr;
		result = Point3f(u, w, -a*(u*u + w * w) + 3*a/2);
		result.Normalize();
	}
	if(result == Point3f(0, 0, 1))
		return Similarityf().SetIdentity();

	Point3f axis = Point3f(0, 0, 1)^result; /* Axis of rotation */
	axis.Normalize();
	Point3f d = result - Point3f(0, 0, 1);
	float t = d.Norm() / 2.0f;
	if(t > thr)
		t += (t - thr) * 0.7f;
	if (t > 1.0f) t = 1.0f;
	if (t < -1.0f) t = -1.0f;
  float phi = 2 * math::Asin(t);
	//return Similarityf().SetRotate(phi * 180/(float)M_PI, axis);
  return Similarityf().SetRotate(phi, axis);
}

Similarityf PlaneMode::Apply(const Point3f &p, const Similarityf &a) {
  return Similarityf(Point3f(p[0], p[1], 0));
  Point3f r = x * a;
	Point3f u = y * a;
	int leading = 0; //leadiing x.
	if(fabs(u[2]) < fabs(r[2])) //sceglie l'asse principale: quello che piu' e' parallelo al piano di vista.
		leading = 1;
	r[2] = 0;
	u[2] = 0;
	if(r == Point3f(0,0,0)) //casi degeneri: un asse e' perpendicolare al piano di vista.
		r = Point3f(0, 1, 0);
	if(u == Point3f(0,0,0))
		u = Point3f(0, 1, 0);
	r.Normalize();
	u.Normalize();
	float cu, cr;
	if(leading == 0) { //leading x
		if(u == r || u == -r) { //caso degenere: i due assi si proiettano sullo stesso
			u[0] = -r[1];
			u[1] = r[0];
		}
		u = u - r * (r * u);
		u.Normalize();
	} else {
		if(r == u || r == -u) { //caso degenere: i due assi si proiettano sullo stesso
			r[0] = -u[1];
			r[1] = u[0];
		}
		r = r - u * (u * r);
		r.Normalize();
	}
	cr = r * p;
	cu = u * p;
	return Similarityf(x * cr + y * cu);
}