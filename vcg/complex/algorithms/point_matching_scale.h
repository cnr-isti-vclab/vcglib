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

#ifndef VCG_POINT_MATCHING_SCALE
#define VCG_POINT_MATCHING_SCALE

#include <vcg/math/matrix44.h>
#include <vcg/space/point3.h>
#include <vcg/space/box3.h>
#include <wrap/newuoa/include/newuoa.h>

namespace vcg {

template <class Scalar>
struct RotoTranslation
{
	RotoTranslation(){}
	Scalar _v[6];
	void toMatrix(vcg::Matrix44<Scalar> & m)
	{
		vcg::Matrix44<Scalar> rot,tra;
		rot.FromEulerAngles(_v[0],_v[1],_v[2]);
		tra.SetTranslate(vcg::Point3<Scalar>(_v[3],_v[4],_v[5]));
		m = tra * rot;
	}
};

class PointMatchingScale {
private:
	static std::vector<vcg::Point3d> *fix;
	static std::vector<vcg::Point3d> *mov;
	static vcg::Box3d b;

public:
	/**
	 * Compute a scaling transformation that bring PMov point as close as possible to Pfix
	 */
	static void computeScalingMatchMatrix(
			vcg::Matrix44d &res,
			std::vector<vcg::Point3d> &Pfix,
			std::vector<vcg::Point3d> &Pmov)
	{
		fix = &Pfix;
		mov = &Pmov;
		b.SetNull();
		for(std::vector<vcg::Point3d>::iterator i = Pmov.begin(); i != Pmov.end(); ++i)
			b.Add(*i);

		double scale = 1.0;
		min_newuoa(1,&scale,errorScale);

		res.SetTranslate( b.Center()*(1.0-scale));
		res[0][0] = res[1][1] = res[2][2] = scale;
	}

	/**
	 * Compute a rototranslation + scaling transformation that bring PMov point as close as possible to Pfix
	 */
	static void computeRotoTranslationScalingMatchMatrix(
			vcg::Matrix44d &res,
			std::vector<vcg::Point3d> &Pfix,
			std::vector<vcg::Point3d> &Pmov)
	{
		fix = &Pfix;
		mov = &Pmov;
		b.SetNull();
		for(std::vector<vcg::Point3d>::iterator i = Pmov.begin(); i != Pmov.end(); ++i)
			b.Add(*i);

		double x[7]={1.0,0.0,0.0,0.0,0.0,0.0,0.0};
		min_newuoa(7,&x[0],errorRotoTranslationScale);

		// rtm = rototranslation
		RotoTranslation<double> rt;
		vcg::Matrix44d rtm;
		for (unsigned int i = 0; i < 6; ++i)
			rt._v[i] = x[i+1];
		rt.toMatrix(rtm);

		// res= scaling w.r.t. barycenter
		res.SetTranslate( b.Center()*(1.0-x[0]));
		res[0][0] = res[1][1] = res[2][2] = x[0];
		res = rtm*res;
	}

	static double errorScale(int n, double *x)
	{
		assert(n==1); (void)n;
		double dist = 0;
		std::vector<vcg::Point3d>::iterator i = mov->begin();
		std::vector<vcg::Point3d>::iterator ifix = fix->begin();
		for(; i !=  mov->end(); ++i,++ifix)
			dist += vcg::SquaredDistance(((*i)-b.Center())*(*x)+b.Center() , *ifix);

		return dist;
	}

	static double errorRotoTranslationScale(int n, double *x) {
		assert(n==7); (void)n;
		double dist = 0;
		std::vector<vcg::Point3d>::iterator i = mov->begin();
		std::vector<vcg::Point3d>::iterator ifix = fix->begin();

		RotoTranslation<double> rt;
		vcg::Matrix44d m;
		for (unsigned int i = 0; i < 6; ++i)
			rt._v[i] = x[i+1];
		rt.toMatrix(m);

		for(; i !=  mov->end(); ++i,++ifix) {
			dist += vcg::SquaredDistance(	 m*(((*i)-b.Center())*(x[0])+b.Center()),*ifix);
		}
		return dist;
	}

};

}

#endif
