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

****************************************************************************/

#ifndef __VCGLIB_FACEFUNCTORS
#define __VCGLIB_FACEFUNCTORS

// vcg headers
#include <vcg/space/point3.h>
#include <vcg/space/box3.h>
#include <vcg/space/line3.h>
#include <vcg/space/intersection3.h>
#include <vcg/simplex/face/distance.h>

namespace vcg {

class FaceBarycenterFunctor {
public:
	template <class FACETYPE, class COORDTYPE>
	inline void operator () (const FACETYPE & f, COORDTYPE & b) {
		b = COORDTYPE::Construct(f.Barycenter());
	}
};

class FaceBoxFunctor {
public:
	template <class FACETYPE, class BOXTYPE>
	inline void operator () (const FACETYPE & f, BOXTYPE & b) {
		b.Set(Point3<BOXTYPE::ScalarType>::Construct(f.P(0)));
		b.Add(Point3<BOXTYPE::ScalarType>::Construct(f.P(1)));
		b.Add(Point3<BOXTYPE::ScalarType>::Construct(f.P(2)));
	}
};

class FacePointDistanceFunctor {
public:
	template <class FACETYPE, class COORDTYPE>
	inline bool operator () (const FACETYPE & f, const COORDTYPE & p, typename COORDTYPE::ScalarType & minDist, COORDTYPE & q) {
		const Point3<typename FACETYPE::ScalarType> fp = Point3<typename FACETYPE::ScalarType>::Construct(p);
		Point3<typename FACETYPE::ScalarType> fq;
		typename FACETYPE::ScalarType md = (typename FACETYPE::ScalarType)(minDist);
		const bool ret = face::PointDistance(f, fp, md, fq);
		minDist = (typename COORDTYPE::ScalarType)(md);
		q = COORDTYPE::Construct(fq);
		return (ret);
	}
};

template <bool BACKFACETEST = true>
class FaceRayIntersectFunctor {
public:
	template <class FACETYPE, class COORDTYPE>
	inline bool operator () (const FACETYPE & f, const COORDTYPE & rayOrigin, const COORDTYPE & rayDirection, typename COORDTYPE::ScalarType & t, COORDTYPE & q) {
		typedef typename COORDTYPE::ScalarType ScalarType;
		Line3<ScalarType, false> ln(rayOrigin, rayDirection);
		ScalarType a;
		ScalarType b;

		bool bret = Intersection(ln, f.P(0), f.P(1), f.P(2), a, b, t);
		if (BACKFACETEST) {
			if (!bret) {
				bret = Intersection(ln, f.P(0), f.P(2), f.P(1), a, b, t);
			}
		}

		if (bret) {
			q = rayOrigin + rayDirection * t;
			return (true);
		}
		return (false);
	}
};

} // end namespace vcg

#endif // #ifndef __VCGLIB_FACEFUNCTORS
