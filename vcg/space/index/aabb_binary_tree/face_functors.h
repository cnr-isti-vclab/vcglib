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
Revision 1.1  2005/09/26 18:33:16  m_di_benedetto
First Commit.


****************************************************************************/

#ifndef __VCGLIB_FACEFUNCTORS
#define __VCGLIB_FACEFUNCTORS

// vcg headers
#include <vcg/space/intersection3.h>

namespace vcg {

template <bool BACKFACETEST = true>
class FaceRayIntersectFunctor {
public:
	template <class FACETYPE, class SCALARTYPE>
	inline bool operator () (const FACETYPE & f, const Ray3<SCALARTYPE> & ray, SCALARTYPE & t) {
		typedef typename SCALARTYPE ScalarType;
		ScalarType a;
		ScalarType b;

		bool bret = Intersection(ray, f.P(0), f.P(1), f.P(2), a, b, t);
		if (BACKFACETEST) {
			if (!bret) {
				bret = Intersection(ray, f.P(0), f.P(2), f.P(1), a, b, t);
			}
		}
		return (bret);
	}
};

} // end namespace vcg

#endif // #ifndef __VCGLIB_FACEFUNCTORS
