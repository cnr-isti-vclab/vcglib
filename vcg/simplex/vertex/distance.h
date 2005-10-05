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

****************************************************************************/

#ifndef __VCGLIB_VERTEX_DISTANCE
#define __VCGLIB_VERTEX_DISTANCE

#include <vcg/math/base.h>
#include <vcg/space/point3.h>


namespace vcg {
	namespace vertex{

	class PointDistanceFunctor {
	public:
		template <class VERTEXTYPE, class SCALARTYPE>
		inline bool operator () (const VERTEXTYPE & v, const Point3<SCALARTYPE> & p, SCALARTYPE & minDist, Point3<SCALARTYPE> & q) {
			const Point3<typename VERTEXTYPE::ScalarType> fp = Point3<typename VERTEXTYPE::ScalarType>::Construct(p);
			typename VERTEXTYPE::ScalarType md;
			md=(v.P()-fp).Norm();
			bool ret = (md<=minDist);
			if (ret) minDist = (SCALARTYPE)(md);
			q = v.P();
			return (ret);
		}
	};


}	 // end namespace face
	
}	 // end namespace vcg


#endif

