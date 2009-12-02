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
// marco :
// comments
// corrected bug


/****************************************************************************
  History

****************************************************************************/

#ifndef __VCGLIB_VERTEX_DISTANCE
#define __VCGLIB_VERTEX_DISTANCE

#include <vcg/math/base.h>
#include <vcg/space/point3.h>


namespace vcg {
	namespace vertex{

template <class SCALARTYPE>
	class PointDistanceFunctor {
	public:
		typedef Point3<SCALARTYPE> QueryType;
		static inline const Point3<SCALARTYPE> & Pos(const QueryType & qt)  {return qt;}
		template <class VERTEXTYPE>

		/*
		 * @param v [IN] is a reference to the current object being tested,
		 * @param p [IN] is the query point,
		 * @param minDist [IN/OUT] is in input the reject distance and in output the closest distance,
	     * @param q [OUT] is the closest point.
		 * 
		 * @remarks The operator returns true if the closest distance is less than input reject distance.
		 *
		 */
        inline bool operator ()  (const VERTEXTYPE & v, const Point3<SCALARTYPE> & p, SCALARTYPE & minDist, Point3<SCALARTYPE> & q) const
		{
			// convert the coordinates of p from SCALARTYPE to VERTEXTYPE::ScalarType type
			const Point3<typename VERTEXTYPE::ScalarType> fp = Point3<typename VERTEXTYPE::ScalarType>::Construct(p);

			typename VERTEXTYPE::ScalarType md;		// distance between v and fp
			md = (v.P() - fp).Norm();

			if (md <= minDist) 
			{
				minDist = (SCALARTYPE)(md);		// minDist is updated to the closest distance
				q = v.P();						// q is the current closest point

				return true;
			}

			return false;
		}
	};

template <class VERTYPE>
class PointNormalDistanceFunctor {
	public:
		typedef VERTYPE QueryType; 
		typedef typename VERTYPE::ScalarType ScalarType;
		static inline const Point3<typename VERTYPE::ScalarType> &  Pos(const QueryType & qt)  {return qt.P();}

		static ScalarType & Alpha(){static ScalarType alpha = 1.0; return alpha;}
		static ScalarType & Beta(){static ScalarType beta= 1.0; return beta;}
		static ScalarType & Gamma(){static ScalarType gamma= 1.0; return gamma;}
		static ScalarType & InterPoint (){static ScalarType interpoint= 1.0; return interpoint;} 

		template <class VERTEXTYPE, class SCALARTYPE>
		inline bool operator () (const VERTEXTYPE & v, const VERTEXTYPE & vp, SCALARTYPE & minDist, Point3<SCALARTYPE> & q) {

			float  h = vcg::Distance(v.cP(),vp.P())  ;
			float  dev = InterPoint() * ( pow((ScalarType) (1-v.cN().dot(vp.cN())), (ScalarType)Beta()) / (Gamma()*h +0.1));
			if(h+dev < minDist){
				minDist = h+dev;
				q = v.P(); 
				return true;
			}

		//	minDist = h +0.0* (1-v.cN()*vp.cN()) / (h + 0.1);

			return false;
		}
	};

}	 // end namespace vertex
	
}	 // end namespace vcg


#endif

