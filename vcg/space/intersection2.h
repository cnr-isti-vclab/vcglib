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


****************************************************************************/



#ifndef __VCGLIB_INTERSECTION_2
#define __VCGLIB_INTERSECTION_2

#include <vcg/space/point2.h>
#include <vcg/space/triangle2.h>



namespace vcg {
/** \addtogroup space */
/*@{*/
/** 
    Function computing the intersection between couple of geometric primitives in
    2 dimension
*/

/// return true if the algle is convex (right rotation)
template<class SCALAR_TYPE>
    inline bool Convex(const Point2<SCALAR_TYPE> & p0,const Point2<SCALAR_TYPE> & p1,const Point2<SCALAR_TYPE> & p2)
{
  return (((p0-p1)^(p2-p1))<=0);
}

/// interseciton between point and triangle
template<class SCALAR_TYPE>
    inline bool Intersection( const Triangle2<SCALAR_TYPE> & t,const Point2<SCALAR_TYPE> & p)
{
  Point2<SCALAR_TYPE> p0=t.P0(0);
  Point2<SCALAR_TYPE> p1=t.P0(1);
  Point2<SCALAR_TYPE> p2=t.P0(2);
  if (!Convex(p0,p1,p2))
    std::swap<Point2<SCALAR_TYPE> >(p1,p2);
  return((Convex(p,p0,p1))&&(Convex(p,p1,p2))&&(Convex(p,p2,p0)));
  //return((Convex(p,p0,p1))&&(Convex(p,p1,p2))&&(Convex(p,p2,p0)));
}
/*@}*/
} // end namespace
#endif
