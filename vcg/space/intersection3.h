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



#ifndef __VCGLIB_INTERSECTION_3
#define __VCGLIB_INTERSECTION_3

#include <vcg/space/point3.h>

namespace vcg {
// sphere line

template<class T>
inline bool Intersection( const Sphere3<T> & sp, const Line3<T> & li, Point3<T> & p0,Point3<T> & p1 ){

	// Per prima cosa si sposta il sistema di riferimento 
	// fino a portare il centro della sfera nell'origine
	Point3<T> neworig=li.Origin()-sp.Center();
	// poi si risolve il sistema di secondo grado (con maple...)
	 T t1 = li.Direction().x()*li.Direction().x();
	 T t2 = li.Direction().y()*li.Direction().y();
   T t3 = li.Direction().z()*li.Direction().z();
   T t6 = neworig.y()*li.Direction().y();
   T t7 = neworig.x()*li.Direction().x();
   T t8 = neworig.z()*li.Direction().z();
   T t15 = sp.Radius()*sp.Radius();
   T t17 = neworig.z()*neworig.z();
   T t19 = neworig.y()*neworig.y();
   T t21 = neworig.x()*neworig.x();
   T t28 = 2.0*t7*t6+2.0*t6*t8+2.0*t7*t8+t1*t15-t1*t17-t1*t19-t2*t21+t2*t15-t2*t17-t3*t21+t3*t15-t3*t19;
	if(t28<0) return false;
	T t29 = sqrt(t28);      
	T val0 = 1/(t1+t2+t3)*(-t6-t7-t8+t29); 
	T val1 = 1/(t1+t2+t3)*(-t6-t7-t8-t29);

	p0=li.P(val0);
	p1=li.P(val1);
return true;
}


} // end namespace
#endif