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
Revision 1.5  2004/05/05 08:21:55  cignoni
syntax errors in inersection plane line.

Revision 1.4  2004/05/04 02:37:58  ganovelli
Triangle3<T> replaced by TRIANGLE
Segment<T> replaced by EDGETYPE

Revision 1.3  2004/04/29 10:48:44  ganovelli
error in plane segment corrected

Revision 1.2  2004/04/26 12:34:50  ganovelli
plane line
plane segment
triangle triangle added

Revision 1.1  2004/04/21 14:22:27  cignoni
Initial Commit


****************************************************************************/



#ifndef __VCGLIB_INTERSECTION_3
#define __VCGLIB_INTERSECTION_3

#include <vcg/space/point3.h>
#include <vcg/space/line3.h>
#include <vcg/space/plane3.h>
#include <vcg/space/segment3.h>
#include <vcg/space/sphere3.h>
#include <vcg/space/triangle3.h>
#include <vcg/space/intersection/triangle_triangle3.h>



/** \addtogroup space */
/*@{*/
/** 
    Function computing the intersection between couple of geometric primitives in
    3 dimension
*/

namespace vcg {
  /// interseciton between sphere and line
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

  /// intersection between line and plane
  template<class T>
    inline bool Intersection( const Plane3<T> & pl, const Line3<T> & li, Point3<T> & po){
    const T epsilon = T(1e-8);

    T k = pl.Direction() * li.Direction();						// Compute 'k' factor
    if( (k > -epsilon) && (k < epsilon))
      return false;
    T r = (pl.Offset() - pl.Direction()*li.Origin())/k;	// Compute ray distance
    po = li.Origin() + li.Direction()*r;
    return true;
  }

  /// intersection between segment and plane
  template<typename SEGMENTTYPE>
    inline bool Intersection( const Plane3<typename SEGMENTTYPE::ScalarType> & pl, 
			      const SEGMENTTYPE & sg, 
			      Point3<typename SEGMENTTYPE::ScalarType> & po){
    typedef typename SEGMENTTYPE::ScalarType T;
    const T epsilon = T(1e-8);

    T k = pl.Direction() * (sg.P1()-sg.P0());
    if( (k > -epsilon) && (k < epsilon))
      return false;
    T r = (pl.Offset() - pl.Direction()*sg.P0())/k;	// Compute ray distance
    if( (r<0) || (r > 1.0))
      return false;
    po = sg.P0()*(1-r)+sg.P1() * r;
    return true;
  }

  /// intersection between plane and triangle 
  // not optimal: uses plane-segment intersection (and the fact the two or none edges can be intersected)
  template<typename TRIANGLETYPE> 
    inline bool Intersection( const Plane3<typename TRIANGLETYPE::ScalarType> & pl, 
			      const  TRIANGLETYPE & tr, 
			      Segment3<typename TRIANGLETYPE::ScalarType> & sg){
    typedef typename TRIANGLETYPE::ScalarType T;
    if(Intersection(pl,Segment3<T>(tr.P(0),tr.P(1)),sg.P0())){
      if(Intersection(pl,Segment3<T>(tr.P(0),tr.P(2)),sg.P1()))
	return true;
      else
	{
	  Intersection(pl,Segment3<T>(tr.P(1),tr.P(2)),sg.P1());
	  return true;
	}
    }else
      {
	if(Intersection(pl,Segment3<T>(tr.P(1),tr.P(2)),sg.P0()))
	  {
	    Intersection(pl,Segment3<T>(tr.P(0),tr.P(2)),sg.P1());
	    return true;
	  }
      }
    return false;
  }

  /// intersection between two triangles
  template<typename TRIANGLETYPE> 
    inline bool Intersection(const TRIANGLETYPE & t0,const TRIANGLETYPE & t1){
    return NoDivTriTriIsect(t0.P0(0),t0.P0(1),t0.P0(2),
			    t1.P0(0),t1.P0(1),t1.P0(2));
  }
  template<class T>
    inline bool Intersection(Point3<T> V0,Point3<T> V1,Point3<T> V2,
			     Point3<T> U0,Point3<T> U1,Point3<T> U2){
    return NoDivTriTriIsect(V0,V1,V2,U0,U1,U2);
  }

  template<class T>
    inline bool Intersection(Point3<T> V0,Point3<T> V1,Point3<T> V2,
			     Point3<T> U0,Point3<T> U1,Point3<T> U2,int *coplanar,
			     Point3<T> &isectpt1,Point3<T> &isectpt2){

    return tri_tri_intersect_with_isectline(V0,V1,V2,U0,U1,U2,
					    coplanar,isectpt1,isectpt2);
  }
  template<typename TRIANGLETYPE,typename SEGMENTTYPE >
    inline bool Intersection(const TRIANGLETYPE & t0,const TRIANGLETYPE & t1,bool &coplanar,
			     SEGMENTTYPE  & sg){
    Point3<SEGMENTTYPE::PointType> ip0,ip1; 
    return  tri_tri_intersect_with_isectline(t0.P0(0),t0.P0(1),t0.P0(2),
					     t1.P0(0),t1.P0(1),t1.P0(2),
					     coplanar,sg.P0(),sg.P1()
					     );              
  }



} // end namespace
#endif
