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
Revision 1.5  2004/09/28 10:23:28  ponchio
Various generic changes.

Revision 1.4  2004/05/12 20:55:18  ponchio
*** empty log message ***

Revision 1.3  2004/03/31 15:06:41  ponchio
#include <camera> -> #include <view>

Revision 1.2  2004/03/25 14:55:25  ponchio
Adding copyright.


****************************************************************************/

#ifndef FRUSTUM_H
#define FRUSTUM_H

#include <wrap/gui/view.h>
#include <vcg/space/plane3.h>
#include <vcg/space/line3.h>

#include <iostream>
using namespace std;

namespace vcg {

template <class T> class Frustum: public View<T> {
public:                                            
  void GetView();  
  Point3<T> ViewPoint();
  T Resolution(float dist = 1);
  bool IsOutside(Point3<T> &point);
  bool IsOutside(Point3<T> &point, T radius);
  T Distance(Point3<T> &point, int plane);  
  
protected:
  T resolution;  
  Plane3<T> planes[6];  
  Point3<T> view_point;  
};


typedef Frustum<float> Frustumf;
typedef Frustum<double> Frustumd;


//Implementation
template <class T> Point3<T> Frustum<T>::ViewPoint() {
  return view_point;  
}

template <class T> T Frustum<T>::Resolution(float dist) {
  return resolution * dist;  
}

template <class T> bool Frustum<T>::IsOutside(Point3<T> &point) {
  Point3<T> r = Project(point);
  if(r[0] < viewport[0] || r[0] > viewport[0]+viewport[2] ||
     r[1] < viewport[1] || r[1] > viewport[1]+viewport[3]) 
    return true;
  return false;
}

template <class T> bool Frustum<T>::IsOutside(Point3<T> &point, T radius) {
  for(int i = 0; i < 4; i++) {
    T dist = Distance(point, i);   
    if(dist < -radius) 
      return true;  
  }
  return false;
}

template <class T> T Frustum<T>::Distance(Point3<T> &point, int plane) {    
  return vcg::Distance(planes[plane], point);  
}

template <class T> void Frustum<T>::GetView() {
  View<T>::GetView();
  
  int t = viewport[1] + viewport[3];
  int b = viewport[1];
  int r = viewport[0] + viewport[2];
  int l = viewport[0];
  
  Point3<T> nw = UnProject(Point3<T>(l, b, 0));
  Point3<T> sw = UnProject(Point3<T>(l, t, 0));
  Point3<T> ne = UnProject(Point3<T>(r, b, 0));
  Point3<T> se = UnProject(Point3<T>(r, t, 0));
  Point3<T> NW = UnProject(Point3<T>(l, b, 1));
  Point3<T> SW = UnProject(Point3<T>(l, t, 1));
  Point3<T> NE = UnProject(Point3<T>(r, b, 1));
  Point3<T> SE = UnProject(Point3<T>(r, t, 1));

  view_point = View<T>::ViewPoint();  	

  planes[0].Init(nw, NW, NE);  
  planes[1].Init(ne, NE, SE);
  planes[2].Init(se, SE, SW);
  planes[3].Init(sw, SW, NW);
  planes[4].Init(se, sw, nw);
  planes[5].Init(SW, SE, NE);   

  //compute resolution: sizeo of a pixel unitary distance from view_point
  resolution = ((ne + NE) - (nw + NW)).Norm() /
               (viewport[2] * ((ne + NE) - view_point).Norm());
}

}//namespace

#endif

  
	

  
