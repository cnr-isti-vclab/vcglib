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
  return Distance(point, planes[plane]);  
}

template <class T> void Frustum<T>::GetView() {
  View<T>::GetView();
  
  //  Point3d NE, SE, SW, NW, ne, se, sw, nw;
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

  /*  gluUnProject(l,b,0,model,proj,viewport,&nw[0], &nw[1], &nw[2]);
  gluUnProject(l,t,0,model,proj,viewport,&sw[0], &sw[1], &sw[2]);
  gluUnProject(r,b,0,model,proj,viewport,&ne[0], &ne[1], &ne[2]);
  gluUnProject(r,t,0,model,proj,viewport,&se[0], &se[1], &se[2]);
  gluUnProject(l,b,1,model,proj,viewport,&NW[0], &NW[1], &NW[2]);
  gluUnProject(l,t,1,model,proj,viewport,&SW[0], &SW[1], &SW[2]);
  gluUnProject(r,b,1,model,proj,viewport,&NE[0], &NE[1], &NE[2]);
  gluUnProject(r,t,1,model,proj,viewport,&SE[0], &SE[1], &SE[2]);*/
  
  view_point = View<T>::ViewPoint();  	
  
  planes[0].Init(view_point, nw, ne);  
  planes[1].Init(view_point, ne, se);
  planes[2].Init(view_point, se, sw);
  planes[3].Init(view_point, sw, nw);
  planes[4].Init(se, sw, nw);
  planes[5].Init(SW, SE, NE);   

  for(int i = 0; i < 6; i++)
    planes[i].Normalize();

  //calcoliamo la risoluzione: dimenzione di un pixel a distanza 1 dal view_point
  resolution = ((ne + NE) - (nw + NW)).Norm() /( viewport[2] * ((ne + NE) - view_point).Norm());
}

}//namespace

#endif

  
	

  
