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
Revision 1.2  2004/03/25 14:55:25  ponchio
Adding copyright.


****************************************************************************/

#ifndef FRUSTUM_H
#define FRUSTUM_H

#include <wrap/gui/view.h>
#include <vcg/space/plane3.h>
#include <vcg/space/line3.h>


namespace vcg {

  template <class T> class Frustum: public Camera {
public:                                            
  void GetView();  
  bool IsOutside(Point3<T> &point);
  bool IsOutside(Point3<T> &point, T radius);
  T Distance(Point3<T> &point, int plane);  
  Point3<T> ViewPoint();

protected:
  T resolution;  
  Plane3<T> planes[6];  
  Point3<T> view_point;  
};


//Implementation
template <class T> Point3<T> Frustum<T>::ViewPoint() {
  return view_point;  
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
  return Distance<T>(point, planes[plane]);  
}

template <class T> void Frustum<T>::GetView() {
  Camera::GetView();
  
  Point3d NE, SE, SW, NW, ne, se, sw, nw;
  int t = viewport[1] + viewport[3];
  int b = viewport[1];
  int r = viewport[0] + viewport[2];
  int l = viewport[0];
	
  Point3d NE, SE, SW, NW, ne, se, sw, nw;
	gluUnProject(l, b, 0, model_matrix, proj_matrix, viewport, &nw[0], &nw[1], &nw[2]);
	gluUnProject(l, t, 0, model_matrix, proj_matrix, viewport, &sw[0], &sw[1], &sw[2]);
	gluUnProject(r, b, 0, model_matrix, proj_matrix, viewport, &ne[0], &ne[1], &ne[2]);
	gluUnProject(r, t, 0, model_matrix, proj_matrix, viewport, &se[0], &se[1], &se[2]);
	gluUnProject(l, b, 1, model_matrix, proj_matrix, viewport, &NW[0], &NW[1], &NW[2]);
	gluUnProject(l, t, 1, model_matrix, proj_matrix, viewport, &SW[0], &SW[1], &SW[2]);
	gluUnProject(r, b, 1, model_matrix, proj_matrix, viewport, &NE[0], &NE[1], &NE[2]);
	gluUnProject(r, t, 1, model_matrix, proj_matrix, viewport, &SE[0], &SE[1], &SE[2]);

  view_point = Camera::ViewPoint();  	
  
  planes[0].init(view_point,             Point3<T>().Import(nw), Point3<T>().Import(ne));  
  planes[1].init(view_point,             Point3<T>().Import(ne), Point3<T>().Import(se));
  planes[2].init(view_point,             Point3<T>().Import(se), Point3<T>().Import(sw));
  planes[3].init(view_point,             Point3<T>().Import(sw), Point3<T>().Import(nw));
  planes[4].init(Point3<T>().Import(se), Point3<T>().Import(sw), Point3<T>().Import(nw));
  planes[5].init(Point3<T>().Import(SW), Point3<T>().Import(SE), Point3<T>().Import(NE));   

  for(int i = 0; i < 6; i++)
    planes[i].Normalize();

  //calcoliamo la risoluzione: dimenzione di un pixel a distanza 1 dal view_point
  resolution =  (T)((ne + NE) - (nw + NW)).Norm() /( viewport[2] * ((ne + NE) - (nw + NW)).Norm());
}

}//namespace

#endif

  
	

  
