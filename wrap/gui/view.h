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

#ifndef VCG_CAMERA_H
#define VCG_CAMERA_H

#include <vcg/space/point3.h>
#include <vcg/math/Matrix44.h>

#include <windows.h>
#include <GL/GL.h>

namespace vcg {

template <class T> class View {
public:
  void GetView();
  void SetView();                      
  Point3<T> Project(const Point3<T> &p) const;
  Point3<T> UnProject(const Point3<T> &p) const;
  Point3<T> ViewPoint();
  //convert coordinates range 0-1 to range 0 viewport[2]
  Point3<T> ScreenToViewport(const Point3<T> &p) const;
  //viceversa
  Point3<T> ViewportToScreen(const Point3<T> &p) const;

  Matrix44<T> proj;
  Matrix44<T> model;
  Matrix44<T> matrix;
  Matrix44<T> inverse;  
  int viewport[4];    
};

template <class T> void View<T>::GetView() {
  double m[16];
  glGetDoublev(GL_PROJECTION_MATRIX, m);
  for(int i = 0; i < 16; i++)
    proj[i] = (T)m[i];
	glGetDoublev(GL_MODELVIEW_MATRIX, m);
  for(int i = 0; i < 16; i++)
    model[i] = (T)m[i];
	
	glGetIntegerv(GL_VIEWPORT, viewport);
  
  matrix = model * proj;
  inverse = matrix;
  Invert(inverse);
}

template <class T> void View<T>::SetView() {
  
}

template <class T> Point3<T> View<T>::ViewPoint() {
  return inverse * Point3<T>(0, 0, 0);
  /*Matrix44d model(model_matrix);
	model.Invert();
  Point3d view = model * Point3d(0, 0, 0);	   
	return Point3<T>(view[0], view[1], view[2]); */
}

template <class T> Point3<T> View<T>::Project(const Point3<T> &p) const {
  Point3<T> r;
  r = p * matrix;  	
	r[0] = (r[0]+1)*(viewport[2]/(T)2.0)+viewport[0];
	r[1] =(r[1]+1)*(viewport[3]/(T)2.0)+viewport[1];
  r[1] = viewport[3]-r[1];
  return r;
  /*double r[3];
  gluProject(p[0], p[1], p[2], model_matrix, proj_matrix, viewport, &r[0], &r[1], &r[2]);
  return Point3<T>((T)r[0], (T)r[1], (T)r[2]);*/
}

template <class T> Point3<T> View<T>::UnProject(const Point3<T> &p) const {
  Point3<T> s = p;	
	s[0] = (p[0]- viewport[0])/ (viewport[2]/(T)2.0) - 1;
	s[1] = (p[1]- viewport[1])/ (viewport[3]/(T)2.0) - 1;
  s[1] = -s[1];
	s[2] = p[2];
  //s[1] = -s[1];   // pezza aggiunta per il tan2....  ?????????????
  
  s = s * inverse;
  return s;	    
  /*double r[3];  
  gluUnProject(p[0], p[1], p[2], model_matrix, proj_matrix, viewport, &r[0], &r[1], &r[2]);
  return Point3<T>((T)r[0], (T)r[1], (T)r[2]);*/
}

template <class T> Point3<T> View<T>::ScreenToViewport(const Point3<T> &p) const {
  Point3<T> a;
  a[0] = (p[0]+1)*(viewport[2]/(T)2.0)+viewport[0];
	a[1] =(p[1]+1)*(viewport[3]/(T)2.0)+viewport[1];
  a[1] = viewport[3] - a[1];
  a[2] = p[2];
  return a;
}
  //viceversa
template <class T> Point3<T> View<T>::ViewportToScreen(const Point3<T> &p) const {
  Point3<T> a;
  a[0] = (p[0]- viewport[0])/ (viewport[2]/(T)2.0) - 1;
	a[1] = (p[1]- viewport[1])/ (viewport[3]/(T)2.0) - 1;
  a[1] = -a[1];
	a[2] = p[2];
	//a[1] = -a[1];   // pezza aggiunta per il tan2....  ?????????????
  return a;
}


}//namespace

#endif