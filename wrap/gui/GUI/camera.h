#ifndef CAMERA_H
#define CAMERA_H

#include <vcg/space/point3.h>
#include <vcg/math/Matrix44.h>

namespace vcg {

template <class T> class ViewGL {
public:
  void GetView();
  void SetView();                      
  Point3<T> Project(const Point3<T> &p) const;
  Point3<T> UnProject(const Point3<T> &p) const;
  Point3<T> ViewPoint();

  Matrix44<T> proj;
  Matrix44<T> model;
  Matrix44<T> matrix;
  Matrix44<T> inverse;  
  int viewport[4];    
};

template <class T> void Camera<T>::GetView() {
  glGetDoublev(GL_PROJECTION_MATRIX, (double *)&proj);
	glGetDoublev(GL_MODELVIEW_MATRIX, (double *)&model);
	glGetIntegerv(GL_VIEWPORT, viewport);

  proj.Transpose();
  model.Transpose();

  matrix.Import(proj * model);
  inverse = matrix;
  Invert(inverse);
}

template <class T> Point3<T> Camera<T>::ViewPoint() {
  return inverse * Point3<T>(0, 0, 0);
  /*Matrix44d model(model_matrix);
	model.Invert();
  Point3d view = model * Point3d(0, 0, 0);	   
	return Point3<T>(view[0], view[1], view[2]); */
}

template <class T> Point3<T> Camera<T>::Project(const Point3<T> &p) const {
  Point3<T> r;
  r = matrix * p;  	
	r[0] = (r[0]+1)*(viewport[2]/2.0)+viewport[0];
	r[1] =(r[1]+1)*(viewport[3]/2.0)+viewport[1];
  return r;
  /*double r[3];
  gluProject(p[0], p[1], p[2], model_matrix, proj_matrix, viewport, &r[0], &r[1], &r[2]);
  return Point3<T>((T)r[0], (T)r[1], (T)r[2]);*/
}

template <class T> Point3<T> Camera<T>::UnProject(const Point3<T> &p) const {
  Point3<T> s;	
	s[0] = (p[0]- viewport[0])/ (viewport[2]/2.0) - 1;
	s[1] = (p[1]- viewport[1])/ (viewport[3]/2.0) - 1;
	s[2] = p[2];
	s[1] = -s[1];   // pezza aggiunta per il tan2....  ?????????????
  s = inverse * s; 
  return s;	    
  /*double r[3];  
  gluUnProject(p[0], p[1], p[2], model_matrix, proj_matrix, viewport, &r[0], &r[1], &r[2]);
  return Point3<T>((T)r[0], (T)r[1], (T)r[2]);*/
}


}//namespace

#endif