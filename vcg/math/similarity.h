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

$LOG$

****************************************************************************/

#ifndef SIMILARITY_H
#define SIMILARITY_H

#include <vcg/math/quaternion.h>
#include <vcg/math/matrix44.h>

namespace vcg {

template <class S> class Similarity {
public:
  Similarity() {}
  
  Similarity operator*(const Similarity &affine) const;
  Similarity &operator*=(const Similarity &affine);
  Point3<S> operator*(const Point3<S> &p) const;
  
  
  Similarity &SetIdentity();
  Similarity &SetScale(const S s);
	Similarity &SetTranslate(const Point3<S> &t);	
  ///use radiants for angle.
  void SetRotate(S angle, const Point3<S> & axis); 

  Matrix44<S> Matrix() const;

  Quaternion<S> rot;
  Point3<S> tra;
  S sca;  
};

template <class S> Similarity<S> Similarity<S>::operator*(const Similarity &a) const {
  Similarity<S> r;
  r.rot = rot * a.rot;
  r.sca = sca * a.sca;
  r.tra = (rot.Rotate(a.tra)) * sca + tra;
  return r;
}

template <class S> Similarity<S> &Similarity<S>::operator*=(const Similarity &a) {  
  rot = rot * a.rot;
  sca = sca * a.sca;
  tra = (rot.Rotate(a.tra)) * sca + tra;
  return *this;
}

template <class S> Point3<S> Similarity<S>::operator*(const Point3<S> &p) const {
  Point3<S> r = rot.Rotate(p);
  r *= sca;
  r += tra;
  return r;
}
  
template <class S> Similarity<S> &Similarity<S>::SetIdentity() {
  rot.FromAxis(0, Point3<S>(1, 0, 0));
  tra = Point3<S>(0, 0, 0);
  sca = 1;
  return *this;
}

template <class S> Similarity<S> &Similarity<S>::SetScale(const S s) {
  SetIdentity();
  sca = s;
  return *this;
}

template <class S> Similarity<S> &Similarity<S>::SetTranslate(const Point3<S> &t) {
  SetIdentity();
  tra = t;
  return *this;
}

template <class S> void Similarity<S>::SetRotate(S angle, const Point3<S> &axis) {
  SetIdentity();
  rot.FromAxis(angle, axis);
}

template <class S> Matrix44<S> Similarity<S>::Matrix() const {
  Matrix44<S> r;
  rot.ToMatrix(r);
  r *= Matrix44<S>().SetScale(sca, sca, sca);
  Matrix44<S> t = Matrix44<S>().SetTranslate(tra[0], tra[1], tra[2]);
  r *= t;
  return r;
}

template <class S> Similarity<S> &invert(Similarity<S> &a) {
  //WARNING:: TEST THIS!!!
  a.rot.Invert();
  a.sca = 1/a.sca;
  a.tra = a.rot.Rotate(-a.tra)*a.sca;
  return a;
}

template <class S> Similarity<S> interpolate(const Similarity<S> &a, const Similarity<S> &b, const S t) {
  Similarity<S> r;
  r.rot = interpolate(a.rot, b.rot, t);
  r.tra = t * a.tra + (1-t) * b.tra;
  r.sca = t * a.sca + (1-t) * b.sca;
  return r;
}

typedef Similarity<float> Similarityf;
typedef Similarity<double>Similarityd;

} //namespace

#endif