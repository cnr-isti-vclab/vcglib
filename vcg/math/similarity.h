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
Revision 1.8  2004/05/07 10:09:13  cignoni
missing final newline

Revision 1.7  2004/05/04 23:23:45  cignoni
unified to the gl stlyle matix*vector. removed vector*matrix operator

Revision 1.6  2004/03/25 14:57:49  ponchio
Microerror. ($LOG$ -> $Log: not supported by cvs2svn $
Microerror. ($LOG$ -> Revision 1.8  2004/05/07 10:09:13  cignoni
Microerror. ($LOG$ -> missing final newline
Microerror. ($LOG$ ->
Microerror. ($LOG$ -> Revision 1.7  2004/05/04 23:23:45  cignoni
Microerror. ($LOG$ -> unified to the gl stlyle matix*vector. removed vector*matrix operator
Microerror. ($LOG$ ->


****************************************************************************/

#ifndef SIMILARITY_H
#define SIMILARITY_H

#include <vcg/math/quaternion.h>
#include <vcg/math/matrix44.h>

namespace vcg {

template <class S> class Similarity {
public:
  Similarity() {}
  Similarity(const Quaternion<S> &q) { SetRotate(q); }
  Similarity(const Point3<S> &p) { SetTranslate(p); }
  Similarity(S s) { SetScale(s); }
  
  Similarity operator*(const Similarity &affine) const;
  Similarity &operator*=(const Similarity &affine);
  //Point3<S> operator*(const Point3<S> &p) const;
  
  
  Similarity &SetIdentity();
  Similarity &SetScale(const S s);
	Similarity &SetTranslate(const Point3<S> &t);	
  ///use radiants for angle.
  Similarity &SetRotate(S angle, const Point3<S> & axis); 
  Similarity &SetRotate(const Quaternion<S> &q);

  Matrix44<S> Matrix() const;
  Matrix44<S> InverseMatrix() const;
  void FromMatrix(const Matrix44<S> &m);

  Quaternion<S> rot;
  Point3<S> tra;
  S sca;  
};

template <class S> Similarity<S> &Invert(Similarity<S> &m);
template <class S> Similarity<S> Inverse(const Similarity<S> &m);
template <class S> Point3<S> operator*(const Similarity<S> &m, const Point3<S> &p);


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

template <class S> Similarity<S> &Similarity<S>::SetRotate(S angle, const Point3<S> &axis) {
  SetIdentity();
  rot.FromAxis(angle, axis);
  return *this;
}

template <class S> Similarity<S> &Similarity<S>::SetRotate(const Quaternion<S> &q) {
  SetIdentity();
  rot = q;  
  return *this;
}


template <class S> Matrix44<S> Similarity<S>::Matrix() const {
  Matrix44<S> r;
  rot.ToMatrix(r);
  Matrix44<S> s = Matrix44<S>().SetScale(sca, sca, sca);
  Matrix44<S> t = Matrix44<S>().SetTranslate(tra[0], tra[1], tra[2]);
  return s*r*t;  // scale * rot * trans
}

template <class S> Matrix44<S> Similarity<S>::InverseMatrix() const {
  return Inverse(Matrix());
}


template <class S> void Similarity<S>::FromMatrix(const Matrix44<S> &m) {
  sca = (S)pow(m.Determinant(), 1/3);  
  assert(sca != 0);
  Matrix44<S> t = m * Matrix44<S>().SetScale(1/sca, 1/sca, 1/sca);
  rot.FromMatrix(t);
  tra[0] = t.element(3, 0);
  tra[1] = t.element(3, 1);
  tra[2] = t.element(3, 2);
}

template <class S> Similarity<S> &Invert(Similarity<S> &a) {  
  a.rot.Invert();
  a.sca = 1/a.sca;
  a.tra = a.rot.Rotate(-a.tra)*a.sca;
  return a;
}

template <class S> Similarity<S> Inverse(const Similarity<S> &m) {
  Similarity<S> a = m;
  return Invert(a);
}


template <class S> Similarity<S> Interpolate(const Similarity<S> &a, const Similarity<S> &b, const S t) {
  Similarity<S> r;
  r.rot = interpolate(a.rot, b.rot, t);
  r.tra = t * a.tra + (1-t) * b.tra;
  r.sca = t * a.sca + (1-t) * b.sca;
  return r;
}

template <class S> Point3<S> operator*(const Similarity<S> &m, const Point3<S> &p) {
  Point3<S> r = m.rot.Rotate(p);
  r *= m.sca;
  r += m.tra;
  return r;
}

typedef Similarity<float> Similarityf;
typedef Similarity<double>Similarityd;

} //namespace

#endif
