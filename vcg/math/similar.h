#ifndef SIMILAR_H
#define SIMILAR_H

#include <vcg/math/quaternion.h>
#include <vcg/math/matrix44.h>

namespace vcg {

template <class S> class Similar {
public:
  Similar() {}
  
  Similar operator*(const Similar &affine) const;
  Similar &operator*=(const Similar &affine);
  Point3<S> operator*(const Point3<S> &p) const;
  
  
  void SetIdentity();
  void SetScale(const S s);
	Similar &SetTranslate(const Point3<S> &t);	
  ///use radiants for angle.
  void SetRotate(S angle, const Point3<S> & axis); 

  Matrix44<S> Matrix() const;

  Quaternion<S> rot;
  Point3<S> tra;
  S sca;  
};

template <class S> Similar<S> Similar<S>::operator*(const Similar &a) const {
  Similar<S> r;
  r.rot = rot * a.rot;
  r.sca = sca * a.sca;
  r.tra = (rot.Rotate(a.tra)) * sca + tra;
  return r;
}

template <class S> Similar<S> &Similar<S>::operator*=(const Similar &a) {  
  rot = rot * a.rot;
  sca = sca * a.sca;
  tra = (rot.Rotate(a.tra)) * sca + tra;
  return *this;
}

template <class S> Point3<S> Similar<S>::operator*(const Point3<S> &p) const {
  Point3<S> r = rot.Rotate(p);
  r *= sca;
  r += tra;
  return r;
}
  
template <class S> void Similar<S>::SetIdentity() {
  rot.FromAxis(0, Point3<S>(1, 0, 0));
  tra = Point3<S>(0, 0, 0);
  sca = 1;
}

template <class S> void Similar<S>::SetScale(const S s) {
  SetIdentity();
  sca = s;
}

template <class S> Similar<S> &Similar<S>::SetTranslate(const Point3<S> &t) {
  SetIdentity();
  tra = t;
  return *this;
}

template <class S> void Similar<S>::SetRotate(S angle, const Point3<S> &axis) {
  SetIdentity();
  rot.FromAxis(angle, axis);
}

template <class S> Matrix44<S> Similar<S>::Matrix() const {
  Matrix44<S> r;
  rot.ToMatrix(r);
  r *= Matrix44<S>().SetScale(sca, sca, sca);
  Matrix44<S> t = Matrix44<S>().SetTranslate(tra[0], tra[1], tra[2]);
  r *= t;
  return r;
}

template <class S> Similar<S> &invert(Similar<S> &a) {
  //WARNING:: TEST THIS!!!
  a.rot.Invert();
  a.sca = 1/a.sca;
  a.tra = a.rot.Rotate(-a.tra)*a.sca;
  return a;
}

template <class S> Similar<S> interpolate(const Similar<S> &a, const Similar<S> &b, const S t) {
  Similar<S> r;
  r.rot = interpolate(a.rot, b.rot, t);
  r.tra = t * a.tra + (1-t) * b.tra;
  r.sca = t * a.sca + (1-t) * b.sca;
  return r;
}

typedef Similar<float> Similarf;
typedef Similar<double>Similard;

} //namespace

#endif