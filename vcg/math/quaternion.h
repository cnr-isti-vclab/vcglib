


#ifndef QUATERNION_H
#define QUATERNION_H

#include <vcg/space/point3.h>
#include <vcg/space/point4.h>
#include <vcg/math/base.h>
#include <vcg/math/matrix44.h>

namespace vcg {

	/** Classe quaternion.
	    A quaternion is a point in the unit sphere in four dimension: all
			rotations in three-dimensional space can be represented by a quaternion.
		*/
template<class S> class Quaternion: public Point4<S> {
public:

	Quaternion() {}
	Quaternion(const S v0, const S v1, const S v2, const S v3): Point4<S>(v0,v1,v2,v3){}	
	Quaternion(const Point4<S> p) : Point4<S>(p)	{}
  Quaternion(const S phi, const Point3<S> &a);

  Quaternion operator*(const S &s) const;
  //Quaternion &operator*=(S d);
  Quaternion operator*(const Quaternion &q) const;
  Quaternion &operator*=(const Quaternion &q);
  void Invert();

  void FromAxis(const S phi, const Point3<S> &a);
  void ToAxis(S &phi, Point3<S> &a ) const;

  void FromMatrix(Matrix44<S> &m);
  void ToMatrix(Matrix44<S> &m) const;
  
  Point3<S> Rotate(const Point3<S> vec) const;  
};

template <class S> Quaternion<S> interpolate(const Quaternion<S> a, const Quaternion<S> b, double t);


//Implementation
	
template <class S> Quaternion<S>::Quaternion(const S phi, const Point3<S> &a) {
  FromAxis(phi, a);
}
 	 

template <class S> Quaternion<S> Quaternion<S>::operator*(const S &s) const {
		return (Quaternion(V(0)*s,V(1)*s,V(2)*s,V(3)*s));
}
 
template <class S> Quaternion<S> Quaternion<S>::operator*(const Quaternion &q) const {		
	Point3<S> t1(V(1), V(2), V(3));
  Point3<S> t2(q.V(1), q.V(2), q.V(3));    
		
  S d  = t2 * t1;
	Point3<S> t3 = t1 ^ t2;
		
  t1 *= q.V(0);
	t2 *= V(0);

	Point3<S> tf = t1 + t2 +t3;

   Quaternion<S> t;
	t.V(0) = V(0) * q.V(0) - d;
	t.V(1) = tf[0];
	t.V(2) = tf[1];
	t.V(3) = tf[2];
	return t;
}

template <class S> Quaternion<S> &Quaternion<S>::operator*=(const Quaternion &q) {
  S ww = V(0) * q.V(0) - V(1) * q.V(1) - V(2) * q.V(2) - V(3) * q.V(3);
	S xx = V(0) * q.V(1) + V(1) * q.V(0) + V(2) * q.V(3) - V(3) * q.V(2);
	S yy = V(0) * q.V(2) - V(1) * q.V(3) + V(2) * q.V(0) + V(3) * q.V(1);
	    	
	V(0) = ww; 
  V(1) = xx; 
  V(2) = yy;
  V(3) = V(0) * q.V(3) + V(1) * q.V(2) - V(2) * q.V(1) + V(3) * q.V(0);
	return *this;
}	 

template <class S> void Quaternion<S>::Invert() {
	V(1)*=-1;
	V(2)*=-1;
	V(3)*=-1;
}


template <class S> void Quaternion<S>::FromAxis(const S phi, const Point3<S> &a) {
  Point3<S> b = a;
  b.Normalize();
  S s = math::Sin(phi/(S(2.0)));

  V(0) = math::Cos(phi/(S(2.0)));
	V(1) = b[0]*s;
	V(2) = b[1]*s;
	V(3) = b[2]*s;
}

template <class S> void Quaternion<S>::ToAxis(S &phi, Point3<S> &a) const {
  S s = math::Asin(V(0))*S(2.0);
  phi = math::Acos(V(0))*S(2.0);

	if(s < 0) 
		phi = - phi;

  a.V(0) = V(1);
	a.V(1) = V(2);
	a.V(2) = V(3);
  a.Normalize();
}


template <class S> Point3<S> Quaternion<S>::Rotate(const Point3<S> p) const {
		Quaternion<S> co = *this;
		co.Invert();

    Quaternion<S> tmp(0, p.V(0), p.V(1), p.V(2));		

		tmp = (*this) * tmp * co;
		return 	Point3<S>(tmp.V(1), tmp.V(2), tmp.V(3));
	}


template <class S> void Quaternion<S>::ToMatrix(Matrix44<S> &m) const	{
	S q00 = V(1)*V(1);
	S q01 = V(1)*V(2);
	S q02 = V(1)*V(3);
	S q03 = V(1)*V(0);
	S q11 = V(2)*V(2);
	S q12 = V(2)*V(3);
	S q13 = V(2)*V(0);
	S q22 = V(3)*V(3);
	S q23 = V(3)*V(0);

  m.element(0, 0) = (S)(1.0-(q11 + q22)*2.0);
  m.element(1, 0) = (S)((q01 - q23)*2.0);
  m.element(2, 0) = (S)((q02 + q13)*2.0);
  m.element(3, 0) = (S)0.0;

  m.element(0, 1) = (S)((q01 + q23)*2.0);
  m.element(1, 1) = (S)(1.0-(q22 + q00)*2.0);
  m.element(2, 1) = (S)((q12 - q03)*2.0);
  m.element(3, 1) = (S)0.0;

  m.element(0, 2) = (S)((q02 - q13)*2.0);
  m.element(1, 2) = (S)((q12 + q03)*2.0);
  m.element(2, 2) = (S)(1.0-(q11 + q00)*2.0);
  m.element(3, 2) = (S)0.0;

  m.element(0, 3) = (S)0.0;
  m.element(1, 3) = (S)0.0;
  m.element(2, 3) = (S)0.0;
  m.element(3, 3) = (S)1.0;
}
	

///warning m deve essere una matrice di rotazione pena il disastro.
template <class S> void Quaternion<S>::FromMatrix(Matrix44<S> &m) {		
	double Sc;
	double t = (m[0] + m[5] + m[10] + 1);	
	if(t > 0) {
      Sc = 0.5 / sqrt(t);
      V(0) = 0.25 / Sc;
      V(1) = ( m[9] - m[6] ) * Sc;
      V(2) = ( m[2] - m[8] ) * Sc;
      V(3) = ( m[4] - m[1] ) * Sc;
	}	else {	
		if(m[0] > m[5] && m[0] > m[10]) {
			Sc  = sqrt( 1.0 + m[0] - m[5] - m[10] ) * 2;
			V(1) = 0.5 / Sc;
			V(2) = (m[1] + m[4] ) / Sc;
			V(3) = (m[2] + m[8] ) / Sc;
			V(0) = (m[6] + m[9] ) / Sc;	
		}	else if( m[5] > m[10]) {
			Sc  = sqrt( 1.0 + m[5] - m[0] - m[10] ) * 2;
			V(1) = (m[1] + m[4] ) / Sc;
			V(2) = 0.5 / Sc;
			V(3) = (m[6] + m[9] ) / Sc;
			V(0) = (m[2] + m[8] ) / Sc;
		}	else {
			Sc  = sqrt( 1.0 + m[10] - m[0] - m[5] ) * 2;
			V(1) = (m[2] + m[8] ) / Sc;
			V(2) = (m[6] + m[9] ) / Sc;
			V(3) = 0.5 / Sc;
			V(0) = (m[1] + m[4] ) / Sc;
		}
	}
}
template <class S> Quaternion<S> interpolate(const Quaternion<S> a, const Quaternion<S> b, double t) {
		double v = a.V(0) * b.V(0) + a.V(1) * b.V(1) + a.V(2) * b.V(2) + a.V(3) * b.V(3);
		double phi = Acos(v);
		if(phi > 0.01) {
			a = a * (Sin(phi *(1-t))/Sin(phi));
			b = b * (Sin(phi * t)/Sin(phi));	
		}
		
		Quaternion<S> c;
		c.V(0) = a.V(0) + b.V(0);
		c.V(1) = a.V(1) + b.V(1);
		c.V(2) = a.V(2) + b.V(2);
		c.V(3) = a.V(3) + b.V(3);
		
		if(v < -0.999) { //almost opposite
			double d = t * (1 - t);
			if(c.V(0) == 0)
				c.V(0) += d;
			else
				c.V(1) += d;		
		}
		c.Normalize();
		return c;
	}
		
	

typedef Quaternion<float>  Quaternionf;
typedef Quaternion<double> Quaterniond;

} // end namespace


#endif
