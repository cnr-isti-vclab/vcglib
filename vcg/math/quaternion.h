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
Revision 1.13  2005/04/15 09:19:50  ponchio
Typo: Point3 -> Point4

Revision 1.12  2005/04/14 17:22:34  ponchio
*** empty log message ***

Revision 1.11  2005/04/14 11:35:09  ponchio
*** empty log message ***

Revision 1.10  2004/12/15 18:45:50  tommyfranken
*** empty log message ***

Revision 1.9  2004/10/22 14:35:42  ponchio
m.element(x, y) -> m[x][y]

Revision 1.8  2004/10/07 13:54:03  ganovelli
added SetIdentity

Revision 1.7  2004/04/07 10:48:37  cignoni
updated access to matrix44 elements through V() instead simple []

Revision 1.6  2004/03/25 14:57:49  ponchio
Microerror. ($LOG$ -> $Log: not supported by cvs2svn $
Microerror. ($LOG$ -> Revision 1.13  2005/04/15 09:19:50  ponchio
Microerror. ($LOG$ -> Typo: Point3 -> Point4
Microerror. ($LOG$ ->
Microerror. ($LOG$ -> Revision 1.12  2005/04/14 17:22:34  ponchio
Microerror. ($LOG$ -> *** empty log message ***
Microerror. ($LOG$ ->
Microerror. ($LOG$ -> Revision 1.11  2005/04/14 11:35:09  ponchio
Microerror. ($LOG$ -> *** empty log message ***
Microerror. ($LOG$ ->
Microerror. ($LOG$ -> Revision 1.10  2004/12/15 18:45:50  tommyfranken
Microerror. ($LOG$ -> *** empty log message ***
Microerror. ($LOG$ ->
Microerror. ($LOG$ -> Revision 1.9  2004/10/22 14:35:42  ponchio
Microerror. ($LOG$ -> m.element(x, y) -> m[x][y]
Microerror. ($LOG$ ->
Microerror. ($LOG$ -> Revision 1.8  2004/10/07 13:54:03  ganovelli
Microerror. ($LOG$ -> added SetIdentity
Microerror. ($LOG$ ->
Microerror. ($LOG$ -> Revision 1.7  2004/04/07 10:48:37  cignoni
Microerror. ($LOG$ -> updated access to matrix44 elements through V() instead simple []
Microerror. ($LOG$ ->


****************************************************************************/


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
  
	
	void SetIdentity();
	

  void FromAxis(const S phi, const Point3<S> &a);
  void ToAxis(S &phi, Point3<S> &a ) const;

  void FromMatrix(Matrix44<S> &m);
  void ToMatrix(Matrix44<S> &m) const;
  
  Point3<S> Rotate(const Point3<S> vec) const; 
  //duplicated ... because of gcc new confoming to ISO template derived classes
  //do no 'see' parent members (unless explicitly specified) 
  const S & V ( const int i ) const	{ assert(i>=0 && i<4); return Point4<S>::V(i); }
  S & V ( const int i )	{ assert(i>=0 && i<4); return Point4<S>::V(i); }
};

template <class S> Quaternion<S> Interpolate(  Quaternion<S>   a,   Quaternion<S>   b, double t);
template <class S> Quaternion<S> &Invert(Quaternion<S> &q);
template <class S> Quaternion<S> Inverse(const Quaternion<S> &q);


//Implementation
template <class S> 
void Quaternion<S>::SetIdentity(){
	FromAxis(0, Point3<S>(1, 0, 0));
}
	

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

/*	<<<<<<< quaternion.h
  m[0][ 0] = (S)(1.0-(q11 + q22)*2.0);
  m[1][ 0] = (S)((q01 - q23)*2.0);
  m[2][ 0] = (S)((q02 + q13)*2.0);
  m[3][ 0] = (S)0.0;

  m[0][ 1] = (S)((q01 + q23)*2.0);
  m[1][ 1] = (S)(1.0-(q22 + q00)*2.0);
  m[2][ 1] = (S)((q12 - q03)*2.0);
  m[3][ 1] = (S)0.0;

  m[0][ 2] = (S)((q02 - q13)*2.0);
  m[1][ 2] = (S)((q12 + q03)*2.0);
  m[2][ 2] = (S)(1.0-(q11 + q00)*2.0);
  m[3][ 2] = (S)0.0;

  m[0][ 3] = (S)0.0;
  m[1][ 3] = (S)0.0;
  m[2][ 3] = (S)0.0;
  m[3][ 3] = (S)1.0;
=======*/
  m[0][0] = (S)(1.0-(q11 + q22)*2.0);
  m[1][0] = (S)((q01 - q23)*2.0);
  m[2][0] = (S)((q02 + q13)*2.0);
  m[3][0] = (S)0.0;

  m[0][1] = (S)((q01 + q23)*2.0);
  m[1][1] = (S)(1.0-(q22 + q00)*2.0);
  m[2][1] = (S)((q12 - q03)*2.0);
  m[3][1] = (S)0.0;

  m[0][2] = (S)((q02 - q13)*2.0);
  m[1][2] = (S)((q12 + q03)*2.0);
  m[2][2] = (S)(1.0-(q11 + q00)*2.0);
  m[3][2] = (S)0.0;

  m[0][3] = (S)0.0;
  m[1][3] = (S)0.0;
  m[2][3] = (S)0.0;
  m[3][3] = (S)1.0;
//>>>>>>> 1.9
}
	

///warning m deve essere una matrice di rotazione pena il disastro.
template <class S> void Quaternion<S>::FromMatrix(Matrix44<S> &m) {		
	S Sc;
	S t = (m.V()[0] + m.V()[5] + m.V()[10] + (S)1.0);	
	if(t > 0) {
      Sc = (S)0.5 / math::Sqrt(t);
      V(0) = (S)0.25 / Sc;
      V(1) = ( m.V()[9] - m.V()[6] ) * Sc;
      V(2) = ( m.V()[2] - m.V()[8] ) * Sc;
      V(3) = ( m.V()[4] - m.V()[1] ) * Sc;
	}	else {	
		if(m.V()[0] > m.V()[5] && m.V()[0] > m.V()[10]) {
      Sc  = math::Sqrt( (S)1.0 + m.V()[0] - m.V()[5] - m.V()[10] ) * (S)2.0;
			V(1) = (S)0.5 / Sc;
			V(2) = (m.V()[1] + m.V()[4] ) / Sc;
			V(3) = (m.V()[2] + m.V()[8] ) / Sc;
			V(0) = (m.V()[6] + m.V()[9] ) / Sc;	
		}	else if( m.V()[5] > m.V()[10]) {
      Sc  = math::Sqrt( (S)1.0 + m.V()[5] - m.V()[0] - m.V()[10] ) * (S)2.0;
			V(1) = (m.V()[1] + m.V()[4] ) / Sc;
			V(2) = (S)0.5 / Sc;
			V(3) = (m.V()[6] + m.V()[9] ) / Sc;
			V(0) = (m.V()[2] + m.V()[8] ) / Sc;
		}	else {
      Sc  = math::Sqrt( (S)1.0 + m.V()[10] - m.V()[0] - m.V()[5] ) * (S)2.0;
			V(1) = (m.V()[2] + m.V()[8] ) / Sc;
			V(2) = (m.V()[6] + m.V()[9] ) / Sc;
			V(3) = (S)0.5 / Sc;
			V(0) = (m.V()[1] + m.V()[4] ) / Sc;
		}
	}
}
template <class S> Quaternion<S> &Invert(Quaternion<S> &m) {
  m.Invert();
  return m;
}

template <class S> Quaternion<S> Inverse(const Quaternion<S> &m) {
  Quaternion<S> a = m;
  a.Invert();
  return a;
}

template <class S> Quaternion<S> Interpolate(   Quaternion<S>   a ,    Quaternion<S>   b , double t) {
 	 
		double v = a.V(0) * b.V(0) + a.V(1) * b.V(1) + a.V(2) * b.V(2) + a.V(3) * b.V(3);
		double phi = math::Acos(v);
		if(phi > 0.01) {
			a = a  * (math::Sin(phi *(1-t))/math::Sin(phi));
			b = b  * (math::Sin(phi * t)/math::Sin(phi));	
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
