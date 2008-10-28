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

#ifndef VCG_USE_EIGEN
#include "deprecated_point3.h"
#else

#ifndef __VCGLIB_POINT3
#define __VCGLIB_POINT3

#include "../math/eigen.h"
#include <vcg/math/base.h>

namespace vcg{
template<class Scalar> class Point3;
}

namespace Eigen{
template<typename Scalar>
struct ei_traits<vcg::Point3<Scalar> > : ei_traits<Eigen::Matrix<Scalar,3,1> > {};

template<typename Scalar>
struct NumTraits<vcg::Point3<Scalar> > : NumTraits<Scalar>
{
  enum {
    ReadCost = 3,
    AddCost = 3,
    MulCost = 3
  };
};

}

namespace vcg {

/** \addtogroup space */
/*@{*/
    /**
        The templated class for representing a point in 3D space.
        The class is templated over the ScalarType class that is used to represent coordinates. All the usual
        operator overloading (* + - ...) is present. 
     */
template <class T> class Box3;

template <class _Scalar> class Point3 : public Eigen::Matrix<_Scalar,3,1>
{
	typedef Eigen::Matrix<_Scalar,3,1> _Base;
	using _Base::coeff;
	using _Base::coeffRef;
	using _Base::setZero;
	using _Base::data;
	using _Base::V;

public:

	_EIGEN_GENERIC_PUBLIC_INTERFACE(Point3,_Base);
	typedef Scalar ScalarType;

	VCG_EIGEN_INHERIT_ASSIGNMENT_OPERATORS(Point3)
	
	enum {Dimension = 3};
	
	
//@{

  /** @name Standard Constructors and Initializers 
   No casting operators have been introduced to avoid automatic unattended (and costly) conversion between different point types
   **/

  inline Point3 () {}
	inline Point3 ( const Scalar nx, const Scalar ny, const Scalar nz ) : Base(nx,ny,nz) {}
	inline Point3 ( Point3 const & p ) : Base(p) {}
	inline Point3 ( const Scalar nv[3] ) : Base(nv) {}
	template<typename OtherDerived>
	inline Point3(const Eigen::MatrixBase<OtherDerived>& other) : Base(other) {}


	template<class OtherDerived>
	inline void Import( const Eigen::MatrixBase<OtherDerived>& b )
	{
		EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(OtherDerived,3);
		data()[0] = Scalar(b[0]);
		data()[1] = Scalar(b[1]);
		data()[2] = Scalar(b[2]);
	}

  template <class Q> 
  static inline Point3 Construct( const Point3<Q> & b )
  {
    return Point3(Scalar(b[0]),Scalar(b[1]),Scalar(b[2]));
  }

  template <class Q> 
  static inline Point3 Construct( const Q & P0, const Q & P1, const Q & P2)
  {
    return Point3(Scalar(P0),Scalar(P1),Scalar(P2));
  }

  static inline Point3 Construct( const Point3<ScalarType> & b )
  {
    return b;
  }

//@}

//@{

  /** @name Data Access. 
   access to data is done by overloading of [] or explicit naming of coords (x,y,z)**/

  inline const Scalar &X() const { return data()[0]; }
	inline const Scalar &Y() const { return data()[1]; }
	inline const Scalar &Z() const { return data()[2]; }
	inline Scalar &X() { return data()[0]; }
	inline Scalar &Y() { return data()[1]; }
	inline Scalar &Z() { return data()[2]; }
	// overloaded to return a const reference
	inline const Scalar & V( const int i ) const
	{
		assert(i>=0 && i<3);
		return data()[i];
	}
//@}
//@{

  /** @name Classical overloading of operators 
  Note   
  **/

	// Scalatura differenziata
	inline Point3 & Scale( const Scalar sx, const Scalar sy, const Scalar sz )
	{
		data()[0] *= sx;
		data()[1] *= sy;
		data()[2] *= sz;
		return *this;
	}
	inline Point3 & Scale( const Point3 & p )
	{
		data()[0] *= p.data()[0];
		data()[1] *= p.data()[1];
		data()[2] *= p.data()[2];
		return *this;
	}
	
	/** 
	 * Convert to polar coordinates from cartesian coordinates.
	 *
	 * Theta is the azimuth angle and ranges between [0, 360) degrees.
	 * Phi is the elevation angle (not the polar angle) and ranges between [-90, 90] degrees.
	 *
	 * /note Note that instead of the classical polar angle, which ranges between 
	 *       0 and 180 degrees we opt for the elevation angle to obtain a more 
	 *       intuitive spherical coordinate system.
	 */
	void ToPolar(Scalar &ro, Scalar &theta, Scalar &phi) const
	{
		ro = this->norm();
		theta = (Scalar)atan2(data()[2], data()[0]);
		phi   = (Scalar)asin(data()[1]/ro);
	}

	/**
	 * Convert from polar coordinates to cartesian coordinates.
	 *
	 * Theta is the azimuth angle and ranges between [0, 360) degrees.
	 * Phi is the elevation angle (not the polar angle) and ranges between [-90, 90] degrees.
	 *
	 * \note Note that instead of the classical polar angle, which ranges between 
	 *       0 and 180 degrees, we opt for the elevation angle to obtain a more 
	 *       intuitive spherical coordinate system.
	 */
  void FromPolar(const Scalar &ro, const Scalar &theta, const Scalar &phi)
	{
    data()[0]= ro*cos(theta)*cos(phi);
    data()[1]= ro*sin(phi);
    data()[2]= ro*sin(theta)*cos(phi);
	}
	
  Box3<_Scalar> GetBBox(Box3<_Scalar> &bb) const;
//@}

}; // end class definition


// versione uguale alla precedente ma che assume che i due vettori sono unitari
template <class Scalar>
inline Scalar AngleN( Point3<Scalar> const & p1, Point3<Scalar> const & p2 )
{
	Scalar w = p1*p2;
	if(w>1) 
		w = 1;
	else if(w<-1) 
		w=-1;
  return (Scalar) acos(w);
}


template <class Scalar>
inline Point3<Scalar> & Normalize( Point3<Scalar> & p )
{
    p.Normalize();
    return p;
}

	// Dot product preciso numericamente (solo double!!)
	// Implementazione: si sommano i prodotti per ordine di esponente
	// (prima le piu' grandi)
template<class Scalar>
double stable_dot ( Point3<Scalar> const & p0, Point3<Scalar> const & p1 )
{
	Scalar k0 = p0.data()[0]*p1.data()[0];
	Scalar k1 = p0.data()[1]*p1.data()[1];
	Scalar k2 = p0.data()[2]*p1.data()[2];

	int exp0,exp1,exp2;

	frexp( double(k0), &exp0 );
	frexp( double(k1), &exp1 );
	frexp( double(k2), &exp2 );

	if( exp0<exp1 )
	{
		if(exp0<exp2)
			return (k1+k2)+k0;
		else
			return (k0+k1)+k2;
	}
	else
	{
		if(exp1<exp2)
			return(k0+k2)+k1;
		else
			return (k0+k1)+k2;
	}
}  



/// Point(p) Edge(v1-v2) dist, q is the point in v1-v2 with min dist
template<class Scalar>
Scalar PSDist( const Point3<Scalar> & p,
			         const Point3<Scalar> & v1,
					 const Point3<Scalar> & v2,
			         Point3<Scalar> & q )
{
    Point3<Scalar> e = v2-v1;
    Scalar  t = ((p-v1).dot(e))/e.SquaredNorm();
    if(t<0)      t = 0;
	else if(t>1) t = 1;
	q = v1+e*t;
    return Distance(p,q);
}


template <class Scalar>
void GetUV( Point3<Scalar> &n,Point3<Scalar> &u, Point3<Scalar> &v, Point3<Scalar> up=(Point3<Scalar>(0,1,0)) )
{
	n.Normalize();
	const double LocEps=double(1e-7);
	u=n^up;
  double len = u.Norm();
 	if(len < LocEps)
	{
		if(fabs(n[0])<fabs(n[1])){
			if(fabs(n[0])<fabs(n[2])) up=Point3<Scalar>(1,0,0); // x is the min
			                         else up=Point3<Scalar>(0,0,1); // z is the min
		}else {
			if(fabs(n[1])<fabs(n[2])) up=Point3<Scalar>(0,1,0); // y is the min
			                         else up=Point3<Scalar>(0,0,1); // z is the min
		}
		u=n^up;
	}
	u.Normalize();
	v=n^u;
	v.Normalize();
	Point3<Scalar> uv=u^v;
}


template <class SCALARTYPE>
inline Point3<SCALARTYPE> Abs(const Point3<SCALARTYPE> & p) {
	return (Point3<SCALARTYPE>(math::Abs(p[0]), math::Abs(p[1]), math::Abs(p[2])));
}
// probably a more uniform naming should be defined...
template <class SCALARTYPE>
inline Point3<SCALARTYPE> LowClampToZero(const Point3<SCALARTYPE> & p) {
	return (Point3<SCALARTYPE>(math::Max(p[0], (SCALARTYPE)0), math::Max(p[1], (SCALARTYPE)0), math::Max(p[2], (SCALARTYPE)0)));
}

typedef Point3<short>  Point3s;
typedef Point3<int>	   Point3i;
typedef Point3<float>  Point3f;
typedef Point3<double> Point3d;

/*@}*/

} // end namespace

#endif

#endif
