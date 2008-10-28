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
#include "deprecated_point2.h"
#else

#ifndef __VCGLIB_POINT2
#define __VCGLIB_POINT2

#include "../math/eigen.h"
#include <vcg/math/base.h>

namespace vcg{
template<class Scalar> class Point2;
}

namespace Eigen{
template<typename Scalar>
struct ei_traits<vcg::Point2<Scalar> > : ei_traits<Eigen::Matrix<Scalar,2,1> > {};
}

namespace vcg {

/** \addtogroup space */
/*@{*/
    /**
        The templated class for representing a point in 2D space.
        The class is templated over the Scalar class that is used to represent coordinates.
				All the usual operator overloading (* + - ...) is present.
     */
template <class _Scalar> class Point2 : public Eigen::Matrix<_Scalar,2,1>
{
	typedef Eigen::Matrix<_Scalar,2,1> _Base;
	using _Base::coeff;
	using _Base::coeffRef;
	using _Base::setZero;
	using _Base::data;
	using _Base::V;

public:

	_EIGEN_GENERIC_PUBLIC_INTERFACE(Point2,_Base);
	typedef Scalar ScalarType;

	VCG_EIGEN_INHERIT_ASSIGNMENT_OPERATORS(Point2)

	enum {Dimension = 2};

//@{
  /** @name Access to Coords.
   access to coords is done by overloading of [] or explicit naming of coords (X,Y,)
	 ("p[0]" or "p.X()" are equivalent) **/
	inline const Scalar &X() const {return data()[0];}
	inline const Scalar &Y() const {return data()[1];}
	inline Scalar &X() {return data()[0];}
	inline Scalar &Y() {return data()[1];}

	// overloaded to return a const reference
	inline const Scalar & V( const int i ) const
	{
		assert(i>=0 && i<2);
		return data()[i];
	}
//@}

	/// empty constructor (does nothing)
	inline Point2 () { }
	/// x,y constructor
	inline Point2 ( const Scalar nx, const Scalar ny ) : Base(nx,ny) {}
	/// copy constructor
	inline Point2(Point2 const & p) : Base(p) {}
	template<typename OtherDerived>
	inline Point2(const Eigen::MatrixBase<OtherDerived>& other) : Base(other) {}

	/// cross product
	// hm.. this is not really a cross product
	inline Scalar operator ^ ( Point2 const & p ) const
	{
		return data()[0]*p.data()[1] - data()[1]*p.data()[0];
	}

	/// returns the angle with X axis (radiants, in [-PI, +PI] )
	inline Scalar Angle() const
	{
		return math::Atan2(data()[1],data()[0]);
	}
		/// transform the point in cartesian coords into polar coords
	inline Point2 & Cartesian2Polar()
	{
		Scalar t = Angle();
		data()[0] = this->norm();
		data()[1] = t;
		return *this;
	}
		/// transform the point in polar coords into cartesian coords
	inline Point2 & Polar2Cartesian()
	{
		Scalar l = data()[0];
		data()[0] = (Scalar)(l*math::Cos(data()[1]));
		data()[1] = (Scalar)(l*math::Sin(data()[1]));
		return *this;
	}
		/// rotates the point of an angle (radiants, counterclockwise)
	inline Point2 & Rotate( const Scalar rad )
	{
		Scalar t = data()[0];
		Scalar s = math::Sin(rad);
		Scalar c = math::Cos(rad);

		data()[0] = data()[0]*c - data()[1]*s;
		data()[1] =   t *s + data()[1]*c;

		return *this;
	}

		/// imports from 2D points of different types
	template <class T>
	inline void Import( const Point2<T> & b )
	{
		data()[0] = b.X(); data()[1] = b.Y();
	}
		/// constructs a 2D points from an existing one of different type
	template <class T>
	static Point2 Construct( const Point2<T> & b )
	{
    return Point2(b.X(),b.Y());
	}
}; // end class definition

typedef Point2<short>  Point2s;
typedef Point2<int>	   Point2i;
typedef Point2<float>  Point2f;
typedef Point2<double> Point2d;

/*@}*/
} // end namespace
#endif

#endif
