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

#ifndef EIGEN_VCGLIB
#define EIGEN_VCGLIB

#define EIGEN_MATRIXBASE_PLUGIN <vcg/math/eigen_vcgaddons.h>

#include "../Eigen/LU"
#include "../Eigen/Geometry"
#include "../Eigen/Array"
#include "../Eigen/Core"
#include "base.h"

namespace Eigen {
template<> struct NumTraits<unsigned char>
{
  typedef unsigned char Real;
  typedef float FloatingPoint;
  enum {
    IsComplex = 0,
    HasFloatingPoint = 0,
    ReadCost = 1,
    AddCost = 1,
    MulCost = 1
  };
};

template<> struct NumTraits<short int>
{
  typedef short int Real;
  typedef float FloatingPoint;
  enum {
    IsComplex = 0,
    HasFloatingPoint = 0,
    ReadCost = 1,
    AddCost = 1,
    MulCost = 1
  };
};

}

#define VCG_EIGEN_INHERIT_ASSIGNMENT_OPERATOR(Derived, Op) \
	template<typename OtherDerived> \
	Derived& operator Op(const Eigen::MatrixBase<OtherDerived>& other) \
	{ \
		Base::operator Op(other.derived());  return *this;\
	} \
	Derived& operator Op(const Derived& other) \
	{ \
		Base::operator Op(other); return *this;\
	}

	#define VCG_EIGEN_INHERIT_SCALAR_ASSIGNMENT_OPERATOR(Derived, Op) \
	template<typename Other> \
	Derived& operator Op(const Other& scalar) \
	{ \
		Base::operator Op(scalar); return *this;\
	}

#define VCG_EIGEN_INHERIT_ASSIGNMENT_OPERATORS(Derived) \
	VCG_EIGEN_INHERIT_ASSIGNMENT_OPERATOR(Derived, =) \
	VCG_EIGEN_INHERIT_ASSIGNMENT_OPERATOR(Derived, +=) \
	VCG_EIGEN_INHERIT_ASSIGNMENT_OPERATOR(Derived, -=) \
	VCG_EIGEN_INHERIT_SCALAR_ASSIGNMENT_OPERATOR(Derived, *=) \
	VCG_EIGEN_INHERIT_SCALAR_ASSIGNMENT_OPERATOR(Derived, /=)


namespace vcg {

template<typename Derived1, typename Derived2>
typename Eigen::ei_traits<Derived1>::Scalar
Angle(const Eigen::MatrixBase<Derived1>& p1, const Eigen::MatrixBase<Derived2> & p2)
{
	EIGEN_STATIC_ASSERT_VECTOR_ONLY(Derived1)
	EIGEN_STATIC_ASSERT_VECTOR_ONLY(Derived2)
	EIGEN_STATIC_ASSERT_FIXED_SIZE(Derived1)
	EIGEN_STATIC_ASSERT_FIXED_SIZE(Derived2)
	EIGEN_STATIC_ASSERT_SAME_VECTOR_SIZE(Derived1,Derived2)
	typedef typename Eigen::ei_traits<Derived1>::Scalar Scalar;

	Scalar w = p1.norm()*p2.norm();
	if(w==0) return Scalar(-1);
	Scalar t = (p1.dot(p2))/w;
	if(t>1) t = 1;
	else if(t<-1) t = -1;
	return vcg::math::Acos(t);
}

template<typename Derived1>
inline typename Eigen::ei_traits<Derived1>::Scalar Norm( const Eigen::MatrixBase<Derived1>& p)
{ return p.norm(); }

template<typename Derived1>
inline typename Eigen::ei_traits<Derived1>::Scalar SquaredNorm( const Eigen::MatrixBase<Derived1>& p)
{ return p.norm2(); }

template<typename Derived1, typename Derived2>
inline typename Eigen::ei_traits<Derived1>::Scalar
Distance(const Eigen::MatrixBase<Derived1>& p1, const Eigen::MatrixBase<Derived2> & p2)
{ return (p1-p2).norm(); }

template<typename Derived1, typename Derived2>
inline typename Eigen::ei_traits<Derived1>::Scalar
SquaredDistance(const Eigen::MatrixBase<Derived1>& p1, const Eigen::MatrixBase<Derived2> & p2)
{ return (p1-p2).norm2(); }

}

#endif
