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

// TODO enable the vectorization
#define EIGEN_DONT_VECTORIZE
#define EIGEN_MATRIXBASE_PLUGIN <vcg/math/eigen_matrixbase_addons.h>
#define EIGEN_MATRIX_PLUGIN <vcg/math/eigen_matrix_addons.h>

// forward declarations
namespace Eigen {
template<typename Derived1, typename Derived2, int Size> struct ei_lexi_comparison;
}

#include "base.h"
#include "../Eigen/LU"
#include "../Eigen/Geometry"
#include "../Eigen/Array"
#include "../Eigen/Core"

// add support for unsigned char and short int
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

// WARNING this is a default version provided so that Intersection() stuff can compile.
// Indeed, the compiler try to instanciate all versions of Intersection() leading to
// the instanciation of Eigen::Matrix<Face,...> !!!
template<typename T> struct NumTraits
{
	struct wrong_type
	{
		wrong_type() { assert(0 && "Eigen: you are using a wrong scalar type" ); }
	};

	typedef wrong_type Real;
	typedef wrong_type FloatingPoint;
	enum {
		IsComplex = 0,
		HasFloatingPoint = 0,
		ReadCost = 0,
		AddCost = 0,
		MulCost = 0
	};
};

// implementation of Lexicographic order comparison
// TODO should use meta unrollers
template<typename Derived1, typename Derived2> struct ei_lexi_comparison<Derived1,Derived2,2>
{
	inline static bool less(const Derived1& a, const Derived2& b) {
		return (a.coeff(1)!=b.coeff(1))?(a.coeff(1)< b.coeff(1)) : (a.coeff(0)<b.coeff(0));
	}

	inline static bool greater(const Derived1& a, const Derived2& b) {
		return (a.coeff(1)!=b.coeff(1))?(a.coeff(1)> b.coeff(1)) : (a.coeff(0)>b.coeff(0));
	}

	inline static bool lessEqual(const Derived1& a, const Derived2& b) {
		return (a.coeff(1)!=b.coeff(1))?(a.coeff(1)< b.coeff(1)) : (a.coeff(0)<=b.coeff(0));
	}

	inline static bool greaterEqual(const Derived1& a, const Derived2& b) {
		return (a.coeff(1)!=b.coeff(1))?(a.coeff(1)> b.coeff(1)) : (a.coeff(0)>=b.coeff(0));
	}
};

template<typename Derived1, typename Derived2> struct ei_lexi_comparison<Derived1,Derived2,3>
{
	inline static bool less(const Derived1& a, const Derived2& b) {
		return	(a.coeff(2)!=b.coeff(2))?(a.coeff(2)< b.coeff(2)):
				    (a.coeff(1)!=b.coeff(1))?(a.coeff(1)< b.coeff(1)) : (a.coeff(0)<b.coeff(0));
	}

	inline static bool greater(const Derived1& a, const Derived2& b) {
		return	(a.coeff(2)!=b.coeff(2))?(a.coeff(2)> b.coeff(2)):
				    (a.coeff(1)!=b.coeff(1))?(a.coeff(1)> b.coeff(1)) : (a.coeff(0)>b.coeff(0));
	}

	inline static bool lessEqual(const Derived1& a, const Derived2& b) {
		return	(a.coeff(2)!=b.coeff(2))?(a.coeff(2)< b.coeff(2)):
				    (a.coeff(1)!=b.coeff(1))?(a.coeff(1)< b.coeff(1)) : (a.coeff(0)<=b.coeff(0));
	}

	inline static bool greaterEqual(const Derived1& a, const Derived2& b) {
		return	(a.coeff(2)!=b.coeff(2))?(a.coeff(2)> b.coeff(2)):
				    (a.coeff(1)!=b.coeff(1))?(a.coeff(1)> b.coeff(1)) : (a.coeff(0)>=b.coeff(0));
	}
};

template<typename Derived1, typename Derived2> struct ei_lexi_comparison<Derived1,Derived2,4>
{
	inline static bool less(const Derived1& a, const Derived2& b) {
		return	(a.coeff(3)!=b.coeff(3))?(a.coeff(3)< b.coeff(3)) : (a.coeff(2)!=b.coeff(2))?(a.coeff(2)< b.coeff(2)):
				    (a.coeff(1)!=b.coeff(1))?(a.coeff(1)< b.coeff(1)) : (a.coeff(0)<b.coeff(0));
	}

	inline static bool greater(const Derived1& a, const Derived2& b) {
		return	(a.coeff(3)!=b.coeff(3))?(a.coeff(3)> b.coeff(3)) : (a.coeff(2)!=b.coeff(2))?(a.coeff(2)> b.coeff(2)):
				    (a.coeff(1)!=b.coeff(1))?(a.coeff(1)> b.coeff(1)) : (a.coeff(0)>b.coeff(0));
	}

	inline static bool lessEqual(const Derived1& a, const Derived2& b) {
		return	(a.coeff(3)!=b.coeff(3))?(a.coeff(3)< b.coeff(3)) : (a.coeff(2)!=b.coeff(2))?(a.coeff(2)< b.coeff(2)):
				    (a.coeff(1)!=b.coeff(1))?(a.coeff(1)< b.coeff(1)) : (a.coeff(0)<=b.coeff(0));
	}

	inline static bool greaterEqual(const Derived1& a, const Derived2& b) {
		return	(a.coeff(3)!=b.coeff(3))?(a.coeff(3)> b.coeff(3)) : (a.coeff(2)!=b.coeff(2))?(a.coeff(2)> b.coeff(2)):
				    (a.coeff(1)!=b.coeff(1))?(a.coeff(1)> b.coeff(1)) : (a.coeff(0)>=b.coeff(0));
	}
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

template<typename Derived1, typename Derived2>
typename Eigen::ei_traits<Derived1>::Scalar
AngleN(const Eigen::MatrixBase<Derived1>& p1, const Eigen::MatrixBase<Derived2> & p2)
{
	EIGEN_STATIC_ASSERT_VECTOR_ONLY(Derived1)
	EIGEN_STATIC_ASSERT_VECTOR_ONLY(Derived2)
	EIGEN_STATIC_ASSERT_FIXED_SIZE(Derived1)
	EIGEN_STATIC_ASSERT_FIXED_SIZE(Derived2)
	EIGEN_STATIC_ASSERT_SAME_VECTOR_SIZE(Derived1,Derived2)
	typedef typename Eigen::ei_traits<Derived1>::Scalar Scalar;

	Scalar t = (p1.dot(p2));
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

template<typename Derived>
inline const Eigen::CwiseUnaryOp<Eigen::ei_scalar_abs_op<typename Eigen::ei_traits<Derived>::Scalar>, Derived>
Abs(const Eigen::MatrixBase<Derived>& p)
{ return p.cwise().abs(); }

template<typename Derived>
inline const Eigen::CwiseBinaryOp<Eigen::ei_scalar_max_op<typename Eigen::ei_traits<Derived>::Scalar>,
																	Derived,
																	Eigen::NestByValue<typename Derived::ConstantReturnType> >
LowClampToZero(const Eigen::MatrixBase<Derived>& p)
{ return p.cwise().max(Derived::Zero().nestByValue()); }

}

#endif
