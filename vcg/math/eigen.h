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

#include "../Eigen/Array"
#include "../Eigen/Core"

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

#endif
