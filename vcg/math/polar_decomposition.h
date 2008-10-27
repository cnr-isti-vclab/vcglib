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


****************************************************************************/
#ifndef POLAR_DECOMPOSITION_VCG
#define POLAR_DECOMPOSITION_VCG


#include <vcg/math/matrix33.h>
#include <vcg/math/matrix44.h>
#include <vcg/math/lin_algebra.h>

namespace vcg{

	/** \addtogroup math */
	/* @{ */
	/// Extract the rotational part of a matrix by polar decomposition
	/// m = R*S. Polar decomposition is computed by taking r = m * sqrt(m^t*m)^{-1}
template <class S>
void RotationalPartByPolarDecomposition( const vcg::Matrix33<S> & m, vcg::Matrix33<S> &r ){

	vcg::Matrix33<S> tmp,s;

	r.SetZero();
	s.SetZero();

	tmp = m*m.transpose();

	Matrix33<S> res;
	Point3<S> e;

	bool ss = SingularValueDecomposition<vcg::Matrix33<S> >(tmp,&e[0],res);

	e[0]=math::Sqrt(e[0]);
	e[1]=math::Sqrt(e[1]);
	e[2]=math::Sqrt(e[2]);
	#ifdef VCG_USE_EIGEN
	tmp = tmp*e.asDiagonal()*res.transpose();
	#else
	tmp = tmp*Matrix33Diag<S>(e)*res.transpose();
	#endif

	bool s1 = SingularValueDecomposition<vcg::Matrix33<S> >(tmp,&e[0],res.transpose());
	e[0]=1/e[0];
	e[1]=1/e[1];
	e[2]=1/e[2];

	#ifdef VCG_USE_EIGEN
	tmp = res*e.asDiagonal()*tmp.transpose();
	#else
	tmp = res*Matrix33Diag<S>(e)*tmp.transpose();
	#endif

	r = m*tmp;
}

/*! @} */
};
#endif
