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

#warning You are including deprecated math stuff

enum {Dimension = SizeAtCompileTime};
typedef vcg::VoidType   ParamType;
typedef Matrix          PointType;
using Base::V;

/// importer for points with different scalar type and-or dimensionality
// FIXME the Point3/Point4 specialization were only for same sizes ??
// while the Point version was generic like this one
template<typename OtherDerived>
inline void Import(const MatrixBase<OtherDerived>& b)
{
	EIGEN_STATIC_ASSERT_VECTOR_ONLY(Matrix);
	EIGEN_STATIC_ASSERT_FIXED_SIZE(Matrix);
	EIGEN_STATIC_ASSERT_VECTOR_ONLY(OtherDerived);
	EIGEN_STATIC_ASSERT_FIXED_SIZE(OtherDerived);
	enum { OtherSize = OtherDerived::SizeAtCompileTime };
	assert(SizeAtCompileTime<=4 && OtherSize<=4);
	data()[0] = Scalar(b[0]);
	if (SizeAtCompileTime>1) { if (OtherSize>1) data()[1] = Scalar(b[1]); else data()[1] = 0; }
	if (SizeAtCompileTime>2) { if (OtherSize>2) data()[2] = Scalar(b[2]); else data()[2] = 0; }
	if (SizeAtCompileTime>3) { if (OtherSize>3) data()[3] = Scalar(b[3]); else data()[3] = 0; }
}

/// importer for homogeneous points
template<typename OtherDerived>
inline void ImportHomo(const MatrixBase<OtherDerived>& b)
{
	EIGEN_STATIC_ASSERT_VECTOR_ONLY(Matrix);
	EIGEN_STATIC_ASSERT_FIXED_SIZE(Matrix);
	EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(OtherDerived,SizeAtCompileTime-1);
	
	this->template start<SizeAtCompileTime-1> = b;
	data()[SizeAtCompileTime-1] = Scalar(1.0);
}

/// constructor for points with different scalar type and-or dimensionality
template<typename OtherDerived>
static inline Matrix Construct(const MatrixBase<OtherDerived>& b)
{
	Matrix p; p.Import(b);
	return p;
}

	/// constructor for homogeneus point.
template<typename OtherDerived>
static inline Matrix ConstructHomo(const MatrixBase<OtherDerived>& b)
{
	Matrix p; p.ImportHomo(b);
	return p;
}

inline const Scalar &X() const { return data()[0]; }
inline const Scalar &Y() const { return data()[1]; }
inline const Scalar &Z() const { assert(SizeAtCompileTime>2); return data()[2]; }
inline Scalar &X() { return data()[0]; }
inline Scalar &Y() { return data()[1]; }
inline Scalar &Z() { assert(SizeAtCompileTime>2); return data()[2]; }

// note, W always returns the last entry
inline Scalar& W() { return data()[SizeAtCompileTime-1]; }
inline const Scalar& W() const { return data()[SizeAtCompileTime-1]; }

Scalar* V() { return data(); }
const Scalar* V() const { return data(); }

// overloaded to return a const reference
inline const Scalar& V( const int i ) const
{
	assert(i>=0 && i<SizeAtCompileTime);
	return data()[i];
}

//--------------------------------------------------------------------------------
// SPACE
//--------------------------------------------------------------------------------

/** Local to Glocal
  * (provided for uniformity with other spatial classes. trivial for points) */
inline Matrix LocalToGlobal(ParamType p) const { return *this; }
/** Glocal to Local
  * (provided for uniformity with other spatial classes. trivial for points) */
inline ParamType GlobalToLocal(PointType /*p*/) const { return ParamType(); }

/** 
	* Convert to polar coordinates from cartesian coordinates.
	*
	* Theta is the azimuth angle and ranges between [0, 360) degrees.
	* Phi is the elevation angle (not the polar angle) and ranges between [-90, 90] degrees.
	*
	* \note Note that instead of the classical polar angle, which ranges between
	*       0 and 180 degrees we opt for the elevation angle to obtain a more
	*       intuitive spherical coordinate system.
	*/
void ToPolar(Scalar &ro, Scalar &theta, Scalar &phi) const
{
	EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(Matrix,3);
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
	EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(Matrix,3);
	data()[0]= ro*ei_cos(theta)*ei_cos(phi);
	data()[1]= ro*ei_sin(phi);
	data()[2]= ro*ei_sin(theta)*ei_cos(phi);
}
