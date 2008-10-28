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
#include "deprecated_point.h"
#else

#ifndef __VCGLIB_POINT
#define __VCGLIB_POINT

#include "../math/eigen.h"
#include <vcg/math/base.h>
#include <vcg/space/space.h>

namespace vcg{
template<class Scalar> class Point4;
}

namespace Eigen{
template<typename Scalar,int Size>
struct ei_traits<vcg::Point<Scalar,Size> > : ei_traits<Eigen::Matrix<Scalar,Size,1> > {};
}

namespace vcg {
namespace ndim{


/** \addtogroup space */
/*@{*/
/**
		The templated class for representing a point in R^N space.
		The class is templated over the ScalarType class that is used to represent coordinates.
		PointBase provides the interface and the common operators for points
		of any dimensionality.
	*/
template <int N, class S> class Point : public Eigen::Matrix<S,N,1>
{
	typedef Eigen::Matrix<T,N,1> _Base;
	using _Base::coeff;
	using _Base::coeffRef;
	using _Base::setZero;
	using _Base::data;
	using _Base::V;

public:

	_EIGEN_GENERIC_PUBLIC_INTERFACE(Point,_Base);
	
	typedef S ScalarType;
	typedef VoidType   ParamType;
	typedef Point<N,S> PointType;
	enum {Dimension = N};

	VCG_EIGEN_INHERIT_ASSIGNMENT_OPERATORS(Point)

//@{

	/** @name Standard Constructors and Initializers
		No casting operators have been introduced to avoid automatic unattended (and costly) conversion between different PointType types
		**/
	inline Point() : Base() {}
	template<typename OtherDerived>
 	inline Point(const Eigen::MatrixBase<OtherDerived>& other) : Base(other) {}
 	
  /// Padding function: give a default 0 value to all the elements that are not in the [0..2] range.
  /// Useful for managing in a consistent way object that could have point2 / point3 / point4
	inline S Ext( const int i ) const
	{
		if(i>=0 && i<=N) return data()[i];
		else             return 0;
	}

	/// importer for points with different scalar type and-or dimensionality
	template <int N2, class S2>
	inline void Import( const Point<N2,S2> & b )
	{
		data()[0] = ScalarType(b[0]);
		data()[1] = ScalarType(b[1]);
		if (N>2) { if (N2>2) data()[2] = ScalarType(b[2]); else data()[2] = 0;};
		if (N>3) { if (N2>3) data()[3] = ScalarType(b[3]); else data()[3] = 0;};
	}

	/// constructor for points with different scalar type and-or dimensionality
	template <int N2, class S2>
  static inline PointType Construct( const Point<N2,S2> & b )
  {
		PointType p; p.Import(b);
    return p;
  }

	  /// importer for homogeneous points
	template <class S2>
	inline void ImportHomo( const Point<N-1,S2> & b )
	{
		data()[0] = ScalarType(b[0]);
		data()[1] = ScalarType(b[1]);
		if (N>2) { data()[2] = ScalarType(data()[2]); };
		data()[N-1] = 1.0;
	}

		/// constructor for homogeneus point.
	template <int N2, class S2>
  static inline PointType Construct( const Point<N-1,S2> & b )
  {
		PointType p; p.ImportHomo(b);
    return p;
  }

//@}

//@{

  /** @name Data Access.
   access to data is done by overloading of [] or explicit naming of coords (x,y,z)**/

  inline const S &X() const { return data()[0]; }
	inline const S &Y() const { return data()[1]; }
	inline const S &Z() const { static_assert(N>2); return data()[2]; }
	 /// W is in any case the last coordinate.
	 /// (in a 2D point, W() == Y(). In a 3D point, W()==Z()
	 ///  in a 4D point, W() is a separate component)
	inline const S &W() const { return data()[N-1]; }
	inline S &X() { return data()[0]; }
	inline S &Y() { return data()[1]; }
	inline S &Z() { static_assert(N>2); return data()[2]; }
	inline S &W() { return data()[N-1]; }
//@}


//@{

  /** @name Dot products (cross product "%" is defined olny for 3D points)
  **/

	/// slower version, more stable (double precision only)
	inline S StableDot ( const PointType & p ) const;

//@}

//@{

  /** @name Norms
  **/

		/// Euclidean norm, static version
	template <class PT> static S Norm(const PT &p );
		/// Squared Euclidean norm, static version
	template <class PT> static S SquaredNorm(const PT &p );
		/// Normalization (division by norm), static version
	template <class PT> static PointType & Normalize(const PT &p);

//@}

		/// Signed area operator
	  /// a % b returns the signed area of the parallelogram inside a and b
	// inline S operator % ( PointType const & p ) const;

	 /// Convert to polar coordinates
	void ToPolar( S & ro, S & tetha, S & fi ) const
	{
		ro = Norm();
		tetha = (S)atan2( data()[1], data()[0] );
		fi    = (S)acos( data()[2]/ro );
	}

//@{

  /** @name Comparison Operators.
   Lexicographic order.
   **/

	inline bool operator == ( PointType const & p ) const;
	inline bool operator != ( PointType const & p ) const;
	inline bool operator <  ( PointType const & p ) const;
	inline bool operator >  ( PointType const & p ) const;
	inline bool operator <= ( PointType const & p ) const;
	inline bool operator >= ( PointType const & p ) const;
 //@}

//@{

  /** @name
	Glocal to Local and viceversa
	(provided for uniformity with other spatial classes. trivial for points)
   **/

	inline PointType LocalToGlobal(ParamType p) const { return *this; }
	inline ParamType GlobalToLocal(PointType p) const { ParamType p(); return p; }
//@}

}; // end class definition


// workaround the lack of template typedef (the next c++ standard will support them :) )

template <typename S>
struct Point2:public Point<2,S>{
	typedef Point<3,S> Base;
	inline Point2() : Base() {};
	inline Point2(const Point2& p):Base(p){};
	inline Point2(S a, S b):Base(a,b){};
	template<typename OtherDerived>
	inline Point2(const Eigen::MatrixBase<OtherDerived>& other) : Base(other) {}
};

template <typename S>
struct Point3:public Point<3,S> {
	typedef Point<3,S> Base;
	inline Point3() : Base() {};
	inline Point3(const Point3& p):Base(p){}
	inline Point3(S a, S b, S c):Base(a,b,c){};
	template<typename OtherDerived>
 	inline Point3(const Eigen::MatrixBase<OtherDerived>& other) : Base(other) {}
};


template <typename S>
struct Point4:public Point<4,S>{
	typedef Point<4,S> Base;
	inline Point4() : Base() {};
	inline Point4(const Point4& p):Base(p) {}
	inline Point4(S a, S b, S c, S d):Base(a,b,c,d){};
	template<typename OtherDerived>
 	inline Point4(const Eigen::MatrixBase<OtherDerived>& other) : Base(other) {}
};

typedef Point<2,short>  Point2s;
typedef Point<2,int>	  Point2i;
typedef Point<2,float>  Point2f;
typedef Point<2,double> Point2d;
typedef Point<2,short>  Vector2s;
typedef Point<2,int>	  Vector2i;
typedef Point<2,float>  Vector2f;
typedef Point<2,double> Vector2d;

typedef Point<3,short>  Point3s;
typedef Point<3,int>	  Point3i;
typedef Point<3,float>  Point3f;
typedef Point<3,double> Point3d;
typedef Point<3,short>  Vector3s;
typedef Point<3,int>	  Vector3i;
typedef Point<3,float>  Vector3f;
typedef Point<3,double> Vector3d;


typedef Point<4,short>  Point4s;
typedef Point<4,int>	  Point4i;
typedef Point<4,float>  Point4f;
typedef Point<4,double> Point4d;
typedef Point<4,short>  Vector4s;
typedef Point<4,int>	  Vector4i;
typedef Point<4,float>  Vector4f;
typedef Point<4,double> Vector4d;

/*@}*/

} // end namespace ndim
} // end namespace vcg
#endif

#endif
