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



#ifndef __VCGLIB_RAY3
#define __VCGLIB_RAY3

#include <vcg/space/point3.h>

namespace vcg {

/** \addtogroup space */
/*@{*/
/** 
Templated class for 3D ray.
  This is the class for a semi-line in 3D space. A ray is stored just as origin and direction (Point3).
	@param RayScalarType (template parameter) Specifies the type of scalar used to represent coords.
*/
template <class RayScalarType >
class Ray3
{
public: 

	/// The scalar type
	typedef RayScalarType ScalarType;

	/// The point type
	typedef Point3<RayScalarType> PointType;

	/// The point type
	typedef Ray3<RayScalarType> RayType;

private:

	/// Origin
	PointType _ori;

	/// Direction (not necessarily normalized)
	PointType _dir;

public:

		/// Members to access the origin, direction
  inline const PointType &Ori() const { return _ori; } 
  inline const PointType &Dir() const { return _dir; } 
  inline PointType &Ori() { return _ori; } 
  inline PointType &Dir() { return _dir; } 
		/// The empty constructor
	Ray3() {};
		/// The (a,b) constructor
	RayType(const PointType &a, const PointType &b) { _p0=a; _p1=b; };
		/// Operator to compare two rays
	inline bool operator == ( LineType const & p ) const
	{	return _ori==p._ori && _dir==p._dir; }
		/// Operator to dispare two rays
	inline bool operator != ( LineType const & p ) const
	{	return _ori!=p._ori || _dir!=p._dir; }
		/// Projects a point on the ray (returns distance)
	  /// projection exists iff result > 0
	inline ScalarType Projection( const  PointType &p ) const
	{ ScalarType l = dire.SquaredNorm();
		return ScalarType((p-_ori)*_dir/l); }
		///set up of the line.
	void Set( const PointType & ori, const PointType & dir )
	{	_ori = ori; _dir=dir }
	  /// calculates the point of parameter t on the ray.
	  /// if t>0, point is on the ray
	inline PointType P( const ScalarType t ) const
	{ return orig + dire * t; }
		/// normalizes direction field
	LineType &Normalize()
	{ _dir.Normalize(); return *this;}
	static LineType &Normalize(LineType &p)
	{ p.Normalize(); return p;}
	  /// importer for different Ray types
	template <class Q>
	inline void Import( const Ray3<Q> & b )
	{ _ori.Import( b._ori);	_dir.Import( b._dir);
	}

}; // end class definition



typedef Ray3<short>  Ray3s;
typedef Ray3<int>	   Ray3i;
typedef Ray3<float>  Ray3f;
typedef Ray3<double> Ray3d;


/*@}*/

} // end namespace
#endif