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
Revision 1.1  2004/03/08 19:46:47  tarini
First Version (tarini)



****************************************************************************/



#ifndef __VCGLIB_RAY3
#define __VCGLIB_RAY3

#include <vcg/space/point3.h>

namespace vcg {

/** \addtogroup space */
/*@{*/
/** 
Templated class for 3D rays.
  This is the class for infinite rays in 3D space. A Ray is stored just as two Point3:
	an origin and a direction (not necessarily normalized).
	@param RayScalarType (template parameter) Specifies the type of scalar used to represent coords.
	@param NORM: if on, the direction is always Normalized
*/
template <class RayScalarType, bool NORM=false>
class Ray3
{
public: 

	/// The scalar type
	typedef RayScalarType ScalarType;

	/// The point type
	typedef Point3<RayScalarType> PointType;

	/// The point type
	typedef Ray3<RayScalarType,NORM> RayType;

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
  inline PointType &Dir() { 
		assert(NORM); // Direction can't be set for NORMALIZED Rays! Use SetDir instead!
		return _dir; 
	} 
		/// The empty constructor
	Ray3() {};
		/// The (origin, direction) constructor
	RayType(const PointType &ori, const PointType &dir) {SetOri(ori); SetDir(dir);};
		/// Operator to compare two rays
	inline bool operator == ( RayType const & p ) const
	{	return _ori==p._ori && _dir==p._dir; }
		/// Operator to dispare two rays
	inline bool operator != ( RayType const & p ) const
	{	return _ori!=p._ori || _dir!=p._dir; }
		/// Projects a point on the ray
	inline ScalarType Projection( const  PointType &p ) const
	{ if (NORM) return ScalarType((p-_ori)*_dir); 
		else      return ScalarType((p-_ori)*_dir/_dir.SquaredNorm()); 
	}
	inline bool IsNorm() const {return NORM;};
		///set the origin
	inline void SetOri( const PointType & ori )
	{	_ori=ori; }
		///set the direction
	inline void SetDir( const PointType & dir)
	{	_dir=dir; if (NORM) _dir.Normalize();  }
		///set both the origina and direction.
	inline void Set( const PointType & ori, const PointType & dir )
	{	SetOri(ori); SetDir(dir); }
	  /// calculates the point of parameter t on the ray.
	inline PointType P( const ScalarType t ) const
	{ return _ori + _dir * t; }
		/// normalizes direction field (returns a Normalized Ray)
	inline Ray3<ScalarType,true> &Normalize()
	{ if (!NORM) _dir.Normalize(); return *((Ray3<ScalarType,true>*)this);}
		/// normalizes direction field (returns a Normalized Ray) - static version
	static Ray3<ScalarType,true> &Normalize(RayType &p)
	{ p.Normalize(); return *((Ray3<ScalarType,true>*)(&p));}
	  /// importer for different ray types
	template <class Q, bool K>
	inline void Import( const Ray3<Q,K> & b )
	{ _ori.Import( b.Ori() );	_dir.Import( b.Dir() ); 
	  if ((NORM) && (!K)) _dir.Normalize();
	}
	PointType ClosestPoint(const PointType & p) const{
	return P(Projection(p));
	}
	  /// flips the ray
	inline void Flip(){
		_dir=-_dir;
	};
}; // end class definition

typedef Ray3<short>  Ray3s;
typedef Ray3<int>    Ray3i;
typedef Ray3<float>  Ray3f;
typedef Ray3<double> Ray3d;

typedef Ray3<short ,true> Ray3sN;
typedef Ray3<int   ,true> Ray3iN;
typedef Ray3<float ,true> Ray3fN;
typedef Ray3<double,true> Ray3dN;

	  /// returns closest point
template <class ScalarType, bool NORM> 
Point3<ScalarType> ClosestPoint( Ray3<ScalarType,NORM> r, const Point3<ScalarType> & p) 
{
	ScalarType t = r.Projection(p); 
	if (t<0) t=0;
	return r.P(t); 
}

/*@}*/

} // end namespace
#endif