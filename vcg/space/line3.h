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

Revision 1.1  2004/03/08 16:15:48  tarini
first version (tarini)



****************************************************************************/



#ifndef __VCGLIB_LINE3
#define __VCGLIB_LINE3

#include <vcg/space/point3.h>

namespace vcg {

/** \addtogroup space */
/*@{*/
/** 
Templated class for 3D lines.
  This is the class for infinite lines in 3D space. A Line is stored just as two Point3:
	an origin and a direction (not necessarily normalized).
	@param LineScalarType (template parameter) Specifies the type of scalar used to represent coords.
	@param NORM: if on, the direction is always Normalized
*/
template <class LineScalarType, bool NORM=false>
class Line3
{
public: 

	/// The scalar type
	typedef LineScalarType ScalarType;

	/// The point type
	typedef Point3<LineScalarType> PointType;

	/// The point type
	typedef Line3<LineScalarType,NORM> LineType;

private:

	/// Origingin
	PointType _ori;

	/// Directionection (not necessarily normalized)
	PointType _dir;

public:

		/// Members to access the origin, direction
  inline const PointType &Origin() const { return _ori; } 
  inline const PointType &Direction() const { return _dir; } 
  inline PointType &Origin() { return _ori; } 
  inline PointType &Direction() {
		assert(!IsNormalized()); // Directionection can't be set for NORMALIZED Lines! Use SetDirection instead!
		return _dir; 
	} 
		/// The empty constructor
	Line3() {};
		/// The (origin, direction) constructor
	LineType(const PointType &ori, const PointType &dir) {SetOrigin(ori); SetDirection(dir);};
		/// Operator to compare two lines
	inline bool operator == ( LineType const & p ) const
	{	return _ori==p._ori && _dir==p._dir; }
		/// Operator to dispare two lines
	inline bool operator != ( LineType const & p ) const
	{	return _ori!=p._ori || _dir!=p._dir; }
		/// Projects a point on the line
	inline ScalarType Projection( const  PointType &p ) const
	{ if (NORM) return ScalarType((p-_ori)*_dir); 
		else      return ScalarType((p-_ori)*_dir/_dir.SquaredNorm()); 
	}
	inline bool IsNormalized() const {return NORM;};
		///set the origin
	inline void SetOrigin( const PointType & ori )
	{	_ori=ori; }
		///set the direction
	inline void SetDirection( const PointType & dir)
	{	_dir=dir; if (NORM) _dir.Normalize();  }
		///set both the origina and direction.
	inline void Set( const PointType & ori, const PointType & dir )
	{	SetOrigin(ori); SetDirection(dir); }
	  /// calculates the point of parameter t on the line.
	inline PointType P( const ScalarType t ) const
	{ return _ori + _dir * t; }
		/// normalizes direction field (returns a Normalized Line)
	Line3<ScalarType,true> &Normalize()
	{ if (!NORM) _dir.Normalize(); return *((Line3<ScalarType,true>*)this);}
		/// normalizes direction field (returns a Normalized Line) - static version
	static Line3<ScalarType,true> &Normalize(LineType &p)
	{ p.Normalize(); return *((Line3<ScalarType,true>*)(&p));}
	  /// importer for different line types
	template <class Q, bool K>
	inline void Import( const Line3<Q,K> & b )
	{ _ori.Import( b.Origin());	_dir.Import( b.Direction()); 
	  if ((NORM) && (!K)) _dir.Normalize();
	}
	PointType ClosestPoint(const PointType & p) const{
	return P(Projection(p));
	}
	  /// flips the line
	inline void Flip(){
		_dir=-_dir;
	};
}; // end class definition

typedef Line3<short>  Line3s;
typedef Line3<int>	  Line3i;
typedef Line3<float>  Line3f;
typedef Line3<double> Line3d;

typedef Line3<short ,true> Line3sN;
typedef Line3<int   ,true> Line3iN;
typedef Line3<float ,true> Line3fN;
typedef Line3<double,true> Line3dN;

	  /// returns closest point
template <class ScalarType, bool NORM> 
Point3<ScalarType> ClosestPoint( Line3<ScalarType,NORM> l, const Point3<ScalarType> & p) 
{
	return l.P(l.Projection(p)); 
}

/*@}*/

} // end namespace
#endif