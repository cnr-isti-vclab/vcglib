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
*/
template <class LineScalarType >
class Line3
{
public: 

	/// The scalar type
	typedef LineScalarType ScalarType;

	/// The point type
	typedef Point3<LineScalarType> PointType;

	/// The point type
	typedef Line3<LineScalarType> LineType;

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
	Line3() {};
		/// The (origin, direction) constructor
	LineType(const PointType &ori, const PointType &dir) { _ori=ori; _dir=dir; };
		/// Operator to compare two bounding box
	inline bool operator == ( LineType const & p ) const
	{	return _ori==p._ori && _dir==p._dir; }
		/// Operator to dispare two bounding box
	inline bool operator != ( LineType const & p ) const
	{	return _ori!=p._ori || _dir!=p._dir; }
		/// initializes the bounding box
	inline PointType Projection( const  PointType &p ) const
	{ ScalarType l = dire.SquaredNorm();
		return ScalarType((p-_ori)*_dir/l); }
		///set up of the line.
	void Set( const PointType & ori, const PointType & dir )
	{	_ori = ori; _dir=dir }
	  /// calculates the point of parameter t on the line.
	inline PointType P( const ScalarType t ) const
	{ return orig + dire * t; }
		/// normalizes direction field
	LineType &Normalize()
	{ _dir.Normalize(); return *this;}
	static LineType &Normalize(LineType &p)
	{ p.Normalize(); return p;}

}; // end class definition



typedef Line3<short>  Line3s;
typedef Line3<int>	  Line3i;
typedef Line3<float>  Line3f;
typedef Line3<double> Line3d;


/*@}*/

} // end namespace
#endif