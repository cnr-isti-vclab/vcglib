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
Revision 1.2  2004/03/08 19:46:47  tarini
First Version (tarini)

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
	Line3() {};
		/// The (origin, direction) constructor
	LineType(const PointType &ori, const PointType &dir) {SetOri(ori); SetDir(dir);};
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
	bool IsNorm() const {return NORM;};
		///set the origin
	void SetOri( const PointType & ori )
	{	_ori=ori; }
		///set the direction
	void SetDir( const PointType & dir)
	{	_dir=dir; if (NORM) _dir.Normalize();  }
		///set both the origina and direction.
	void Set( const PointType & ori, const PointType & dir )
	{	SetOri(ori); SetDir(dir); }
	  /// calculates the point of parameter t on the line.
	inline PointType P( const ScalarType t ) const
	{ return orig + dire * t; }
		/// normalizes direction field
	Line3<ScalarType,true> &Normalize()
	{ if (!NORM) _dir.Normalize(); return *((Line3<ScalarType,true>*)this);}
	static Line3<ScalarType,true> &Normalize(LineType &p)
	{ p.Normalize(); return *((Line3<ScalarType,true>*)(&p));}
	  /// importer for different line types
	template <class Q, bool K>
	inline void Import( const Line3<Q,K> & b )
	{ _ori.Import( b._ori);	_dir.Import( b._dir); 
	  if ((NORM) && (!K)) _dir.Normalize();
	}
}; // end class definition



typedef Line3<short>  Line3s;
typedef Line3<int>	  Line3i;
typedef Line3<float>  Line3f;
typedef Line3<double> Line3d;

typedef Line3<short ,true> Line3sN;
typedef Line3<int   ,true> Line3iN;
typedef Line3<float ,true> Line3fN;
typedef Line3<double,true> Line3dN;


/*@}*/

} // end namespace
#endif