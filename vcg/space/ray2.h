/****************************************************************************
* VCGLib                                                            o o     *
* Visual and Computer Graphics Library                            o     o   *
*                                                                _   O  _   *
* Copyright(C) 2004-2016                                           \/)\/    *
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



#ifndef __VCGLIB_RAY2
#define __VCGLIB_RAY2

#include <vcg/space/point2.h>

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
class Ray2
{
public: 

	/// The scalar type
	typedef RayScalarType ScalarType;

	/// The point type
	typedef Point2<RayScalarType> PointType;

	/// The ray type
	typedef Ray2<RayScalarType,NORM> RayType;

private:

	/// Origin
	PointType _ori;

	/// Direction (not necessarily normalized, unless so specified by NORM)
	PointType _dir;

public:

//@{
	 /** @name Members to access the origin or direction
	   Direction() cannot be assigned directly.
		 Use SetDirection() or Set() instead.
	**/
		/// 
  inline const PointType &Origin() const { return _ori; } 
  inline PointType &Origin() { return _ori; } 
  inline const PointType &Direction() const { return _dir; } 
		/// sets the origin
	inline void SetOrigin( const PointType & ori )
	{	_ori=ori; }
		/// sets the direction
	inline void SetDirection( const PointType & dir)
	{	_dir=dir; if (NORM) _dir.Normalize();  }
		/// sets origin and direction.
	inline void Set( const PointType & ori, const PointType & dir )
	{	SetOrigin(ori); SetDirection(dir); }
//@}

//@{
	 /** @name Constructors 
	**/
 		/// The empty constructor
	Ray2() {};
		/// The (origin, direction) constructor
	Ray2(const PointType &ori, const PointType &dir) {SetOrigin(ori); SetDirection(dir);};
//@}

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
	  /// returns whether this type is normalized or not
	static bool IsNormalized() {return NORM;};
	  /// calculates the point of parameter t on the ray.
	inline PointType P( const ScalarType t ) const
	{ return _ori + _dir * t; }
		/// normalizes direction field (returns a Normalized Ray)
	inline Ray2<ScalarType,true> &Normalize()
	{ if (!NORM) _dir.Normalize(); return *((Ray2<ScalarType,true>*)this);}
		/// normalizes direction field (returns a Normalized Ray) - static version
	static Ray2<ScalarType,true> &Normalize(RayType &p)
	{ p.Normalize(); return *((Ray2<ScalarType,true>*)(&p));}
	  /// importer for different ray types (with any scalar type or normalization beaviour)
	template <class Q, bool K>
	inline void Import( const Ray2<Q,K> & b )
	{ _ori.Import( b.Origin() );	_dir.Import( b.Direction() ); 
	  if ((NORM) && (!K)) _dir.Normalize(); 
		//printf("(=)%c->%c ",(!NORM)?'N':'n', NORM?'N':'n');
	}
		/// constructs a new ray importing it from an existing one
	template <class Q, bool K>
	static RayType Construct( const Ray2<Q,K> & b )
	{ RayType res; res.Import(b);  return res;
	}	
	PointType ClosestPoint(const PointType & p) const{
	return P(Projection(p));
	}
	  /// flips the ray
	inline void Flip(){
		_dir=-_dir;
	};

//@{
	 /** @name Linearity for 3d rays 
   (operators +, -, *, /) so a ray can be set as a linear combination
	 of several rays. Note that the result of any operation returns 
	 a non-normalized ray; however, the command r0 = r1*a + r2*b is licit 
	 even if r0,r1,r2 are normalized rays, as the normalization will
	 take place within the final assignement operation. 
	**/
	inline Ray2<ScalarType,false> operator + ( RayType const & p) const
	{return Ray2<ScalarType,false> ( _ori+p.Origin(), _dir+p.Direction() );}
	inline Ray2<ScalarType,false> operator - ( RayType const & p) const
	{return Ray2<ScalarType,false> ( _ori-p.Origin(), _dir-p.Direction() );}
	inline Ray2<ScalarType,false> operator * ( const ScalarType s ) const
	{return Ray2<ScalarType,false> ( _ori*s, _dir*s );}
	inline Ray2<ScalarType,false> operator / ( const ScalarType s ) const
	{ScalarType s0=((ScalarType)1.0)/s; return RayType( _ori*s0, _dir*s0 );}
//@}


//@{
	 /** @name Automatic normalized to non-normalized
	 "Ray2dN r0 = r1" is equivalent to
	 "Ray2dN r0 = r1.Normalize()" if r1 is a Ray2d
	**/
		/// copy constructor that takes opposite beaviour
	Ray2(const Ray2<ScalarType,!NORM > &r) 
	{ Import(r); };
		/// assignment
	inline RayType & operator = ( Ray2<ScalarType,!NORM> const &r) 
	{ Import(r); return *this; };
//@}

}; // end class definition

typedef Ray2<short>  Ray2s;
typedef Ray2<int>    Ray2i;
typedef Ray2<float>  Ray2f;
typedef Ray2<double> Ray2d;

typedef Ray2<short ,true> Ray2sN;
typedef Ray2<int   ,true> Ray2iN;
typedef Ray2<float ,true> Ray2fN;
typedef Ray2<double,true> Ray2dN;

	  /// returns closest point
template <class ScalarType, bool NORM> 
Point2<ScalarType> ClosestPoint( Ray2<ScalarType,NORM> r, const Point2<ScalarType> & p)
{
	ScalarType t = r.Projection(p); 
	if (t<0) return r.Origin();
	return r.P(t); 
}

/*@}*/

} // end namespace
#endif
