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

#ifndef __VCGLIB_POINT2
#define __VCGLIB_POINT2

//#include <limits>
#include <assert.h>
#include <vcg/math/base.h>

namespace vcg {

template <class FLTYPE> class Point2
{
protected:

	FLTYPE _v[2];

	typedef FLTYPE scalar;

	inline const FLTYPE &X() const {return v[0];} 
	inline const FLTYPE &Y() const {return v[1];}
	inline FLTYPE &X() {return v[0];}
	inline FLTYPE &Y() {return v[1];}
	inline const FLTYPE & operator [] ( const int i ) const
	{
		assert(i>=0 && i<2);
		return v[i];
	}
	inline FLTYPE & operator [] ( const int i )
	{
		assert(i>=0 && i<2);
		return v[i];
	}
	

	inline Point2 () { }
	inline Point2 ( const FLTYPE nx, const FLTYPE ny )
	{
			v[0] = nx; v[1] = ny;
	}
	inline Point2 ( Point2 const & p)
	{   
			v[0]= p.v[0];    v[1]= p.v[1];
	}
	inline Point2 & operator =( Point2 const & p)
	{
			v[0]= p.v[0]; v[1]= p.v[1];
			return *this;
	}

	inline void Zero()
	{
		v[0] = 0;
		v[1] = 0;
	}

	inline Point2 operator + ( Point2 const & p) const
	{ 
			return Point2<FLTYPE>( v[0]+p.v[0], v[1]+p.v[1] );
	}
	inline Point2 operator - ( Point2 const & p) const
	{
			return Point2<FLTYPE>( v[0]-p.v[0], v[1]-p.v[1] );
	}
	inline Point2 operator * ( const FLTYPE s ) const
	{
			return Point2<FLTYPE>( v[0] * s, v[1] * s );
	}
	inline Point2 operator / ( const FLTYPE s ) const
	{
			return Point2<FLTYPE>( v[0] / s, v[1] / s );
	}
	inline FLTYPE operator * ( Point2 const & p ) const
	{
			return ( v[0]*p.v[0] + v[1]*p.v[1] );
	}

	inline FLTYPE operator ^ ( Point2 const & p ) const
	{
			return v[1]*p.v[0] - v[0]*p.v[1];
	} 

	inline Point2 & operator += ( Point2 const & p)
	{
			v[0] += p.v[0];    v[1] += p.v[1];
			return *this;
	}
	inline Point2 & operator -= ( Point2 const & p)
	{
			v[0] -= p.v[0];    v[1] -= p.v[1];
			return *this;
	}
	inline Point2 & operator *= ( const FLTYPE s )
	{
			v[0] *= s;    v[1] *= s;
			return *this;
	}
	inline Point2 & operator /= ( const FLTYPE s )
	{
			v[0] /= s;    v[1] /= s;
			return *this;
	}
	inline FLTYPE Norm( void ) const
	{
			return Sqrt( v[0]*v[0] + v[1]*v[1] );
	}
	inline FLTYPE SquaredNorm( void ) const
	{
			return ( v[0]*v[0] + v[1]*v[1] );
	}
	inline Point2 & Scale( const FLTYPE sx, const FLTYPE sy );

	inline Point2 & Normalize( void )
	{
			FLTYPE n = Sqrt(v[0]*v[0] + v[1]*v[1]);
			if(n>0.0) {	v[0] /= n;	v[1] /= n; }
			return *this;
	}
	inline bool operator == ( Point2 const & p ) const
	{
			return (v[0]==p.v[0] && v[1]==p.v[1]);
	} 
	inline bool operator != ( Point2 const & p ) const
	{
			return ( (v[0]!=p.v[0]) || (v[1]!=p.v[1]) );
	}
	inline bool operator <  ( Point2 const & p ) const
	{
			return	(v[1]!=p.v[1])?(v[1]<p.v[1]):
							(v[0]<p.v[0]);
	}
	inline bool operator >  ( Point2 const & p ) const
	{
			return	(v[1]!=p.v[1])?(v[1]>p.v[1]):
							(v[0]>p.v[0]);
	}

	inline bool operator <= ( Point2 const & p ) const
	{
			return	(v[1]!=p.v[1])?(v[1]< p.v[1]):
							(v[0]<=p.v[0]);
	}

	inline bool operator >= ( Point2 const & p ) const
	{
			return	(v[1]!=p.v[1])?(v[1]> p.v[1]):
							(v[0]>=p.v[0]);
	}
	inline FLTYPE Distance( Point2 const & p ) const
	{
			return Norm(*this-p);
	}

	inline FLTYPE SquaredDistance( Point2 const & p ) const
	{
			return Norm2(*this-p);
	}	

	inline Point2 & Cartesian2Polar()
	{
		FLTYPE t = (FLTYPE)atan2(v[1],v[0]);
		v[0] = Sqrt(v[0]*v[0]+v[1]*v[1]);
		v[1] = t;
		return *this;
	}

	inline Point2 & Polar2Cartesian()
	{
		FLTYPE l = v[0];
		v[0] = (FLTYPE)(l*cos(v[1]));
		v[1] = (FLTYPE)(l*sin(v[1]));
		return *this;
	}


	inline Point2 & rotate( const FLTYPE a )
	{
		FLTYPE t = v[0];
		FLTYPE s = sin(a);
		FLTYPE c = cos(a);

		v[0] = v[0]*c - v[1]*s;
		v[1] =   t *s + v[1]*c;

		return *this;
	}

	/// Questa funzione estende il vettore ad un qualsiasi numero di dimensioni
	/// paddando gli elementi estesi con zeri
	inline FLTYPE Ext( const int i ) const
	{
		if(i>=0 && i<2) return v[i];
		else            return 0;
	}


}; // end class definition


template <class FLTYPE>
inline FLTYPE Angle( Point2<FLTYPE> const & p1, Point2<FLTYPE> const & p2 )
{
	return atan2(p2[1],p2[0]) - atan2(p1[1],p1[0]);
}

template <class FLTYPE>
inline Point2<FLTYPE> operator - ( Point2<FLTYPE> const & p ){
    return Point2<FLTYPE>( -p.v[0], -p.v[1] );
}

template <class FLTYPE>
inline Point2<FLTYPE> operator * ( const FLTYPE s, Point2<FLTYPE> const & p ){
    return Point2<FLTYPE>( p.v[0] * s, p.v[1] * s  );
}

template <class FLTYPE>
inline FLTYPE Norm( Point2<FLTYPE> const & p ){
		return Sqrt( p.v[0]*p.v[0] + p.v[1]*p.v[1] );
}

template <class FLTYPE>
inline FLTYPE Norm2( Point2<FLTYPE> const & p ){
    return ( p.v[0]*p.v[0] + p.v[1]*p.v[1] );
}

template <class FLTYPE>
inline Point2<FLTYPE> & Normalize( Point2<FLTYPE> & p ){
		FLTYPE n = Sqrt( p.v[0]*p.v[0] + p.v[1]*p.v[1] );
    if(n>0.0) p/=n;
    return p;
}

template <class FLTYPE>
inline FLTYPE Distance( Point2<FLTYPE> const & p1,Point2<FLTYPE> const & p2 ){
    return Norm(p1-p2);
}

template <class FLTYPE>
inline FLTYPE SquaredDistance( Point2<FLTYPE> const & p1,Point2<FLTYPE> const & p2 ){
    return Norm2(p1-p2);
}

typedef Point2<short>  Point2s;
typedef Point2<int>	   Point2i;
typedef Point2<float>  Point2f;
typedef Point2<double> Point2d;


} // end namespace
#endif
