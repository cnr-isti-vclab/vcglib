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
Revision 1.1  2004/02/10 01:11:28  cignoni
Edited Comments and GPL license

****************************************************************************/

#ifndef __VCGLIB_POINT4
#define __VCGLIB_POINT4

#include <vcg/space/point3.h>

namespace vcg {



template <class T> class Point4
{
protected:
	T _v[4];

public:
	typedef T scalar;


	inline Point4 () { }
	inline Point4 ( const T nx, const T ny, const T nz , const T nw )
	{
		_v[0] = nx; _v[1] = ny; _v[2] = nz; _v[3] = nw;
	}
	inline Point4 ( const T  p[4] )
	{   
		_v[0] = p[0]; _v[1]= p[1]; _v[2] = p[2]; _v[3]= p[3];
	}
	inline Point4 ( const Point4 & p )
	{   
		_v[0]= p._v[0]; _v[1]= p._v[1]; _v[2]= p._v[2]; _v[3]= p._v[3];
	}
	inline Point4 ( const Point3<T> & p )
	{
		_v[0] = p.V(0);
		_v[1] = p.V(1);
		_v[2] = p.V(2);
		_v[3] = 1.0;
	}
	inline Point4 & operator = ( const Point4 & p )
	{
		_v[0]= p._v[0]; _v[1]= p._v[1]; _v[2]= p._v[2]; _v[3]= p._v[3];
		return *this;
	}
	inline T &x() {return _v[0];}
	inline T &y() {return _v[1];}
	inline T &z() {return _v[2];}
	inline T &w() {return _v[3];}
	inline const T & operator [] ( const int i ) const
	{
		assert(i>=0 && i<4);
		return _v[i];
	}
	inline T & operator [] ( const int i )
	{
		assert(i>=0 && i<4);
		return _v[i];
	}
	inline T const * V() const
	{
		return _v;
	}
	inline const T & V ( const int i ) const
	{
		assert(i>=0 && i<4);
		return _v[i];
	}
	inline T & V ( const int i )
	{
		assert(i>=0 && i<4);
		return _v[i];
	}
	
	inline Point4 operator + ( const Point4 & p) const
	{ 
		return Point4( _v[0]+p._v[0], _v[1]+p._v[1], _v[2]+p._v[2], _v[3]+p._v[3] );
	}
	inline Point4 operator - ( const Point4 & p) const
	{
		return Point4( _v[0]-p._v[0], _v[1]-p._v[1], _v[2]-p._v[2], _v[3]-p._v[3] );
	}
	inline Point4 operator * ( const T s ) const
	{
		return Point4( _v[0]*s, _v[1]*s, _v[2]*s, _v[3]*s );
	}
	inline Point4 operator / ( const T s ) const
	{
		return Point4( _v[0]/s, _v[1]/s, _v[2]/s, _v[3]/s );
	}
	inline T operator * ( const Point4 & p ) const
	{
		return _v[0]*p._v[0] + _v[1]*p._v[1] + _v[2]*p._v[2] + _v[3]*p._v[3];
	} 
	inline Point4 & operator += ( const Point4 & p)
	{
		_v[0] += p._v[0]; _v[1] += p._v[1]; _v[2] += p._v[2]; _v[3] += p._v[3];
		return *this;
	}
	inline Point4 & operator -= ( const Point4 & p )
	{
		_v[0] -= p._v[0]; _v[1] -= p._v[1]; _v[2] -= p._v[2]; _v[3] -= p._v[3];
		return *this;
	}
	inline Point4 & operator *= ( const T s )
	{
		_v[0] *= s; _v[1] *= s; _v[2] *= s; _v[3] *= s;
		return *this;
	}
	inline Point4 & operator /= ( const T s )
	{
		_v[0] /= s; _v[1] /= s; _v[2] /= s; _v[3] /= s;
		return *this;
	}
	inline T Norm() const
	{
		return Sqrt( _v[0]*_v[0] + _v[1]*_v[1] + _v[2]*_v[2] + _v[3]*_v[3] );
	}
	inline T SquaredNorm() const
	{
		return _v[0]*_v[0] + _v[1]*_v[1] + _v[2]*_v[2] + _v[3]*_v[3];
	}

  inline Point4 & Normalize()
	{
		T n = Sqrt(_v[0]*_v[0] + _v[1]*_v[1] + _v[2]*_v[2] + _v[3]*_v[3] );
		if(n>0.0) {	_v[0] /= n;	_v[1] /= n;	_v[2] /= n; _v[3] /= n; }
		return *this;
	}
	inline Point4 operator - () const
	{
		return Point4( -_v[0], -_v[1], -_v[2], -_v[3] );
	}
	inline bool operator == (  const Point4& p ) const
	{
		return _v[0]==p._v[0] && _v[1]==p._v[1] && _v[2]==p._v[2] && _v[3]==p._v[3];
	} 
	inline bool operator != ( const Point4 & p ) const
	{
		return _v[0]!=p._v[0] || _v[1]!=p._v[1] || _v[2]!=p._v[2] || _v[3]!=p._v[3];
	}
	inline bool operator <  ( Point4 const & p ) const
	{
		return	(_v[3]!=p._v[3])?(_v[3]<p._v[3]):
				(_v[2]!=p._v[2])?(_v[2]<p._v[2]):
				(_v[1]!=p._v[1])?(_v[1]<p._v[1]):
				(_v[0]<p._v[0]);
	}
	inline bool operator >  ( const Point4 & p ) const
	{
		return	(_v[3]!=p._v[3])?(_v[3]>p._v[3]):
				(_v[2]!=p._v[2])?(_v[2]>p._v[2]):
				(_v[1]!=p._v[1])?(_v[1]>p._v[1]):
				(_v[0]>p._v[0]);
	}
	inline bool operator <= ( const Point4 & p ) const
	{
		return	(_v[3]!=p._v[3])?(_v[3]< p._v[3]):
				(_v[2]!=p._v[2])?(_v[2]< p._v[2]):
				(_v[1]!=p._v[1])?(_v[1]< p._v[1]):
				(_v[0]<=p._v[0]);
	}
	inline bool operator >= ( const Point4 & p ) const
	{
		return	(_v[3]!=p._v[3])?(_v[3]> p._v[3]):
				(_v[2]!=p._v[2])?(_v[2]> p._v[2]):
				(_v[1]!=p._v[1])?(_v[1]> p._v[1]):
				(_v[0]>=p._v[0]);
	}
		/// Questa funzione estende il vettore ad un qualsiasi numero di dimensioni
		/// paddando gli elementi estesi con zeri
	inline T Ext( const int i ) const
	{
		if(i>=0 && i<=3) return _v[i];
		else             return 0;
	}

	T stable_dot ( const Point4<T> & p ) const
	{
		T k[4];

		k[0] = _v[0]*p._v[0];
		k[1] = _v[1]*p._v[1];
		k[2] = _v[2]*p._v[2];
		k[3] = _v[3]*p._v[3];
    sort(k+0,k+4, math::MagnitudoComparer<T>() );
		T q = k[0];
		q += k[1];
		q += k[2];
		q += k[3];
		return q;
	}  

	template <class Q>
	inline void Import( const Point4<Q> & b )
	{
		_v[0] = T(b[0]);
		_v[1] = T(b[1]);
		_v[2] = T(b[2]);
		_v[3] = T(b[3]);
	}

}; // end class definition

#ifdef __VCG_USE_P4_INTRINSIC__
#include <vcg/p4/point4p4.h>
#endif

template <class T>
T Angle( const Point4<T>& p1, const Point4<T>  & p2 )
{
	T w = p1.Norm()*p2.Norm();
	if(w==0) return -1;
	T t = (p1*p2)/w;
	if(t>1) t=1;
    return T( acos(t) );
}



template <class T>
inline T Norm( const Point4<T> & p )
{
	return p.Norm();
}

template <class T>
inline T SquaredNorm( const Point4<T> & p )
{
    return p.SquaredNorm();
}

/* Deprecato
template <class T>
inline Point4<T> & Normalize( Point4<T> & p ){
		T n = Sqrt( p._v[0]*p._v[0] + p._v[1]*p._v[1] + p._v[2]*p._v[2] + p._v[3]*p._v[3] );
    if(n>0.0) p/=n;
    return p;
}
*/

template <class T>
inline T Distance( const Point4<T> & p1, const Point4<T> & p2 )
{
    return Norm(p1-p2);
}

template <class T>
inline T SquaredDistance( const Point4<T> & p1, const Point4<T> & p2 )
{
    return SquaredNorm(p1-p2);
}


typedef Point4<short>  Point4s;
typedef Point4<int>	   Point4i;
typedef Point4<float>  Point4f;
typedef Point4<double> Point4d;


} // end namespace
#endif
