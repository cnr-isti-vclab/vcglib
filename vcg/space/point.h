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
Revision 1.14  2004/03/05 17:55:01  tarini
errorino: upper case in Zero()

Revision 1.13  2004/03/03 14:22:48  cignoni
Yet against cr lf mismatch

Revision 1.12  2004/02/23 23:42:26  cignoni
Translated comments, removed unusued stuff. corrected linefeed/cr

Revision 1.11  2004/02/19 16:12:28  cignoni
cr lf mismatch 2

Revision 1.10  2004/02/19 16:06:24  cignoni
cr lf mismatch

Revision 1.8  2004/02/19 15:13:40  cignoni
corrected sqrt and added doxygen groups

Revision 1.7  2004/02/17 02:08:47  cignoni
Di prova...

Revision 1.6  2004/02/15 23:35:47  cignoni
Cambiato nome type template in accordo alla styleguide

Revision 1.5  2004/02/10 01:07:15  cignoni
Edited Comments and GPL license

Revision 1.4  2004/02/09 13:48:02  cignoni
Edited doxygen comments
****************************************************************************/

#ifndef __VCGLIB_POINT
#define __VCGLIB_POINT

#include <assert.h>
#include <vcg/math/base.h>
#include <vcg/space/space.h>

namespace vcg {
/** \addtogroup space */
/*@{*/
    /**
        The templated class for representing a point in 3D space.
        The class is templated over the ScalarType class that is used to represent coordinates. All the usual
        operator overloading (* + - ...) is present. 
     */

template <int N, class S> 
class Point : public Space<N,S>, Linear<Point>
{
public:
	typedef S        ScalarType;
	typedef VoidType ParamType;
	typedef Point    PointType;
	enum {Dimension=N};

protected:
  /// The only data member. Hidden to user.
	S _v[N];

public:

//@{

  /** @name Standard Constructors and Initializers 
   No casting operators have been introduced to avoid automatic unattended (and costly) conversion between different point types
   **/

  inline Point () { }
	inline Point ( const S nx, const S ny, const S nz, const S nw )
	{
		static_assert(N==4);
		_v[0] = nx;
		_v[1] = ny;
		_v[2] = nz;
		_v[3] = nw;
	}
	inline Point ( const S nx, const S ny, const S nz)
	{
		static_assert(N==3);
		_v[0] = nx;
		_v[1] = ny;
		_v[2] = nz;
	}
	inline Point ( const S nx, const S ny)
	{
		static_assert(N==2);
		_v[0] = nx;
		_v[1] = ny;
	}
	inline Point ( const S nv[N] )
	{
		_v[0] = nv[0];
		_v[1] = nv[1];
		if (N>2) _v[2] = nv[2];
		if (N>3) _v[3] = nv[3];
	}
  
  /// Padding function: give a default 0 value to all the elements that are not in the [0..2] range. 
  /// Useful for managing in a consistent way object that could have point2 / point3 / point4
	inline S Ext( const int i ) const
	{
		if(i>=0 && i<=N) return _v[i];
		else             return 0;
	}

	template <int N2, class S2>
	inline void Import( const Point<N2,S2> & b )
	{
		_v[0] = ScalarType(b[0]);
		_v[1] = ScalarType(b[1]);
		if (N>2) { if (N2>2) _v[2] = ScalarType(b[2]); else _v[2] = 0};
		if (N>3) { if (N2>3) _v[3] = ScalarType(b[3]); else _v[3] = 0};
	}

  static inline Point Construct( const PointType & b )
  {
		PointType p; p.Import(b);
    return p;
  }

//@}

//@{

  /** @name Data Access. 
   access to data is done by overloading of [] or explicit naming of coords (x,y,z)**/

	inline S & operator [] ( const int i )
	{
		assert(i>=0 && i<N);
		return _v[i];
	}
	inline const S & operator [] ( const int i ) const
	{
		assert(i>=0 && i<3);
		return _v[i];
	}
  inline const S &X() const { return _v[0]; } 
	inline const S &Y() const { return _v[1]; }
	inline const S &Z() const { static_assert(N>2); return _v[2]; }
	inline const S &W() const { static_assert(N>3); return _v[3]; }
	inline S &X() { return _v[0]; }
	inline S &Y() { return _v[1]; }
	inline S &Z() { static_assert(N>2); return _v[2]; }
	inline S &W() { static_assert(N>3); return _v[3]; }
	inline const S * V() const
	{
		return _v;
	}
	inline S & V( const int i )
	{
		assert(i>=0 && i<N);
		return _v[i];
	}
	inline const S & V( const int i ) const
	{
		assert(i>=0 && i<N);
		return _v[i];
	}
//@}
//@{

  /** @name Linearity for points 
  **/

	/// sets a point to Zero
	inline void Zero()
	{
		_v[0] = 0;
		_v[1] = 0;
		if (N>2) _v[2] = 0;
		if (N>3) _v[3] = 0;
	}
	inline Point operator + ( Point const & p) const
	{
		if (N==2) return Point( _v[0]+p._v[0], _v[1]+p._v[1] );
		if (N==3) return Point( _v[0]+p._v[0], _v[1]+p._v[1], _v[2]+p._v[2] );
		if (N==4) return Point( _v[0]+p._v[0], _v[1]+p._v[1], _v[2]+p._v[2], _v[3]+p._v[3] );
	}
	inline Point operator - ( Point const & p) const
	{
		if (N==2) return Point( _v[0]-p._v[0], _v[1]-p._v[1] );
		if (N==3) return Point( _v[0]-p._v[0], _v[1]-p._v[1], _v[2]-p._v[2] );
		if (N==4) return Point( _v[0]-p._v[0], _v[1]-p._v[1], _v[2]-p._v[2], _v[3]-p._v[3] );
	}
	inline Point operator * ( const S s ) const
	{
		if (N==2) return Point( _v[0]*s, _v[1]*s );
		if (N==3) return Point( _v[0]*s, _v[1]*s, _v[2]*s );
		if (N==4) return Point( _v[0]*s, _v[1]*s, _v[2]*s, _v[3]*s );
	}
	inline Point operator / ( const S s ) const
	{
		if (N==2) return Point( _v[0]/s, _v[1]/s );
		if (N==3) return Point( _v[0]/s, _v[1]/s, _v[2]/s );
		if (N==4) return Point( _v[0]/s, _v[1]/s, _v[2]/s, _v[3]/s );
	}
	inline Point & operator += ( Point const & p)
	{
		_v[0] += p._v[0];
		_v[1] += p._v[1];
		if (N>2) _v[2] += p._v[2];
		if (N>3) _v[3] += p._v[3];
		return *this;
	}
	inline Point & operator -= ( Point const & p)
	{
		_v[0] -= p._v[0];
		_v[1] -= p._v[1];
		if (N>2) _v[2] -= p._v[2];
		if (N>3) _v[3] -= p._v[3];
		return *this;
	}
	inline Point & operator *= ( const S s )
	{
		_v[0] *= s;
		_v[1] *= s;
		if (N>2) _v[2] *= s;
		if (N>3) _v[3] *= s;
		return *this;
	}
	inline Point & operator /= ( const S s )
	{
		_v[0] /= s;
		_v[1] /= s;
		if (N>2) _v[2] /= s;
		if (N>3) _v[3] /= s;
		return *this;
	}
	inline Point operator - () const
	{
		if (N==2) return Point ( -_v[0], -_v[1] );
		if (N==3) return Point ( -_v[0], -_v[1], -_v[2] );
		if (N==4) return Point ( -_v[0], -_v[1], -_v[2] , -_v[3] );
	}
//@}
//@{

  /** @name Dot products
  **/
		/// Dot product
	inline S operator * ( Point const & p ) const
	{
		if (N==2) return ( _v[0]*p._v[0] + _v[1]*p._v[1]  );
		if (N==3) return ( _v[0]*p._v[0] + _v[1]*p._v[1] + _v[2]*p._v[2] );
		if (N==4) return ( _v[0]*p._v[0] + _v[1]*p._v[1] + _v[2]*p._v[2] + _v[2]*p._v[2] );
	};
	/// slower version, more stable (double precision only)
	inline S StableDot ( const Point & p ) const
	{
		if (N==2) return _v[0]*p._v[0] + _v[1]*p._v[1];
		if (N==4) {

			S k0=_v[0]*p._v[0],	k1=_v[1]*p._v[1], k2=_v[2]*p._v[2], k3=_v[3]*p._v[3];
			int exp0,exp1,exp2,exp3;

			frexp( double(k0), &exp0 );frexp( double(k1), &exp1 );
			frexp( double(k2), &exp2 );frexp( double(k3), &exp3 );

			if (exp0>exp1) { math::Swap(k0,k1); math::Swap(exp0,exp1); }
			if (exp2>exp3) { math::Swap(k2,k3); math::Swap(exp2,exp3); }
			if (exp0>exp2) { math::Swap(k0,k2); math::Swap(exp0,exp2); }
			if (exp1>exp3) { math::Swap(k1,k3); math::Swap(exp1,exp3); }
			if (exp2>exp3) { math::Swap(k2,k3); math::Swap(exp2,exp3); }

			return ( (k0 + k1) + k2 ) +k3;
		};
		if (N==3) {
			T k0=_v[0]*p._v[0],	k1=_v[1]*p._v[1], k2=_v[2]*p._v[2];
			int exp0,exp1,exp2;

			frexp( double(k0), &exp0 );
			frexp( double(k1), &exp1 );
			frexp( double(k2), &exp2 );

			if( exp0<exp1 ) {
				if(exp0<exp2) return (k1+k2)+k0;
				         else return (k0+k1)+k2;
			} else {
				if(exp1<exp2) return (k0+k2)+k1;
				         else return (k0+k1)+k2;
			}
		};
	}  
//@}
//@{

  /** @name Cross products
  **/
		/// Cross product for 3D Point
	inline Point operator ^ ( Point const & p ) const
	{
		static_assert(N==3);
		return Point <S>
		(
			_v[1]*p._v[2] - _v[2]*p._v[1],
			_v[2]*p._v[0] - _v[0]*p._v[2],
			_v[0]*p._v[1] - _v[1]*p._v[0]
		);
	}
		/// Cross product for 2D Point
	  /// if called from a 3D or 4D points, returns the z component of the cross prod.
	inline S operator % ( Point const & p ) const
	{
		return _v[0]*p._v[1] - _v[1]*p._v[0];
	}
//@}
//@{

  /** @name Norms
  **/


		// Euclidean norm
	inline S Norm() const
	{
		if (N==2) return math::Sqrt( _v[0]*_v[0] + _v[1]*_v[1] );
		if (N==3) return math::Sqrt( _v[0]*_v[0] + _v[1]*_v[1] + _v[2]*_v[2] );
		if (N==4) return math::Sqrt( _v[0]*_v[0] + _v[1]*_v[1] + _v[2]*_v[2] + _v[4]*_v[4] );
	}
		// Squared Euclidean norm
	inline S SquaredNorm() const
	{
		if (N==2) return ( _v[0]*_v[0] + _v[1]*_v[1]  );
		if (N==3) return ( _v[0]*_v[0] + _v[1]*_v[1] + _v[2]*_v[2] );
		if (N==4) return ( _v[0]*_v[0] + _v[1]*_v[1] + _v[2]*_v[2] + _v[3]*_v[3] );
	}
		// Normalization (division by norm)
	inline Point & Normalize()
	{
    S n = Norm();
		if(n>0.0) (*this)/=n;
		return *this;
	}
		/// Homogeneous normalization (division by W)
	inline Point & HomoNormalize(){
		if (_v[3]!=0.0) {	_v[0] /= _v[3];	_v[1] /= _v[3];	_v[2] /= _v[3]; _v[3]=1.0; }
		return *this;
	};

//@}
		// Per component scaling
	inline Point & Scale( const Point & p )
	{
		_v[0] *= p._v[0];
		_v[1] *= p._v[1];
		if (N>2) _v[2] *= p._v[2];
		if (N>3) _v[3] *= p._v[3];
		return *this;
	}


		// Convert to polar coordinates
	void ToPolar( S & ro, S & tetha, S & fi ) const
	{
		ro = Norm();
		tetha = (S)atan2( _v[1], _v[0] );
		fi    = (S)acos( _v[2]/ro );
	}

//@{

  /** @name Comparison Operators. 
   Lexicographical order.
   **/

	inline bool operator == ( Point const & p ) const
	{
		if (N==2) return _v[0]==p._v[0] && _v[1]==p._v[1];
		if (N==3) return _v[0]==p._v[0] && _v[1]==p._v[1] && _v[2]==p._v[2];
		if (N==4) return _v[0]==p._v[0] && _v[1]==p._v[1] && _v[2]==p._v[2] && _v[3]==p._v[3];
	}
	inline bool operator != ( Point const & p ) const
	{
		if (N==2) return _v[0]!=p._v[0] || _v[1]!=p._v[1] ;
		if (N==3) return _v[0]!=p._v[0] || _v[1]!=p._v[1] || _v[2]!=p._v[2];
		if (N==4) return _v[0]!=p._v[0] || _v[1]!=p._v[1] || _v[2]!=p._v[2] || _v[3]!=p._v[3];
	}
	inline bool operator <  ( Point const & p ) const
	{
		return	(_v[2]!=p._v[2])?(_v[2]<p._v[2]):
				(_v[1]!=p._v[1])?(_v[1]<p._v[1]):
						       (_v[0]<p._v[0]);
	}
	inline bool operator >  ( Point const & p ) const
	{
		return	(_v[2]!=p._v[2])?(_v[2]>p._v[2]):
				(_v[1]!=p._v[1])?(_v[1]>p._v[1]):
							   (_v[0]>p._v[0]);
	}
	inline bool operator <= ( Point const & p ) const
	{
		return	(_v[2]!=p._v[2])?(_v[2]< p._v[2]):
				(_v[1]!=p._v[1])?(_v[1]< p._v[1]):
							   (_v[0]<=p._v[0]);
	}
	inline bool operator >= ( Point const & p ) const
	{
		return	(_v[2]!=p._v[2])?(_v[2]> p._v[2]):
				(_v[1]!=p._v[1])?(_v[1]> p._v[1]):
							   (_v[0]>=p._v[0]);
	}

	inline PointType LocalToGlobal(ParamType p) const{
		return *this;
	};
	
 //@}
}; // end class definition


template <class S>
inline S Angle( Point<3,S> const & p1, Point<3,S> const & p2 )
{
	S w = p1.Norm()*p2.Norm();
	if(w==0) return -1;
	S t = (p1*p2)/w;
	if(t>1) t = 1;
	else if(t<-1) t = -1;
    return (S) acos(t);
}

// versione uguale alla precedente ma che assume che i due vettori sono unitari
template <class S>
inline S AngleN( Point<3,S> const & p1, Point<3,S> const & p2 )
{
	S w = p1*p2;
	if(w>1) 
		w = 1;
	else if(w<-1) 
		w=-1;
  return (S) acos(w);
}


template <int N,class S>
inline S Norm( Point<N,S> const & p )
{
	return p.Norm();
}

template <int N,class S>
inline S SquaredNorm( Point<N,S> const & p )
{
    return p.SquaredNorm();
}

template <int N,class S>
inline Point<N,S> & Normalize( Point<N,S> & p )
{
    p.Normalize();
    return p;
}

template <int N, class S>
inline S Distance( Point<N,S> const & p1,Point<N,S> const & p2 )
{
    return (p1-p2).Norm();
}

template <int N, class S>
inline S SquaredDistance( Point<N,S> const & p1,Point<N,S> const & p2 )
{
    return (p1-p2).SquaredNorm();
}

	// Dot product preciso numericamente (solo double!!)
	// Implementazione: si sommano i prodotti per ordine di esponente
	// (prima le piu' grandi)
template<class S>
double StableDot ( Point<3,S> const & p0, Point<3,S> const & p1 )
{

}  

/// Computes a shape quality measure of the triangle composed by points p0,p1,p2
/// It Returns 2*AreaTri/(MaxEdge^2), 
/// the range is range [0.0, 0.866] 
/// e.g. Equilateral triangle sqrt(3)/2, halfsquare: 1/2, ... up to a line that has zero quality.
template<class S>
S Quality( Point<3,S> const &p0, Point<3,S> const & p1,  Point<3,S> const & p2)
{
	PointType<S> d10=p1-p0;
	PointType<S> d20=p2-p0;
	PointType<S> d12=p1-p2;
	PointType<S> x = d10^d20;

	S a = Norm( x );
	if(a==0) return 0; // Area zero triangles have surely quality==0;
	S b = SquaredNorm( d10 );
	S t = b;
	t = SquaredNorm( d20 ); if ( b<t ) b = t;
	t = SquaredNorm( d12 ); if ( b<t ) b = t;
	assert(b!=0.0);
	return a/b;
}

/// Returns the normal to the plane passing through p0,p1,p2
template<class S>
Point<3,S> Normal(const Point<3,S> & p0, const Point<3,S> & p1, const Point<3,S> & p2)
{
	return ((p1 - p0) ^ (p2 - p0));
}

/// Like the above, it returns the normal to the plane passing through p0,p1,p2, but normalized.
template<class S>
Point<3,S> NormalizedNormal(const Point<3,S> & p0, const Point<3,S> & p1, const Point<3,S> & p2)
{
	return ((p1 - p0) ^ (p2 - p0)).Normalize();
}


/// Point(p) Edge(v1-v2) dist, q is the point in v1-v2 with min dist
template<class S>
S PSDist( const Point<3,S> & p,
			         const Point<3,S> & v1,
					 const Point<3,S> & v2,
			         Point<3,S> & q )
{
    Point<3,S> e = v2-v1;
    S  t = ((p-v1)*e)/e.SquaredNorm();
    if(t<0)      t = 0;
	else if(t>1) t = 1;
	q = v1+e*t;
    return Distance(p,q);
}



/*template <class S>
inline Point<2,S>::Point ( const S nx, const S ny )
{_v[0]=nx;_v[1]=ny;};*/


/*template <class S>
inline Point<4,S>::Point ( const S nx, const S ny , const S nz , const S nw )
{_v[0]=nx;_v[1]=ny;_v[2]=nz;_v[3]=nw;};*/

/*template < class S>
	Point<3,S> Point<3,S>::operator * ( const S s ) const
	{
		return Point<3,S>( _v[0]*s, _v[1]*s , _v[2]*s );
	}*/

//template < class S>
	/*Point<3,double> Point<2,double>::operator * ( const double s ) const
	{
		return Point<2,double>( _v[0]*s, _v[1]*s );
	}*/

typedef Point<3,short>  Point3s;
typedef Point<3,int>	  Point3i;
typedef Point<3,float>  Point3f;
typedef Point<3,double> Point3d;
/*@}*/



} // end namespace
#endif

