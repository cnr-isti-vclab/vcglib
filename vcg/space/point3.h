/*#***************************************************************************
 * VCGLib                                                                    *
 *																																					 *
 * Visual Computing Group                                                o>  *
 * IEI Institute, CNUCE Institute, CNR Pisa                             <|   *
 *                                                                      / \  *
 * Copyright(C) 1999 by Paolo Cignoni, Claudio Rocchini                      *
 * All rights reserved.                                                      *
 *																																					 *
 * Permission  to use, copy, modify, distribute  and sell this  software and *
 * its documentation for any purpose is hereby granted without fee, provided *
 * that  the above copyright notice appear  in all copies and that both that *
 * copyright   notice  and  this  permission  notice  appear  in  supporting *
 * documentation. the author makes  no representations about the suitability *
 * of this software for any purpose. It is provided  "as is" without express *
 * or implied warranty.                                                      *
 *					                                                         *
 *****************************************************************************/
/*#**************************************************************************
  History
$Id: point3.h,v 1.2 2004-02-06 02:17:09 cignoni Exp $
$Log: not supported by cvs2svn $

****************************************************************************/

#pragma once
#ifndef __VCGLIB_POINT3
#define __VCGLIB_POINT3

//#include <limits>
#include <assert.h>
#ifndef __VCGLIB_UTILITY
#include <vcg/Utility.h>
#endif

namespace vcg {

    /** The class for representing a 3D point
     *  More details about this class.
     */

template <class T> class Point3
{
protected:
	T _v[3];

public:
	typedef T scalar;

	
		// Costruttori & assegnatori
	inline Point3 () { }
	inline Point3 ( const T nx, const T ny, const T nz )
	{
		_v[0] = nx;
		_v[1] = ny;
		_v[2] = nz;
	}
	inline Point3 ( Point3 const & p )
	{   
		_v[0]= p._v[0];
		_v[1]= p._v[1];
		_v[2]= p._v[2];
	}
	inline Point3 ( const T nv[3] )
	{
		_v[0] = nv[0];
		_v[1] = nv[1];
		_v[2] = nv[2];
	}
	inline Point3 & operator =( Point3 const & p )
	{
			_v[0]= p._v[0]; _v[1]= p._v[1]; _v[2]= p._v[2];
			return *this;
	}
	inline void zero()
	{
		_v[0] = 0;
		_v[1] = 0;
		_v[2] = 0;
	}
		// Accesso alle componenti
	inline const T &x() const { return _v[0]; } 
	inline const T &y() const { return _v[1]; }
	inline const T &z() const { return _v[2]; }
	inline T &x() { return _v[0]; }
	inline T &y() { return _v[1]; }
	inline T &z() { return _v[2]; }
	inline T & operator [] ( const int i )
	{
		assert(i>=0 && i<3);
		return _v[i];
	}
	inline const T & operator [] ( const int i ) const
	{
		assert(i>=0 && i<3);
		return _v[i];
	}
	inline const T * V() const
	{
		return _v;
	}
	inline T & V( const int i )
	{
		assert(i>=0 && i<3);
		return _v[i];
	}
	inline const T & V( const int i ) const
	{
		assert(i>=0 && i<3);
		return _v[i];
	}
		// Operatori matematici di base
	inline Point3 operator + ( Point3 const & p) const
	{
		return Point3<T>( _v[0]+p._v[0], _v[1]+p._v[1], _v[2]+p._v[2] );
	}
	inline Point3 operator - ( Point3 const & p) const
	{
		return Point3<T>( _v[0]-p._v[0], _v[1]-p._v[1], _v[2]-p._v[2] );
	}
	inline Point3 operator * ( const T s ) const
	{
		return Point3<T>( _v[0]*s, _v[1]*s, _v[2]*s );
	}
	inline Point3 operator / ( const T s ) const
	{
		return Point3<T>( _v[0]/s, _v[1]/s, _v[2]/s );
	}
		// dot product
	inline T operator * ( Point3 const & p ) const
	{
		return ( _v[0]*p._v[0] + _v[1]*p._v[1] + _v[2]*p._v[2] );
	}
		// Cross product
	inline Point3 operator ^ ( Point3 const & p ) const
	{
		return Point3 <T>
		(
			_v[1]*p._v[2] - _v[2]*p._v[1],
			_v[2]*p._v[0] - _v[0]*p._v[2],
			_v[0]*p._v[1] - _v[1]*p._v[0]
		);
	}

	inline Point3 & operator += ( Point3 const & p)
	{
		_v[0] += p._v[0];
		_v[1] += p._v[1];
		_v[2] += p._v[2];
		return *this;
	}
	inline Point3 & operator -= ( Point3 const & p)
	{
		_v[0] -= p._v[0];
		_v[1] -= p._v[1];
		_v[2] -= p._v[2];
		return *this;
	}
	inline Point3 & operator *= ( const T s )
	{
		_v[0] *= s;
		_v[1] *= s;
		_v[2] *= s;
		return *this;
	}
	inline Point3 & operator /= ( const T s )
	{
		_v[0] /= s;
		_v[1] /= s;
		_v[2] /= s;
		return *this;
	}
		// Norme
	inline T Norm() const
	{
		return Sqrt( _v[0]*_v[0] + _v[1]*_v[1] + _v[2]*_v[2] );
	}
	inline T SquaredNorm() const
	{
		return ( _v[0]*_v[0] + _v[1]*_v[1] + _v[2]*_v[2] );
	}
		// Scalatura differenziata
	inline Point3 & Scale( const T sx, const T sy, const T sz )
	{
		_v[0] *= sx;
		_v[1] *= sy;
		_v[2] *= sz;
		return *this;
	}
	inline Point3 & Scale( const Point3 & p )
	{
		_v[0] *= p._v[0];
		_v[1] *= p._v[1];
		_v[2] *= p._v[2];
		return *this;
	}

	// Normalizzazione
	inline Point3 & Normalize()
	{
		T n = Sqrt(_v[0]*_v[0] + _v[1]*_v[1] + _v[2]*_v[2]);
		if(n>0.0) {	_v[0] /= n;	_v[1] /= n;	_v[2] /= n;  }
		return *this;
	}
		// Polarizzazione
	void Polar( T & ro, T & tetha, T & fi ) const
	{
		ro = Norm();
		tetha = (T)atan2( _v[1], _v[0] );
		fi    = (T)acos( _v[2]/ro );
	}

		// Operatori di confronto (ordinamento lessicografico)
	inline bool operator == ( Point3 const & p ) const
	{
		return _v[0]==p._v[0] && _v[1]==p._v[1] && _v[2]==p._v[2];
	}
	inline bool operator != ( Point3 const & p ) const
	{
		return _v[0]!=p._v[0] || _v[1]!=p._v[1] || _v[2]!=p._v[2];
	}
	inline bool operator <  ( Point3 const & p ) const
	{
		return	(_v[2]!=p._v[2])?(_v[2]<p._v[2]):
				(_v[1]!=p._v[1])?(_v[1]<p._v[1]):
						       (_v[0]<p._v[0]);
	}
	inline bool operator >  ( Point3 const & p ) const
	{
		return	(_v[2]!=p._v[2])?(_v[2]>p._v[2]):
				(_v[1]!=p._v[1])?(_v[1]>p._v[1]):
							   (_v[0]>p._v[0]);
	}
	inline bool operator <= ( Point3 const & p ) const
	{
		return	(_v[2]!=p._v[2])?(_v[2]< p._v[2]):
				(_v[1]!=p._v[1])?(_v[1]< p._v[1]):
							   (_v[0]<=p._v[0]);
	}
	inline bool operator >= ( Point3 const & p ) const
	{
		return	(_v[2]!=p._v[2])?(_v[2]> p._v[2]):
				(_v[1]!=p._v[1])?(_v[1]> p._v[1]):
							   (_v[0]>=p._v[0]);
	}

	/// Questa funzione estende il vettore ad un qualsiasi numero di dimensioni
	/// paddando gli elementi estesi con zeri
	inline T Ext( const int i ) const
	{
		if(i>=0 && i<=2) return _v[i];
		else             return 0;
	}

	template <class Q>
	inline void Import( const Point3<Q> & b )
	{
		_v[0] = T(b[0]);
		_v[1] = T(b[1]);
		_v[2] = T(b[2]);
	}

	inline Point3 operator - () const
	{
		return Point3<T> ( -_v[0], -_v[1], -_v[2] );
	}

		// Casts
#ifdef __VCG_USE_CAST
inline operator Point3<int>			 (){ return Point3<int>			(_v[0],_v[1],_v[2]); }
inline operator Point3<unsigned int> (){ return Point3<unsigned int>(_v[0],_v[1],_v[2]); }
inline operator Point3<double>		 (){ return Point3<double>		(_v[0],_v[1],_v[2]); }
inline operator Point3<float>		 (){ return Point3<float>		(_v[0],_v[1],_v[2]); }
inline operator Point3<short>		 (){ return Point3<short>		(_v[0],_v[1],_v[2]); }
#endif

}; // end class definition



#ifdef __VCG_USE_P4_INTRINSIC__
#include <vcg/p4/point3p4.h>
#endif

/* Deprecata
template <class T>
inline Point3<T> operator * ( const T s, Point3<T> const & p )
{
    return Point3<T>( p._v[0] * s, p._v[1] * s, p._v[2] * s );
}
*/

#endif // include point3nt

	// ========== Parti comune alla vecchia e alla nuova implementazione ==================

template <class T>
inline T Angle( Point3<T> const & p1, Point3<T> const & p2 )
{
	T w = p1.Norm()*p2.Norm();
	if(w==0) return -1;
	T t = (p1*p2)/w;
	if(t>1) t = 1;
	else if(t<-1) t = -1;
    return (T) acos(t);
}

// versione uguale alla precedente ma che assume che i due vettori sono unitari
template <class T>
inline T AngleN( Point3<T> const & p1, Point3<T> const & p2 )
{
	T w = p1*p2;
	if(w>1) 
		w = 1;
	else if(w<-1) 
		w=-1;
  return (T) acos(w);
}


template <class T>
inline T Norm( Point3<T> const & p )
{
	return p.Norm();
}

template <class T>
inline T SquaredNorm( Point3<T> const & p )
{
    return p.SquaredNorm();
}

template <class T>
inline Point3<T> & Normalize( Point3<T> & p )
{
    p.Normalize();
    return p;
}

template <class T>
inline T Distance( Point3<T> const & p1,Point3<T> const & p2 )
{
    return (p1-p2).Norm();
}

template <class T>
inline T SquaredDistance( Point3<T> const & p1,Point3<T> const & p2 )
{
    return (p1-p2).SquaredNorm();
}

	// Dot product preciso numericamente (solo double!!)
	// Implementazione: si sommano i prodotti per ordine di esponente
	// (prima le piu' grandi)
template<class T>
double stable_dot ( Point3<T> const & p0, Point3<T> const & p1 )
{
	T k0 = p0._v[0]*p1._v[0];
	T k1 = p0._v[1]*p1._v[1];
	T k2 = p0._v[2]*p1._v[2];

	int exp0,exp1,exp2;

	frexp( double(k0), &exp0 );
	frexp( double(k1), &exp1 );
	frexp( double(k2), &exp2 );

	if( exp0<exp1 )
	{
		if(exp0<exp2)
			return (k1+k2)+k0;
		else
			return (k0+k1)+k2;
	}
	else
	{
		if(exp1<exp2)
			return(k0+k2)+k1;
		else
			return (k0+k1)+k2;
	}
}  

	// Returns 2*AreaTri/(MaxEdge^2), range [0.0, 0.866] 
	// e.g. halfsquare: 1/2, Equitri sqrt(3)/2, ecc
	// Modificata il 7/sep/00 per evitare l'allocazione temporanea di variabili
template<class T>
T Quality( Point3<T> const &p0, Point3<T> const & p1,  Point3<T> const & p2)
{
	Point3<T> d10=p1-p0;
	Point3<T> d20=p2-p0;
	Point3<T> d12=p1-p2;
	Point3<T> x = d10^d20;

	T a = Norm( x );
	if(a==0) return 0; // Area zero triangles have surely quality==0;
	T b = SquaredNorm( d10 );
	T t = b;
	t = SquaredNorm( d20 ); if ( b<t ) b = t;
	t = SquaredNorm( d12 ); if ( b<t ) b = t;
	assert(b!=0.0);
	return a/b;
}

		// Return the value of the face normal (internal use only)
template<class T>
Point3<T> Normal(const Point3<T> & p0, const Point3<T> & p1, const Point3<T> & p2)
{
	return ((p1 - p0) ^ (p2 - p0));
}

		// Return the value of the face normal (internal use only)
template<class T>
Point3<T> NormalizedNormal(const Point3<T> & p0, const Point3<T> & p1, const Point3<T> & p2)
{
	return ((p1 - p0) ^ (p2 - p0)).Normalize();
}

template<class T>
Point3<T> Jitter(Point3<T> &n, T RadAngle)
{
	Point3<T> rnd(1.0 - 2.0*T(rand())/RAND_MAX, 1.0 - 2.0*T(rand())/RAND_MAX, 1.0 - 2.0*T(rand())/RAND_MAX);
	rnd*=Sin(RadAngle);
	return (n+rnd).Normalize();
}



	// Point(p) Edge(v1-v2) dist, q is the point in v1-v2 with min dist
template<class T>
T PSDist( const Point3<T> & p,
			         const Point3<T> & v1,
					 const Point3<T> & v2,
			         Point3<T> & q )
{
    Point3<T> e = v2-v1;
    T  t = ((p-v1)*e)/e.SquaredNorm();
    if(t<0)      t = 0;
	else if(t>1) t = 1;
	q = v1+e*t;
    return Distance(p,q);
}

#if defined(FILE)
inline void print(Point3<int>    const & p, FILE * fp = stdout)   { fprintf(fp,"%d %d %d ",p.x(),p.y(),p.z()); }
inline void print(Point3<short>  const & p, FILE * fp = stdout)   { fprintf(fp,"%d %d %d ",p.x(),p.y(),p.z()); }
inline void print(Point3<float>  const & p, FILE * fp = stdout)   { fprintf(fp,"%g %g %g ",p.x(),p.y(),p.z()); }
inline void print(Point3<double> const & p, FILE * fp = stdout)   { fprintf(fp,"%g %g %g ",p.x(),p.y(),p.z()); }
#endif


typedef Point3<short>  Point3s;
typedef Point3<int>	   Point3i;
typedef Point3<float>  Point3f;
typedef Point3<double> Point3d;


} // end namespace
#endif

