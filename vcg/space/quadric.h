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

#ifndef __VCGLIB_QUADRIC
#define __VCGLIB_QUADRIC

#ifndef __VCGLIB_POINT3
#include <vcg/Point3.h>
#endif
#ifndef __VCGLIB_PLANE3
#include <vcg/Plane3.h>
#endif

namespace vcg {

template<class T>
class Quadric
{
public:
	T a[6];		// Matrice 3x3 simmetrica: a11 a12 a13 a22 a23 a33
	T b[3];		// Vettore r3
	T c;			// Fattore scalare (se -1 quadrica nulla)

	inline Quadric() { c = -1; }

	bool IsValid() const { return c>=0; }
	void SetInvalid() { c = -1.0; }

	void ByPlane( const Plane3<T> & p )					// Init dato un piano
	{
		a[0] =  p.n[0]*p.n[0];	// a11
		a[1] =  p.n[1]*p.n[0];	// a12 (=a21)
		a[2] =  p.n[2]*p.n[0];	// a13 (=a31)
		a[3] =  p.n[1]*p.n[1];	// a22
		a[4] =  p.n[2]*p.n[1];	// a23 (=a32)
		a[5] =  p.n[2]*p.n[2];	// a33
		b[0] = (T)(-2.0)*p.d*p.n[0];
		b[1] = (T)(-2.0)*p.d*p.n[1];
		b[2] = (T)(-2.0)*p.d*p.n[2];
		c    =  p.d*p.d;
	}

	void Zero()																// Azzera la quadrica
	{
		a[0] = 0;
		a[1] = 0;
		a[2] = 0;
		a[3] = 0;
		a[4] = 0;
		a[5] = 0;
		b[0] = 0;
		b[1] = 0;
		b[2] = 0;
		c    = 0;
	}

void operator = ( const Quadric & q )			// Assegna una quadrica
	{
		assert( IsValid() );
		assert( q.IsValid() );

		a[0] = q.a[0];
		a[1] = q.a[1];
		a[2] = q.a[2];
		a[3] = q.a[3];
		a[4] = q.a[4];
		a[5] = q.a[5];
		b[0] = q.b[0];
		b[1] = q.b[1];
		b[2] = q.b[2];
		c    = q.c;
	}

  void operator += ( const Quadric & q )			// Somma una quadrica
	{
		assert( IsValid() );
		assert( q.IsValid() );

		a[0] += q.a[0];
		a[1] += q.a[1];
		a[2] += q.a[2];
		a[3] += q.a[3];
		a[4] += q.a[4];
		a[5] += q.a[5];
		b[0] += q.b[0];
		b[1] += q.b[1];
		b[2] += q.b[2];
		c    += q.c;
	}

	T Apply( const Point3<T> & p ) const	// Applica la quadrica al punto p
	{
		assert( IsValid() );

	// Versione Lenta
/*
		Point3d t;
		t[0] = p[0]*a[0] + p[1]*a[1] + p[2]*a[2];
		t[1] = p[0]*a[1] + p[1]*a[3] + p[2]*a[4];
		t[2] = p[0]*a[2] + p[1]*a[4] + p[2]*a[5];
		double k = b[0]*p[0] + b[1]*p[1] + b[2]*p[2];
		double tp =t*p;
		return tp + k + c;
	
*/
	/* Versione veloce */

		return p[0]*p[0]*a[0] + 2*p[0]*p[1]*a[1] + 2*p[0]*p[2]*a[2] + p[0]*b[0] 
			  +   p[1]*p[1]*a[3] + 2*p[1]*p[2]*a[4] + p[1]*b[1]
			 +   p[2]*p[2]*a[5] + p[2]*b[2]	+ c;
	}

/// Draft version. It should be done in a more correctly way by using LRU decomposition.
bool Minimum(Point3<T> &x)
{	
		//T C[3][4];
		//C[0][0]=a[0]; C[0][1]=a[1]; C[0][2]=a[2];
		//C[1][0]=a[1]; C[1][1]=a[3]; C[1][2]=a[4];
		//C[2][0]=a[2]; C[2][1]=a[4]; C[2][2]=a[5];

		//C[0][3]=-b[0]/2;
		//C[1][3]=-b[1]/2;
		//C[2][3]=-b[2]/2;
		//return Gauss33(&(x[0]),C);

  Matrix33<T> mm;
  mm[0][0]=a[0]; mm[0][1]=a[1]; mm[0][2]=a[2];
	mm[1][0]=a[1]; mm[1][1]=a[3]; mm[1][2]=a[4];
	mm[2][0]=a[2]; mm[2][1]=a[4]; mm[2][2]=a[5];

  mm.Invert();
  x=mm*Point3<t>(-b[0]/2,-b[1]/2,-b[2]/2);
 return true; 

}

// determina il punto di errore minimo vincolato nel segmento (a,b)
bool Minimum(Point3<T> &x,Point3<T> &pa,Point3<T> &pb){
T	t1,t2, t4, t5, t8, t9, 
	t11,t12,t14,t15,t17,t18,t25,t26,t30,t34,t35,
	t41,t42,t44,t45,t50,t52,t54,
	t56,t21,t23,t37,t64,lambda;

	  t1 = a[4]*pb.z();      
	  t2 = t1*pa.y();
      t4 = a[1]*pb.y();
      t5 = t4*pa.x();
      t8 = a[1]*pa.y();
      t9 = t8*pa.x();
      t11 = a[4]*pa.z();
      t12 = t11*pa.y();
      t14 = pa.z()*pa.z();
      t15 = a[5]*t14;
      t17 = a[2]*pa.z();
      t18 = t17*pa.x();
      t21 = 2.0*t11*pb.y();
      t23 = a[5]*pb.z()*pa.z();
      t25 = a[2]*pb.z();
      t26 = t25*pa.x();
      t30 = a[0]*pb.x()*pa.x();
      t34 = 2.0*a[3]*pb.y()*pa.y();
      t35 = t17*pb.x();
      t37 = t8*pb.x();
      t41 = pa.x()*pa.x();
      t42 = a[0]*t41;
      t44 = pa.y()*pa.y();
      t45 = a[3]*t44;
      t50 = 2.0*t30+t34+2.0*t35+2.0*t37-(-b[2]/2)*pa.z()-(-b[0]/2)*pa.x()-2.0*t42-2.0*t45+(-b[1]/2)*pb.y()
+(-b[0]/2)*pb.x()-(-b[1]/2)*pa.y();
      t52 = pb.y()*pb.y();
      t54 = pb.z()*pb.z();
      t56 = pb.x()*pb.x();
      t64 = t5+t37-t9+t30-t18+t35+t26-t25*pb.x()+t2-t1*pb.y()+t23;
      lambda = (2.0*t2+2.0*t5+(-b[2]/2)*pb.z()-4.0*t9-4.0*t12-2.0*t15-4.0*t18+t21+2.0*t23+
2.0*t26+t50)/(-t45-a[3]*t52-a[5]*t54-a[0]*t56-t15-t42+t34-2.0*t12+t21-2.0*t4*pb.x()+
2.0*t64)/2.0;

	  if(lambda<0)  lambda=0;  else	  if(lambda>1)   lambda = 1;

		 x = pa*(1.0-lambda)+pb*lambda;		
		 return true;
	}

  void operator *= ( const T & w )			// Amplifica una quadirca
	{
		assert( IsValid() );
		
		a[0] *= w;
		a[1] *= w;
		a[2] *= w;
		a[3] *= w;
		a[4] *= w;
		a[5] *= w;
		b[0] *= w;
		b[1] *= w;
		b[2] *= w;
		c    *= w;
	}


};

typedef Quadric<short>  Quadrics;
typedef Quadric<int>	  Quadrici;
typedef Quadric<float>  Quadricf;
typedef Quadric<double> Quadricd;



} // end namespace

#endif
