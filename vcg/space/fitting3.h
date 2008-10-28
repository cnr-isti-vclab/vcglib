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
Revision 1.2  2005/10/13 14:59:57  ganovelli
versione con svd

Revision 1.1  2005/03/14 17:04:24  ganovelli
created


****************************************************************************/

#ifndef __VCGLIB_FITTING3
#define __VCGLIB_FITTING3

#include <vector>
#include <vcg/space/plane3.h>
#include <vcg/math/matrix44.h>
#include <vcg/math/matrix33.h>
#include <vcg/math/lin_algebra.h>


namespace vcg {

template <class S>
Point3<S> PlaneFittingPoints(  std::vector< Point3<S> > & samples,Plane3<S> &p){

	int j;
	Matrix44<S> m;m.SetZero();
	typename std::vector< Point3<S> > ::iterator i;
	
	Point3<S> c; c.SetZero();
	for(i = samples.begin(); i != samples.end(); ++i)
		c+=*i;
	c/=samples.size();

	for(i = samples.begin(); i != samples.end(); ++i){	
	Point3<S> p = (*i)-c;
	for(j = 0 ; j < 3;++j)
				*(Point3<S>*)&m[j][0] += p * p[j];
	}

	m[0][3]=	m[1][3]=m[2][3]=0.0;
	m[3][3]= 1.0;
	m[3][0]=	m[3][1]=m[3][2]=0.0;

	int n;
	Matrix44<S> res;
	Point4<S> e;
	Point3<S> d;
	Jacobi(m,e,res,n);

  //Sort eigenvalues (tarinisort)
  e[0] = fabs(e[0]);
  e[1] = fabs(e[1]);
  e[2] = fabs(e[2]);
  Point3<S> eval;
	int maxi,mini,medi;
	if (e[1] > e[0]) { maxi=1; mini=0; } else { maxi=0; mini=1;}
	if (e[maxi] < e[2]) maxi=2;	else if(e[mini] > e[2]) mini=2;
  medi = 3 - maxi -mini;
  eval = Point3<S>(e[mini], e[medi], e[maxi]);
  
	d[0]=res[0][mini];
	d[1]=res[1][mini];
	d[2]=res[2][mini];

	p.SetOffset(c.dot(d)/d.Norm());
	p.SetDirection(d/d.Norm());

  return eval;
}

template<class S>
inline double FIT_VExp( const Point3<S> & x, const int i )
{
	assert(i>=0);
	assert(i<4);
	if(i==0) return 1;
	else     return x[i-1];
}

	/** Fitting di piani: trova il piano che meglio approssima
	    l'insieme di punti dato
	 */
template<class S>
bool PlaneFittingPointsOld(  std::vector< Point3<S> > & samples, Plane3<S> & p )
{
	Point3<S> d;

  const int N = 4;
	S P[N][N];		// A = s' . s
	S U[N][N];
	int i,j,k,n;

	n = (int)samples.size();
	if(n<3)
		return false;

	//printf("\n p_prima: %f %f %f %f \n",p.Offset(),p.Direction()[0],p.Direction()[1],p.Direction()[2]);

	for(i=0;i<N;++i)
	{
		for(j=i;j<N;++j)
		{
			P[i][j] = 0;
			for(k=0;k<n;++k)
				P[i][j] += FIT_VExp(samples[k],i) * FIT_VExp(samples[k],j);
		}
		for(j=0;j<i;++j)
			P[i][j] = P[j][i];
	}

	//printf("D \n");
	//for(i=0;i<N;++i){
	//	printf("\n");
 //		for(j=0;j<N;++j)
	//		printf("%2.3f\t",P[i][j]);
	//}
	//
	Matrix44<S> m;
	for(i=0;i<N;++i)
 		for(j=0;j<N;++j)
			m[i][j]=P[i][j];


//	Point4<S> s;s.SetZero();
//
//	s.Normalize();
//	printf("\n RES %f %f %f %f \n",s[0],s[1],s[2],s[3]);
//printf("\n GJ \n");
//	for(i=0;i<N;++i){
//		printf("\n");
// 		for(j=0;j<N;++j)
//			printf("%2.3f\t",m[i][j]);
//	}
	for(i=0;i<N;++i)
	{
		U[i][i] = 1.0;
		for(j=0;j<i;++j)
			U[i][j] = 0.0;
		for(j=i+1;j<N;++j)
		{
			if(P[i][i]==0.0)
				return false;

			U[i][j] = P[i][j]/P[i][i];
			for(k=j;k<N;++k)
				P[j][k] -= U[i][j]*P[i][k];
		}
	}

	//printf("\n U \n");
	//for(i=0;i<N;++i){
	//	printf("\n");
 //		for(j=0;j<N;++j)
	//		printf("%2.3f\t",U[i][j]);
	//}

	
 	S norm = Point3<S>(U[1][2]*U[2][3]-U[1][3],-U[2][3],1).Norm();

 	p.SetDirection(Point3<S>(U[1][2]*U[2][3]-U[1][3],-U[2][3],1));
 	p.SetOffset(-(U[0][2]*U[2][3]-U[0][3]+U[0][1]*U[1][3]-U[0][1]*U[1][2]*U[2][3])/norm);

 
	//printf("\n p: %f %f %f %f \n",p.Offset(),p.Direction()[0],p.Direction()[1],p.Direction()[2]);

	return true;
}
} // end namespace

#endif
