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

#ifndef __VCGLIB_FITTING3
#define __VCGLIB_FITTING3

#include <vector>
#include <vcg/space/plane3.h>


namespace vcg {

	// Funzione di supporto: Ritorna il vettore 1 x y z
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
template<class POINT_TYPE, class PLANE_TYPE>
bool PlaneFittingPoints( const std::vector< POINT_TYPE > & samples, Plane3<typename POINT_TYPE::ScalarType> & p )
{
	typedef typename POINT_TYPE::ScalarType S;
	const int N = 4;

	S P[N][N];		// A = s' . s
	S U[N][N];
	int i,j,k,n;

	n = samples.size();
	if(n<3)
		return false;

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

	p.SetDirection(Point3<S>(U[1][2]*U[2][3]-U[1][3],-U[2][3],1));
	p.SetOffset(-(U[0][2]*U[2][3]-U[0][3]+U[0][1]*U[1][3]-U[0][1]*U[1][2]*U[2][3]));
	return true;
}

} // end namespace

#endif