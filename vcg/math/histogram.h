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

#ifndef __VCG_HISTOGRAM
#define __VCG_HISTOGRAM
//#include <vcg/Box2.h>
namespace vcg {

/**

*/
template <class ScalarType> 
class Histogram {
public:
  std::vector <int> H; /// counters for bins
	std::vector <ScalarType> R; /// Range for bins ( 
	ScalarType minv;
	ScalarType maxv;
	int n;  // numero di intervalli
	int cnt; // numero di campioni accumulati
	ScalarType avg;
	ScalarType rms;

	void SetRange(ScalarType _minv, ScalarType _maxv, int _n);
	void SetRange(ScalarType _minv, ScalarType _maxv, int _n, ScalarType gamma);
	int Interize(ScalarType val);
	void Add(ScalarType v);
	ScalarType Percentile(ScalarType frac) const;
	ScalarType Avg();
	ScalarType RMS();
	ScalarType Variance();
	ScalarType StandardDeviation();
	void Clear();
};

template <class ScalarType> 
void Histogram<ScalarType> ::Clear()
{
	H.clear();
	R.clear();
	cnt=0;
	avg=0;
	rms=0;
	n=0;
	minv=0;
	maxv=1;
}
template <class ScalarType> 
void Histogram<ScalarType>::SetRange(ScalarType _minv, ScalarType _maxv, int _n)
{
	Clear();
	minv=_minv;maxv=_maxv;n=_n;
	H.resize(n+1);
	fill(H.begin(),H.end(),0);
	R.resize(n+1);
	ScalarType dlt=(maxv-minv)/n;
	for(int i=0;i<n+1;++i)
		R[i]=minv+dlt*i;
}

template <class ScalarType> 
void Histogram<ScalarType>::SetRange(ScalarType _minv, ScalarType _maxv, int _n, ScalarType gamma)
{
	Clear();
	minv=_minv;maxv=_maxv;n=_n;
	H.resize(n);
	fill(H.begin(),H.end(),0);
	R.resize(n+1);
	double dlt=(maxv-minv);
	for(int i=0;i<n+1;++i)
		R[i]=minv+dlt*pow(double(i)/n,gamma);
}
template <class ScalarType> 
int Histogram<ScalarType>::Interize(ScalarType val) 
{
	int pos = lower_bound(R.begin(),R.end(),val) - R.begin() - 1;
	if(pos>n) pos=n;
	return pos;
}

template <class ScalarType> 
void Histogram<ScalarType>::Add(ScalarType v){
	int pos= lower_bound(R.begin(),R.end(),v)-R.begin()-1;
	if(pos<=n){
		++H[pos];
		++cnt;
		avg+=v;
		rms += v*v;
	}
}
}// end namespace
#endif	