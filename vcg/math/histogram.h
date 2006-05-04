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
Revision 1.13  2006/03/29 09:25:47  zifnab1974
extra includes necessary for compilation of meshlab on AMD 64 with gcc 3.4.5

Revision 1.12  2006/01/12 13:12:54  callieri
in FileWrite, added file closing after finishing

Revision 1.11  2005/09/16 11:51:23  cignoni
removed signed/unsigned warning

Revision 1.10  2005/06/17 00:54:55  cignoni
Corrected small bug in SetRange (H was resized to n instead of n+1)

Revision 1.9  2005/06/14 14:27:00  ganovelli
added include of algorithm

Revision 1.8  2005/06/10 14:59:39  cignoni
Added include assert.h and cast to ScalarType for a pow in SetRange() function.

Revision 1.7  2005/06/09 14:19:55  cignoni
Added typedef Histogramf and Histogramd

Revision 1.6  2005/06/07 09:37:33  ponchio
Added fabs() to variance, which can sometime be negative in case
of rounding errors (and sqrt chokes on it).

Revision 1.5  2005/06/07 07:44:08  cignoni
Added Percentile and removed small bug in Add

Revision 1.4  2005/04/04 10:48:35  cignoni
Added missing functions Avg, rms etc, now fully (almost) functional

Revision 1.3  2005/03/14 09:23:40  cignoni
Added missing include<vector>

Revision 1.2  2004/08/25 15:15:26  ganovelli
minor changes to comply gcc compiler (typename's and stuff)

Revision 1.1  2004/06/24 09:12:28  cignoni
Initial Release


****************************************************************************/

#ifndef __VCG_HISTOGRAM
#define __VCG_HISTOGRAM

#include <vector>
#include <algorithm>
#include <assert.h>
#include <string>
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
  ScalarType Avg()               { return avg/cnt;              }
  ScalarType RMS()               { return sqrt(rms/double(cnt));}
  ScalarType Variance()          { return fabs(rms/cnt-Avg()*Avg());  }
  ScalarType StandardDeviation() { return sqrt(Variance());     }

  void FileWrite(const std::string &filename);

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
	H.resize(n+1);
	fill(H.begin(),H.end(),0);
	R.resize(n+1);
	double dlt=(maxv-minv);
	for(int i=0;i<n+1;++i)
		R[i]=minv+dlt*pow(ScalarType(i)/n,gamma);
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
	if(pos>=0 && pos<=n){
		++H[pos];
		++cnt;
		avg+=v;
		rms += v*v;
	}
  
}


template <class ScalarType> 
void Histogram<ScalarType>::FileWrite(const std::string &filename)
{
	FILE *fp;
	fp=fopen(filename.c_str(),"w");
	for(unsigned int i=0;i<H.size();i++)
		fprintf (fp,"%12.8lf , %12.8lf \n",R[i],double(H[i])/cnt);
	fclose(fp);
}

template <class ScalarType> 
ScalarType Histogram<ScalarType>::Percentile(ScalarType frac) const
{
	if(H.size()==0 && R.size()==0) return 0;
	assert(frac>=0 && frac<=1);
	ScalarType sum=0,partsum=0;
	int isum=0;
    int i;
    for(i=0;i<n+1;i++)	{ sum+=H[i]; isum+=H[i];}
	
    assert(isum==cnt);
	assert(sum==cnt);

	sum*=frac;
	for(i=0;i<n;i++)	{
			partsum+=H[i];
			if(partsum>=sum) break;
		}
	return R[i+1];
}

typedef Histogram<double> Histogramd ;
typedef Histogram<float> Histogramf ;

}// end namespace
#endif	
