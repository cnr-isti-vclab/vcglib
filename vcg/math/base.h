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
Revision 1.2  2004/02/13 02:18:57  cignoni
Edited Comments and GPL license


****************************************************************************/

#ifndef __VCGLIB_MATH_BASE
#define __VCGLIB_MATH_BASE

#include <math.h>
#include <assert.h>

  #ifdef __BORLANDC__
    float sqrtf (float v) {return sqrt(v);}
    float fabsf (float v) {return fabs(v);}
  #endif     

namespace vcg {

  template<class T>  class Math
  {
  public:
    static T inline Sqrt(const T v);
    static T inline Abs(const T v);
    static T inline Cos(const T v);
    static T inline Sin(const T v);
    
    
    static const T  MaxVal;
    static T ToDeg(const T &a);
    static T ToRad(const T &a);
    // Unspecialized members
    T Clamp( const T & val, const T& minval, const T& maxval);
    class MagnitudoComparer
    {
      public:
	    inline bool operator() ( const T a, const T b ) { return fabs(a)>fabs(b);  }
    };
  };

  float Math<float>::Sqrt(const float v) 
  { return sqrtf(v); }
  float Math<float>::Abs(const float v) 
  { return fabsf(v); }
  float Math<float>::Cos(const float v) 
  { return cosf(v); }
  float Math<float>::Sin(const float v) 
  { return sinf(v); }

  double Math<double>::Sqrt(const double v) 
  { return sqrt(v); }
  double Math<double>::Abs(const double v) 
  { return fabs(v); }
  double Math<double>::Cos(const double v) 
  { return cos(v); }
  double Math<double>::Sin(const double v) 
  { return sin(v); }
  
    
    
  const unsigned char   Math<unsigned char  >::MaxVal = 255;
  const char	          Math<char 	        >::MaxVal = 127; 
  const unsigned short	Math<unsigned short	>::MaxVal = 0xFFFFu; 
  const short	          Math<short	        >::MaxVal = 0x7FFF; 
  const float	          Math<float	        >::MaxVal = 3.4E38F; 
  const int             Math<int            >::MaxVal = 2147483647; 
  const long double     Math<long double    >::MaxVal = 1.2E308; 
  const double          Math<double         >::MaxVal = 1.7E308; 
  const __int64         Math<__int64        >::MaxVal = 9223372036854775807; 

/* Some <math.h> files do not define M_PI... */
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

template <class SCALAR> 
inline SCALAR  Math<SCALAR>::Clamp( const SCALAR & val, const SCALAR& minval, const SCALAR& maxval)
{
	if(val < minval) return minval;
	if(val > maxval) return maxval;
	return val;
}



inline float  Math<float>::ToDeg(const float &a){return a*180.0f/float(M_PI);}
inline float  Math<float>::ToRad(const float &a){return float(M_PI)*a/180.0f;}
inline double  Math<double>::ToDeg(const double &a){return a*180.0/M_PI;}
inline double  Math<double>::ToRad(const double &a){return M_PI*a/180.0;}

}	// End namespace


#endif