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
Revision 1.6  2004/03/08 14:49:37  ponchio
Aggiunti un po di inline davanti alle funzioni

Revision 1.5  2004/03/04 00:21:00  ponchio
added Acos e Asin

Revision 1.4  2004/03/03 22:51:49  cignoni
changed math from class to template

Revision 1.3  2004/02/19 15:28:01  ponchio
*** empty log message ***

Revision 1.2  2004/02/13 02:18:57  cignoni
Edited Comments and GPL license


****************************************************************************/

#ifndef __VCGLIB_MATH_BASE
#define __VCGLIB_MATH_BASE

#include <math.h>
#include <assert.h>
#include <limits.h>

  #ifdef __BORLANDC__
    float sqrtf (float v) {return sqrt(v);}
    float fabsf (float v) {return fabs(v);}
  #endif     

namespace vcg {
namespace math {

 template <class SCALAR> 
 class MagnitudoComparer
    {
      public:
	    inline bool operator() ( const SCALAR a, const SCALAR b ) { return fabs(a)>fabs(b);  }
    };
  

  inline float Sqrt(const float v)   { return sqrtf(v); }
  inline float Abs(const float v)   { return fabsf(v); }
  inline float Cos(const float v)   { return cosf(v); }
  inline float Sin(const float v)   { return sinf(v); }
  inline float Acos(const float v)   { return acosf(v); }
  inline float Asin(const float v)   { return asinf(v); }

  inline double Sqrt(const double v)   { return sqrt(v); }
  inline double Abs(const double v)   { return fabs(v); }
  inline double Cos(const double v)   { return cos(v); }
  inline double Sin(const double v)   { return sin(v); }
  inline double Acos(const double v)   { return acos(v); }
  inline double Asin(const double v)   { return asin(v); }
  

	// max and min values for each scala type 
	// syntax: Max<float>::Value 

  template <class SCALAR> class Min {
		public: static const SCALAR Value;
	};  
	template <class SCALAR> class Max {
		public: static const SCALAR Value;
	};

  const char            Min<char           >::Value = -128; 
  const char            Max<char           >::Value = +127; 
  const unsigned char 	Min<unsigned char  >::Value = 0; 
  const unsigned char  	Max<unsigned char  >::Value = 255; 
  const short	          Min<short	         >::Value = (-32767 -1); 
  const short	          Max<short	         >::Value =  +32767; 
  const unsigned short	Min<unsigned short >::Value = 0; 
  const unsigned short	Max<unsigned short >::Value = 65535; 
  const int             Min<int            >::Value = (-2147483647 - 1); 
  const int             Max<int            >::Value =  +2147483647; 
  const unsigned int    Min<unsigned int   >::Value = 0; 
  const unsigned int    Max<unsigned int   >::Value = 4294967295; 
  const __int64         Min<__int64        >::Value = (-9223372036854775807i64 - 1); 
  const __int64         Max<__int64        >::Value =  +9223372036854775807i64; 
  const long            Min<long           >::Value = (-2147483647L -1L); 
  const long            Max<long           >::Value = 2147483647L; 
  const unsigned long   Min<unsigned long  >::Value = 0; 
  const unsigned long   Max<unsigned long  >::Value = +4294967295; 
  const float	          Min<float	         >::Value = -3.4E38F; 
  const float	          Max<float	         >::Value = +3.4E38F; 
  const long double     Min<long double    >::Value = -1.2E308; //E4931?
  const long double     Max<long double    >::Value = +1.2E308; 
  const double          Min<double         >::Value = -1.7E308; 
  const double          Max<double         >::Value = +1.7E308; 


/* Some <math.h> files do not define M_PI... */
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

template <class SCALAR> 
inline SCALAR  Clamp( const SCALAR & val, const SCALAR& minval, const SCALAR& maxval)
{
	if(val < minval) return minval;
	if(val > maxval) return maxval;
	return val;
}



inline float   ToDeg(const float &a){return a*180.0f/float(M_PI);}
inline float   ToRad(const float &a){return float(M_PI)*a/180.0f;}
inline double  ToDeg(const double &a){return a*180.0/M_PI;}
inline double  ToRad(const double &a){return M_PI*a/180.0;}

}	// End namespace
}	// End namespace


#endif