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
Revision 1.11  2004/03/10 17:37:54  tarini
Added Atan2.
Added common utilities: Max, Min, Swap, Sort(a,b), Sort(a,b,c).
Changed Max values syntax. example:  Value<float>::Max

Revision 1.10  2004/03/10 16:54:57  tarini
Added Atan2.
Added common utilities: Max, Min, Swap, Sort(a,b), Sort(a,b,c).
Changed Max values syntax. example:  Value<float>::Max

Revision 1.7  2004/03/08 19:38:29  tarini
Added Min e Max. usage: Min<float>::Value (tarini)

Revision 1.6  2004/03/08 14:49:37  ponchio
Aggiunti un po di inline davanti alle funzioni

Revision 1.5  2004/03/04 00:21:00  ponchio
added Acos e Asin

Revision 1.4  2004/03/03 22:51:49  cignoni
changed math from class to template

Revision 1.2  2004/02/13 02:18:57  cignoni
Edited Comments and GPL license


****************************************************************************/

#ifndef __VCGLIB_MATH_BASE
#define __VCGLIB_MATH_BASE

#include <math.h>
#include <assert.h>
#include <limits.h>

/// static_assert: implemented as a macro for "assert", but it is separated for clarity.
/// Should be used for checking integrity constraints that can be tested at complile time,
/// as the ones involving templated constants in templated classes.
#define static_assert assert

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
  inline float Atan2(const float v0,const float v1)   { return atan2f(v0,v1); }

  inline double Sqrt(const double v)   { return sqrt(v); }
  inline double Abs(const double v)   { return fabs(v); }
  inline double Cos(const double v)   { return cos(v); }
  inline double Sin(const double v)   { return sin(v); }
  inline double Acos(const double v)   { return acos(v); }
  inline double Asin(const double v)   { return asin(v); }
  inline double Atan2(const double v0,const double v1)   { return atan2(v0,v1); }
  

	/// max and min values for each scala type 
	/// syntax: Max<float>::Value  

  template <class SCALAR> class Value {
		public: static const SCALAR Min; static const SCALAR Max;
	};  

  const char            Value<char           >::Min = -128; 
  const char            Value<char           >::Max = +127; 
  const unsigned char 	Value<unsigned char  >::Min = 0; 
  const unsigned char  	Value<unsigned char  >::Max = 255; 
  const short	          Value<short	         >::Min = (-32767 -1); 
  const short	          Value<short	         >::Max =  +32767; 
  const unsigned short	Value<unsigned short >::Min = 0; 
  const unsigned short	Value<unsigned short >::Max = 65535; 
  const int             Value<int            >::Min = (-2147483647 - 1); 
  const int             Value<int            >::Max =  +2147483647; 
  const unsigned int    Value<unsigned int   >::Min = 0; 
  const unsigned int    Value<unsigned int   >::Max = 4294967295; 
  const __int64         Value<__int64        >::Min = (-9223372036854775807i64 - 1); 
  const __int64         Value<__int64        >::Max =  +9223372036854775807i64; 
  const long            Value<long           >::Min = (-2147483647L -1L); 
  const long            Value<long           >::Max = 2147483647L; 
  const unsigned long   Value<unsigned long  >::Min = 0; 
  const unsigned long   Value<unsigned long  >::Max = +4294967295; 
  const float	          Value<float	         >::Min = -3.4E38F; 
  const float	          Value<float	         >::Max = +3.4E38F; 
  const long double     Value<long double    >::Min = -1.2E308; //E4931?
  const long double     Value<long double    >::Max = +1.2E308; 
  const double          Value<double         >::Min = -1.7E308; 
  const double          Value<double         >::Max = +1.7E308; 

	template<class T> inline Min(const T &a, const T &b){
		if (a<b) return a; else return b;
	};
	template<class T> inline Max(const T &a, const T &b){
		if (a<b) return b; else return a;
	};
	template<class T> inline Swap(T &a, T &b){
		T tmp=a; a=b; b=tmp;
	};
	template<class T> inline Sort(T &a, T &b){
		if (a>b) Swap(a,b);
	};
	template<class T> inline Sort(T &a, T &b, T &c){
		if (a>b) Swap(a,b);
		if (b>c) {Swap(b,c); if (a>b) Swap(a,b);}
	};

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

}	// End math namespace

/// a type that stands for "void". Useful for Parameter type of a point.
class VoidType{ public:
	VoidType(){};
};

}	// End vcg namespace


#endif