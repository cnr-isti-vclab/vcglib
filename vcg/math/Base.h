/*****************************************************************************
 * VCGLib                                                                    *
 *					                                                                 *
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
 *																																					 *
 *****************************************************************************/
/****************************************************************************
  History

 1999 Feb 02 First Draft.

 1999 May 15 Corrected Scope of sqrt.. (added ::)

 2000	Jan 26 inserito include condizionale
						 corretto Distance() e init()
			Jan 28 aggiunti Sin e Cos (per quaternion!) (con assert se uno li usa con int!)			
			Jun 26 Aggiunto Gauss33 (prima stava in lfield3) 
						 Nota: era stato scritto con Fabs invece di Abs....
			Jun 30 Aggiunto TRACE e Definizione M_PI
			Jul  4 aggiunto un cast a FL_TYPE in Gauss33;
			     5 Tolti i parametri inutili (warning 4100)
					 6 Aggiunto include assert.h
 2001 Jul 17 TRACE Release compilabile in unix (CR)
 2002 Jan 17 Aggiunte conversioni Radianti in Gradi
      Feb 27 Aggiunto tipo per Generica funzione di Callback
			Jul 28 Aggiunta seconda callback (utile per progress bar) (pc)
 2003 Jan 09 Aggiunta Clamp (pc)
 2003 Jan 17 Aggiunta MaxVal (mt)
      Sep 10 [BCB] Ridefinite sqrtf e fabsf per C++ Builder 
			    19 Aggiunto suffisso 'u' per evitare un warning (pc)

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
    
    static const T  MaxVal;
    static T ToDeg(const T &a);
    static T ToRad(const T &a);
    // Unspecialized members
    SCALAR Clamp( const SCALAR & val, const SCALAR& minval, const SCALAR& maxval);
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
  double Math<double>::Sqrt(const double v) 
  { return sqrt(v); }
  double Math<double>::Abs(const double v) 
  { return fabs(v); }
    
    
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