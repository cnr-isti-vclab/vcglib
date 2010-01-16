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
Revision 1.15  2007/09/21 11:34:10  ponchio
Just a clarification comment

Revision 1.14  2006/06/08 08:52:02  zifnab1974
gcc 4 needs the extra template keyword

Revision 1.13  2005/06/24 12:21:48  ponchio
Fixed "lerp" function.

Revision 1.12  2005/04/14 11:35:09  ponchio
*** empty log message ***

Revision 1.11  2004/09/09 12:51:28  fasano
corrected ColorRamp code (template specialization)

Revision 1.10  2004/09/09 08:39:33  cignoni
added a 'template<>' to the specialized constructors from a enum

Revision 1.9  2004/09/03 13:58:48  fasano
Corretto errore sintattico nelle specializzazioni parziali (float e char) di due costruttori di Color4

Revision 1.8  2004/07/15 11:01:43  ganovelli
added inclusion of  point3.h

Revision 1.7  2004/06/24 07:55:50  cignoni
Now color ramp can do reverse color ramp

Revision 1.6  2004/05/26 15:10:29  cignoni
Corrected bug in setgrayshade

Revision 1.5  2004/05/07 12:46:55  cignoni
added ifdef for gcc [Bug c++/14479]

Revision 1.4  2004/05/07 10:06:55  cignoni
Corrected template specialization syntax for gcc compiling

Revision 1.3  2004/03/10 21:38:40  cignoni
Written some documentation and added to the space module

Revision 1.2  2004/03/10 00:35:01  cignoni
Removed a wrong (?) copy constructor

Revision 1.1  2004/02/10 01:11:28  cignoni
Edited Comments and GPL license

****************************************************************************/

#ifndef __VCGLIB_COLOR4
#define __VCGLIB_COLOR4

#include <vcg/space/point3.h>
#include <vcg/space/point4.h>

namespace vcg {

/** \addtogroup space */
/*@{*/
    /**
        The templated class for representing 4 entity color.
        The class is templated over the ScalarType.  class that is used to represent color with float or with unsigned chars. All the usual
        operator overloading (* + - ...) is present.
     */
template <class T>
class Color4 : public Point4<T>
{
	typedef Point4<T> Base;
public:
  /// Constant for storing standard colors.
  /// Each color is stored in a simple in so that the bit pattern match with the one of Color4b.
	enum ColorConstant  {
		Black  =0xff000000,
		Gray		=0xff808080,
		White  =0xffffffff,

    Red    =0xff0000ff,
	  Green  =0xff00ff00,
	  Blue   =0xffff0000,

    Cyan   =0xffffff00,
		Yellow =0xff00ffff,
		Magenta=0xffff00ff,

    LightGray		=0xffc0c0c0,
		LightRed		=0xff8080ff,
		LightGreen  =0xff80ff80,
		LightBlue		=0xffff8080,

		DarkGray		=0xff404040,
		DarkRed		  =0xff000040,
		DarkGreen   =0xff004000,
		DarkBlue		=0xff400000
	};

  inline Color4 ( const T nx, const T ny, const T nz , const T nw ) :Point4<T>(nx,ny,nz,nw) {};
 // inline Color4 ( Color4 &c) :Point4<T>(c) {};
  inline Color4 ( const Point4<T> &c) :Point4<T>(c) {};
  inline Color4 (){};
  inline Color4 (ColorConstant cc);
  #ifdef VCG_USE_EIGEN
  template<typename OtherDerived>
  inline Color4 (const Eigen::MatrixBase<OtherDerived>& other) : Base(other)
  {
  	// TODO make sure the types are the same
  }
  #endif

  template <class Q>
	inline void Import(const Color4<Q> & b )
  {
	  (*this)[0] = T(b[0]);
	  (*this)[1] = T(b[1]);
	  (*this)[2] = T(b[2]);
	  (*this)[3] = T(b[3]);
  }

 template <class Q>
	inline void Import(const Point4<Q> & b )
  {
	  (*this)[0] = T(b[0]);
	  (*this)[1] = T(b[1]);
	  (*this)[2] = T(b[2]);
	  (*this)[3] = T(b[3]);
  }

 template <class Q>
  static inline Color4 Construct( const Color4<Q> & b )
  {
    return Color4(T(b[0]),T(b[1]),T(b[2]),T(b[3]));
  }

  //inline void Import(const Color4<float> &b);
  //inline void Import(const Color4<unsigned char> &b);

 inline Color4 operator + ( const Color4 & p) const
	{
		return Color4( (*this)[0]+p.V()[0], (*this)[1]+p.V()[1], (*this)[2]+p.V()[2], (*this)[3]+p.V()[3] );
	}


  inline void lerp(const Color4 &c0, const Color4 &c1, const float x);
	inline void lerp(const Color4 &c0, const Color4 &c1, const Color4 &c2, const Point3f &ip);
  /// given a float and a range set the corresponding color in the well known red->green->blue color ramp. To reverse the direction of the ramp just swap minf and maxf.
	inline void ColorRamp(const float &minf,const float  &maxf ,float v );

	inline void SetRGB( unsigned char r, unsigned char g, unsigned char b )
	{
		(*this)[0] = r;
		(*this)[1] = g;
		(*this)[2] = b;
		(*this)[3] = 0;
	}

	void SetHSVColor( float h, float s, float v){
	float r,g,b;
  if(s==0.0){	// gray color
		r = g = b = v;
		(*this)[0]=(unsigned char)(255*r);
		(*this)[1]=(unsigned char)(255*g);
		(*this)[2]=(unsigned char)(255*b);
		(*this)[3]=255;
		return;
	}
	if(h==1.0) h = 0.0;

	int i   = int( floor(h*6.0) );
	float f = float(h*6.0f - floor(h*6.0f));

	float p = v*(1.0f-s);
	float q = v*(1.0f-s*f);
	float t = v*(1.0f-s*(1.0f-f));

	switch(i){
			case 0: r=v; g=t; b=p; break;
			case 1: r=q; g=v; b=p; break;
			case 2: r=p; g=v; b=t; break;
			case 3: r=p; g=q; b=v; break;
			case 4: r=t; g=p; b=v; break;
			case 5: r=v; g=p; b=q; break;
			default: r=0;g=0;b=0; assert(0);break;
  }
		(*this)[0]=(unsigned char)(255*r);
		(*this)[1]=(unsigned char)(255*g);
		(*this)[2]=(unsigned char)(255*b);
		(*this)[3]=255;
//	V()[0]=r*256;V()[1]=g*256;V()[2]=b*256;
}

inline static Color4 GrayShade(float f)
{
 return Color4(f,f,f,1);
}

inline void SetGrayShade(float f)
{
 Import(Color4<float>(f,f,f,1));
}


/** Given an integer returns a well ordering of colors
// so that every color differs as much as possible form the previous one
// params:
//		n is the maximum expected value (max of the range)
//		v is the requested position
*/
inline static Color4 Scatter(int n, int a,float Sat=.3f,float Val=.9f)
{
  int b, k, m=n;
  int r =n;

    for (b=0, k=1; k<n; k<<=1)
			if (a<<1>=m) {
				if (b==0) r = k;
				b += k;
				a -= (m+1)>>1;
				m >>= 1;
			}
	else m = (m+1)>>1;
	if (r>n-b) r = n-b;

	//TRACE("Scatter range 0..%i, in %i out %i\n",n,a,b);
	Color4 rc;
	rc.SetHSVColor(float(b)/float(n),Sat,Val);
	return rc;
}

};
template <class T>
inline void Color4<T>::lerp(const Color4<T> &c0, const Color4<T> &c1, const float x)
{
	assert(x>=0);
	assert(x<=1);

	(*this)[0]=(T)(c1.V()[0]*x + c0.V()[0]*(1.0f-x));
	(*this)[1]=(T)(c1.V()[1]*x + c0.V()[1]*(1.0f-x));
	(*this)[2]=(T)(c1.V()[2]*x + c0.V()[2]*(1.0f-x));
	(*this)[3]=(T)(c1.V()[3]*x + c0.V()[3]*(1.0f-x));
}

template <class T>
inline void Color4<T>::lerp(const Color4<T> &c0, const Color4<T> &c1, const Color4<T> &c2, const Point3f &ip)
{
	assert(fabs(ip[0]+ip[1]+ip[2]-1)<0.00001);

	(*this)[0]=(T)(c0[0]*ip[0] + c1[0]*ip[1]+ c2[0]*ip[2]);
	(*this)[1]=(T)(c0[1]*ip[0] + c1[1]*ip[1]+ c2[1]*ip[2]);
	(*this)[2]=(T)(c0[2]*ip[0] + c1[2]*ip[1]+ c2[2]*ip[2]);
	(*this)[3]=(T)(c0[3]*ip[0] + c1[3]*ip[1]+ c2[3]*ip[2]);
}


template <class T>
inline void Color4<T>::ColorRamp(const float &minf,const float  &maxf ,float v )
{
  if(minf>maxf) { ColorRamp(maxf,minf,maxf+(minf-v)); return; }
	if(v <  minf ) { *this=Color4<T>(Color4<T>::Red); return; }
  //the case v > maxf is handled automatically at the end of the function

	float step=(maxf-minf)/4;
	v-=minf;
	if(v<step) {lerp(Color4<T>(Color4<T>::Red),  Color4<T>(Color4<T>::Yellow),v/step); return;}
	v-=step;
	if(v<step) {lerp(Color4<T>(Color4<T>::Yellow),Color4<T>(Color4<T>::Green),v/step);return;}
	v-=step;
	if(v<step) {lerp(Color4<T>(Color4<T>::Green),Color4<T>(Color4<T>::Cyan),v/step);  return;}
	v-=step;
	if(v<step) {lerp(Color4<T>(Color4<T>::Cyan),Color4<T>(Color4<T>::Blue),v/step);   return;}

	*this= Color4<T>(Color4<T>::Blue);

}


#if !defined(__GNUC__) || (__GNUC__ > 3)
template <>
#endif
template <>
inline void Color4<float>::Import(const Color4<unsigned char> &b)
{
  (*this)[0]=b[0]/255.0f;
  (*this)[1]=b[1]/255.0f;
  (*this)[2]=b[2]/255.0f;
  (*this)[3]=b[3]/255.0f;
}

#if !defined(__GNUC__) || (__GNUC__ > 3)
template <> // [Bug c++/14479] enum definition in template class with template methods causes error.
#endif
template <>
inline void Color4<unsigned char>::Import(const Color4<float> &b)
{
  (*this)[0]=(unsigned char)(b[0]*255.0f);
  (*this)[1]=(unsigned char)(b[1]*255.0f);
  (*this)[2]=(unsigned char)(b[2]*255.0f);
  (*this)[3]=(unsigned char)(b[3]*255.0f);
}

#if !defined(__GNUC__) || (__GNUC__ > 3)
template <> // [Bug c++/14479] enum definition in template class with template methods causes error.
#endif
template <>
inline void Color4<unsigned char>::Import(const Point4<float> &b)
{
  (*this)[0]=(unsigned char)(b[0]*255.0f);
  (*this)[1]=(unsigned char)(b[1]*255.0f);
  (*this)[2]=(unsigned char)(b[2]*255.0f);
  (*this)[3]=(unsigned char)(b[3]*255.0f);
}

#if !defined(__GNUC__) || (__GNUC__ > 3)
template <>
#endif
template <>
inline Color4<unsigned char> Color4<unsigned char>::Construct( const Color4<float> & b )
{
    return Color4<unsigned char>(
									(unsigned char)(b[0]*255.0f),
									(unsigned char)(b[1]*255.0f),
									(unsigned char)(b[2]*255.0f),
									(unsigned char)(b[3]*255.0f));
}

#if !defined(__GNUC__) || (__GNUC__ > 3)
template <>
#endif
template <>
inline Color4<float> Color4<float>::Construct( const Color4<unsigned char> & b )
{
    return Color4<float>(
									(float)(b[0])/255.0f,
									(float)(b[1])/255.0f,
									(float)(b[2])/255.0f,
									(float)(b[3])/255.0f);
}

//template <class T,class S>
//inline void Color4<T>::Import(const Color4<S> &b)
//{
//	V()[0] = T(b[0]);
//	V()[1] = T(b[1]);
//	V()[2] = T(b[2]);
//	V()[3] = T(b[3]);
//}
//
template<>
inline Color4<unsigned char>::Color4(Color4<unsigned char>::ColorConstant cc)
{
  *((int *)this )= cc;
}

template<>
inline Color4<float>::Color4(Color4<float>::ColorConstant cc)
{
  Import(Color4<unsigned char>((Color4<unsigned char>::ColorConstant)cc));
}

inline Color4<float> Clamp(Color4<float> &c)
{
	c[0]=math::Clamp(c[0],0.0f,1.0f);
	c[1]=math::Clamp(c[1],0.0f,1.0f);
	c[2]=math::Clamp(c[2],0.0f,1.0f);
	c[3]=math::Clamp(c[3],0.0f,1.0f);
	return c;
}

template<>
inline Color4<unsigned char> Color4<unsigned char>::operator + ( const Color4<unsigned char>  & p) const
{
		return Color4<unsigned char>(
									 (unsigned char)(math::Clamp(int((*this)[0])+int(p[0]),0,255)),
									 (unsigned char)(math::Clamp(int((*this)[1])+int(p[1]),0,255)),
									 (unsigned char)(math::Clamp(int((*this)[2])+int(p[2]),0,255)),
									 (unsigned char)(math::Clamp(int((*this)[3])+int(p[3]),0,255))
									 );
}


typedef Color4<unsigned char>  Color4b;
typedef Color4<float>  Color4f;
typedef Color4<double>  Color4d;

/*@}*/


} // end of NameSpace

#endif
