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
Revision 1.1  2004/02/10 01:11:28  cignoni
Edited Comments and GPL license

****************************************************************************/

#ifndef __VCGLIB_COLOR4
#define __VCGLIB_COLOR4

#include <vcg/space/point4.h>

namespace vcg {
template <class T> 
class Color4 : public Point4<T>
{
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
		DarkBlue		=0xff400000,
	};
	
  inline Color4 ( const T nx, const T ny, const T nz , const T nw ) :Point4<T>(nx,ny,nz,nw) {};
 // inline Color4 ( Color4 &c) :Point4<T>(c) {};
  inline Color4 (){};
  inline Color4 (ColorConstant cc);
  template <class Q>  
	inline void Import(const Color4<Q> & b )
  {
	  _v[0] = T(b[0]);
	  _v[1] = T(b[1]);
	  _v[2] = T(b[2]);
	  _v[3] = T(b[3]);
  }

  inline void Import(const Color4<float> &b);
  inline void Import(const Color4<unsigned char> &b);
	
  inline void lerp(const Color4 &c0, const Color4 &c1, const float x);
	inline void lerp(const Color4 &c0, const Color4 &c1, const Color4 &c2, const Point3f &ip);
	inline void ColorRamp(const float &minf,const float  &maxf ,float v );

	inline void SetRGB( unsigned char r, unsigned char g, unsigned char b )
	{
		_v[0] = r;
		_v[1] = g;
		_v[2] = b;
		_v[3] = 0;
	}

	void SetHSVColor( float h, float s, float v){
	float r,g,b;
  if(s==0.0){	// gray color
		r = g = b = v;
		_v[0]=(unsigned char)(255*r);_v[1]=(unsigned char)(255*g);_v[2]=(unsigned char)(255*b);
		_v[3]=255;
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
  }
		_v[0]=(unsigned char)(255*r);_v[1]=(unsigned char)(255*g);_v[2]=(unsigned char)(255*b);
		_v[3]=255;
//	_v[0]=r*256;_v[1]=g*256;_v[2]=b*256;
}

inline static Color4 GrayShade(float f)
{
 return ColorF(f,f,f,1);
}

inline void SetGrayShade(float f)
{
 *this = ColorF(f,f,f,1);
}


// Given an integer returns a well ordering of colors 
// so that every color differs as much as possible form the previous one
// params: 
//		n is the maximum expected value (max of the range)
//		v is the requested position

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

	_v[0]=(T)(c1._v[0]*x + c0._v[0]*(1.0f-x));
	_v[1]=(T)(c1._v[1]*x + c0._v[1]*(1.0f-x));
	_v[2]=(T)(c1._v[2]*x + c0._v[2]*(1.0f-x));
	_v[3]=(T)(c1._v[3]*x + c0._v[3]*(1.0f-x));
}

template <class T>
inline void Color4<T>::lerp(const Color4<T> &c0, const Color4<T> &c1, const Color4<T> &c2, const Point3f &ip)
{
	assert(fabs(ip[0]+ip[1]+ip[2]-1)<0.00001);
	
	c[0]=(T)(c0.c[0]*ip[0] + c1.c[0]*ip[1]+ c2.c[0]*ip[2]);
	c[1]=(T)(c0.c[1]*ip[0] + c1.c[1]*ip[1]+ c2.c[1]*ip[2]);
	c[2]=(T)(c0.c[2]*ip[0] + c1.c[2]*ip[1]+ c2.c[2]*ip[2]);
	c[3]=(T)(c0.c[3]*ip[0] + c1.c[3]*ip[1]+ c2.c[3]*ip[2]);
}


template <class T>
inline void Color4<T>::ColorRamp(const float &minf,const float  &maxf ,float v )
{
	if(v <  minf ) { *this=Color4(Color4::Red); return; }
	
	float step=(maxf-minf)/4;
	v-=minf;
	if(v<step) {lerp(Color4(Color4::Red),  Color4(Color4::Yellow),v/step); return;}
	v-=step;
	if(v<step) {lerp(Color4(Color4::Yellow),Color4(Color4::Green),v/step);return;}
	v-=step;
	if(v<step) {lerp(Color4(Color4::Green),Color4(Color4::Cyan),v/step);  return;}
	v-=step;
	if(v<step) {lerp(Color4(Color4::Cyan),Color4(Color4::Blue),v/step);   return;}

	*this= Color4(Color4::Blue);

}


inline void Color4<float>::Import(const Color4<unsigned char> &b)
{
  _v[0]=b[0]/255.0f;
  _v[1]=b[1]/255.0f;
  _v[2]=b[2]/255.0f;
  _v[3]=b[3]/255.0f;
}

inline void Color4<unsigned char>::Import(const Color4<float> &b)
{
  _v[0]=(unsigned char)(b[0]*255.0f);
  _v[1]=(unsigned char)(b[1]*255.0f);
  _v[2]=(unsigned char)(b[2]*255.0f);
  _v[3]=(unsigned char)(b[3]*255.0f);
}


//template <class T,class S> 
//inline void Color4<T>::Import(const Color4<S> &b)
//{
//	_v[0] = T(b[0]);
//	_v[1] = T(b[1]);
//	_v[2] = T(b[2]);
//	_v[3] = T(b[3]);
//}
//
inline Color4<unsigned char>::Color4<unsigned char>(Color4<unsigned char>::ColorConstant cc)
{
  *((int *)this )= cc; 
}

inline Color4<float>::Color4<float>(Color4<float>::ColorConstant cc)
{
  Import(Color4<unsigned char>((Color4<unsigned char>::ColorConstant)cc)); 
}


typedef Color4<unsigned char>  Color4b;
typedef Color4<float>  Color4f;

} // end of NameSpace

#endif
