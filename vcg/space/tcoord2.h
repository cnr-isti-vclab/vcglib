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

#ifndef __VCGLIB_TCOORD2__
#define __VCGLIB_TCOORD2__

#include <vcg/space/point2.h>

namespace vcg {


template<class T = float,int N = 1>
class TCoord2  
{
private:
	Point2<T> _t[N];
	short     _n[N];
public:
	
	inline T & u() { return _t[0][0]; }
	inline T & v() { return _t[0][1]; }
	inline const T & u() const { return _t[0][0]; }
	inline const T & v() const { return _t[0][1]; }
	inline T & u(const int i) { return _t[i][0]; }
	inline T & v(const int i) { return _t[i][1]; }
	inline const T & u(const int i) const { return _t[i][0]; }
	inline const T & v(const int i) const { return _t[i][1]; }

	inline short     & n() { return _n[0]; }
	inline const short n() const { return _n[0]; }

	inline short     & n(const int i) { return _n[i]; }
	inline const short n(const int i) const { return _n[i]; }
	
	inline Point2<T> & t(const int i) { return _t[i]; }
	inline Point2<T> t(const int i) const { return _t[i]; }

	inline Point2<T> & t() { return _t[0]; }
	inline Point2<T> t() const { return _t[0]; }
	
	inline bool operator == ( TCoord2 const & p ) const
		{
		 for(int i=0;i<N;++i)
			 if(p._t[i] != _t[i] || p._n[i] != _n[i]) return false;
		 return true;
		}
	enum { n_coords=N};
};


template<class T = float>
class TCoordSimple
{
private:
	Point2<T> _t;

	inline short & static_n() const
	{
		static short _n = 0;
		return _n;
	}

public:

	inline T & u() { return _t[0]; }
	inline T & v() { return _t[1]; }
	inline const T & u() const { return _t[0]; }
	inline const T & v() const { return _t[1]; }
	inline T & u(const int i) { assert(i==0); return _t[0]; }
	inline T & v(const int i) { assert(i==0); return _t[1]; }
	inline const T & u(const int i) const { assert(i==0); return _t[0]; }
	inline const T & v(const int i) const { assert(i==0); return _t[1]; }

	inline bool operator == ( TCoordSimple const & p ) const
		{
			return _t==p._t;
		}

	inline Point2<T> & t(const int i) { assert(i==0); return _t; }
	inline Point2<T> t(const int i) const { assert(i==0); return _t; }

	inline Point2<T> & t() { return _t; }
	inline Point2<T> t() const { return _t; }

	inline short & n() { assert(static_n()==0); return static_n(); }
	inline short   n() const { assert(static_n()==0); return 0; }

	inline short & n(const int i) { assert(i==0); return static_n(); }
	inline short   n(const int i) const { assert(i==0); return 0; }

	enum { n_coords=1};

};

#ifdef __GL_H__


//inline void glTexCoord(TCoord2<int> const & p)   { glTexCoord3iv(p.v);}
//inline void glTexCoord(Point3<short> const & p) { glTexCoord3sv(p.v);}
//inline void glTexCoord(Point3<float> const & p) { glTexCoord3fv(p.v);}
//inline void glTexCoord(Point3<double> const & p){ glTexCoord3dv(p.v);}

#endif

}

#endif