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
Revision 1.1  2005/09/28 20:01:35  m_di_benedetto
First Commit.


****************************************************************************/

#ifndef __VCGLIB_WRAPUTILS_H
#define __VCGLIB_WRAPUTILS_H

// vcg headers
#include <vcg/math/base.h>
#include <vcg/space/point3.h>
#include <vcg/space/box3.h>

namespace vcg {

class EmptyClass {
public:
	typedef EmptyClass ClassType;
};

class GetPointerFunctor {
public:
	typedef GetPointerFunctor ClassType;

	template <class T>
	inline T * operator () (T & t) {
		return (&t);
	}

	template <class T>
	inline T * operator () (T * & t) {
		return (t);
	}
};

class GetBox3Functor {
public:
	template <class OBJTYPE, class SCALARTYPE>
	void operator () (const OBJTYPE & obj, Box3<SCALARTYPE> & box) {
		Box3<typename OBJTYPE::ScalarType> tb;
		obj.GetBBox(tb);
		box.Import(tb);
	}
};

class GetBarycenter3Functor {
public:
	template <class OBJTYPE, class SCALARTYPE>
	void operator () (const OBJTYPE & obj, Point3<SCALARTYPE> & bar) {
		bar.Import(obj.Barycenter());
	}
};

template <class SCALARTYPE>
inline Point3<SCALARTYPE> Abs(const Point3<SCALARTYPE> & p) {
	return (Point3<SCALARTYPE>(math::Abs(p[0]), math::Abs(p[1]), math::Abs(p[2])));
}

template <class SCALARTYPE>
inline Point3<SCALARTYPE> LowClampToZero(const Point3<SCALARTYPE> & p) {
	return (Point3<SCALARTYPE>(math::Max(p[0], (SCALARTYPE)0), math::Max(p[1], (SCALARTYPE)0), math::Max(p[2], (SCALARTYPE)0)));
}

} // end namespace vcg

#endif // #ifndef __VCGLIB_WRAPUTILS_H
