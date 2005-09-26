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

#ifndef __VCGLIB_AABBBINARYTREE_UTILS_H
#define __VCGLIB_AABBBINARYTREE_UTILS_H

// vcg headers
#include <vcg/math/base.h>
#include <vcg/space/point3.h>

/***************************************************************************************/

namespace vcg {

class EmptyClass {
public:
	typedef EmptyClass ClassType;
};

class ObjPtrIteratorFunctor {
public:
	typedef ObjPtrIteratorFunctor ClassType;

	template <class T>
	inline T * operator () (T & obj) {
		return (&obj);
	}

	template <class T>
	inline T * operator () (T * & obj) {
		return (obj);
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

#endif // #ifndef __VCGLIB_AABBBINARYTREE_UTILS_H
