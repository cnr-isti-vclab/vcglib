/****************************************************************************
* VCGLib                                                            o o     *
* Visual and Computer Graphics Library                            o     o   *
*                                                                _   O  _   *
* Copyright(C) 2006                                                \/)\/    *
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
Revision 1.2  2007/05/04 16:16:40  ganovelli
standardized to component style

Revision 1.1  2006/10/13 14:11:49  cignoni
first version

****************************************************************************/

#ifndef __VCG_FACE_PLUS_COMPONENT_RT
#define __VCG_FACE_PLUS_COMPONENT_RT

#include <vcg/space/plane3.h>

namespace vcg {
  namespace face {

template <class CoordType>
struct EdgePlaneInfo{
	CoordType edge[3];
	::vcg::Plane3<typename CoordType::ScalarType> plane;
	typename CoordType::ScalarType edgescale;
};

template <class T> class EdgePlane: public T {
public:
	typedef EdgePlaneInfo<typename T::VertexType::CoordType> EdgePlaneType;

  typename T::VertexType::CoordType &Edge(const int j) {
		return _ep.edge[j];
	}
  typename T::VertexType::CoordType  cEdge(const int j)const {
		return _ep.edge[j];
	}

	typename vcg::Plane3<typename T::VertexType::CoordType::ScalarType> &Plane() {
		return _ep.plane;
	}
  typename vcg::Plane3<typename T::VertexType::CoordType::ScalarType>  cPlane()const {
		return _ep.plane;
	}

  static bool HasEdgePlane()   {   return true; }

	static void Name(std::vector<std::string> & name){name.push_back(std::string("EdgePlane"));T::Name(name);}

private:

EdgePlaneType _ep;
};

  } // end namespace face
}// end namespace vcg
#endif
