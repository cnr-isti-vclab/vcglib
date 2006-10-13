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
****************************************************************************/

#ifndef __VCG_FACE_PLUS_COMPONENT_RT
#define __VCG_FACE_PLUS_COMPONENT_RT

#include <vcg/space/plane3.h>

namespace vcg {
  namespace face {

template <class T> class RTInfo: public T {
public:
  typename T::VertexType::CoordType edge[3];
  Plane3<typename T::VertexType::ScalarType> plane;
};

  } // end namespace vert
}// end namespace vcg
#endif
