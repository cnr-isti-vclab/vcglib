/****************************************************************************
* VCGLib                                                            o o     *
* Visual and Computer Graphics Library                            o     o   *
*                                                                _   O  _   *
* Copyright(C) 2004-2016                                           \/)\/    *
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

#ifndef _VCG_TETRA_TOPOLOGY
#define _VCG_TETRA_TOPOLOGY

namespace vcg {
namespace tetrahedron {
/** \addtogroup tetrahedron */
/*@{*/

/** Return a boolean that indicate if the j-th face of the tetra is a border.
    @param j Index of the face
    @return true if j is an face of border, false otherwise
*/
template <class TetraType>
inline bool IsTTBorder(TetraType const & t,  const int j )
{
  if(TetraType::HasTTAdjacency())
      return t.cTTp(j)==&t;
  assert(0);
  return true;
}

}
}

#endif
