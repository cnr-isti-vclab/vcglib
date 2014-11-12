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

#ifndef POLYGON_H
#define POLYGON_H

namespace vcg {

/*
 */

template<class PolygonType>
typename PolygonType::CoordType PolygonBarycenter(PolygonType &F)
{
  typename PolygonType::CoordType bary(0,0,0);
  for (int i=0;i<F.VN();i++)
    bary+=F.V(i)->P();

  bary/=(typename PolygonType::ScalarType)F.VN();
  return bary;
}

template<class PolygonType>
typename PolygonType::CoordType PolygonNormal(PolygonType &F)
{
  typename PolygonType::CoordType n(0,0,0);

  for (int i=0;i<F.VN();i++)
    n+=Normal(F.P0(i),F.P(1),F.P2()).Normalize();

  return n.Normalize();
}


 }
#endif // POLYGON_H
