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
Revision 1.9  2005/02/22 10:38:11  ponchio
Debug, cleaning and optimization.

Revision 1.8  2005/02/21 17:55:36  ponchio
debug debug debug

Revision 1.7  2005/02/20 18:07:01  ponchio
cleaning.

Revision 1.6  2005/02/19 10:45:04  ponchio
Patch generalized and small fixes.

Revision 1.5  2005/02/18 13:04:13  ponchio
Added patch reordering.

Revision 1.4  2005/02/08 12:43:03  ponchio
Added copyright


****************************************************************************/

#ifndef NXS_ALGO_H
#define NXS_ALGO_H

#include <set>
#include <vector>
#include "patch.h"
#include <vcg/space/sphere3.h>

namespace nxs {
  
  class Nexus;
  class Patch;


  struct ZEntry {
    unsigned int id;
    unsigned int pos;
    bool operator<(const ZEntry &e) const { return pos < e.pos; }
  };

  //for every patch return close by (sphere intersecting)
  //threshold is added to the distance to make sure we do not miss anything
  void Connect(Nexus &nexus, std::vector< std::set<unsigned int> > &close,
	       float threshold);

  void ComputeNormals(Nexus &nexus);

  void Unify(std::vector<vcg::Point3f> &points, 
	     std::vector<unsigned short> &faces,
	     std::vector<unsigned int> &vremap, float threshold);
//  void Unify(Nexus &nexus, float threshold);
  void ZSort(Nexus &nexus, std::vector<unsigned int> &forward,
	     std::vector<unsigned int> &backward);
  //  void TightSphere(vcg::Sphere3f &sphere, std::vector<vcg::Point3f> &points);
}

#endif
