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
Revision 1.6  2004/09/21 00:53:23  ponchio
Lotsa changes.

Revision 1.5  2004/08/27 00:39:28  ponchio
Rewrote.

Revision 1.4  2004/07/20 14:17:51  ponchio
*** empty log message ***

Revision 1.3  2004/07/01 21:34:59  ponchio
int -> Key

Revision 1.2  2004/06/25 16:47:13  ponchio
Various debug

Revision 1.1  2004/06/24 14:32:45  ponchio
Moved from wrap/nexus

Revision 1.2  2004/06/24 14:19:20  ponchio
Debugged

Revision 1.1  2004/06/23 17:17:46  ponchio
Created


****************************************************************************/

#ifndef NXS_NET_GRID_H
#define NXS_NET_GRID_H

#include <vector>
#include <set>
#include <string>
#include <stdio.h>

#include <vcg/space/point3.h>
#include <vcg/space/box3.h>
#include <vcg/space/index/grid_static_ptr.h>

#include "crude.h"

//TODO provide a Sort function, to sort spatially the seeds.

namespace nxs {

  class Seed {
  public:
    vcg::Point3f p;  
    float weight;
    typedef float ScalarType;
    bool Dist(const vcg::Point3f & point, float &mindist, vcg::Point3f &res);
    void GetBBox(vcg::Box3f &b) {b.Set(p);}
    bool IsD() { return false; }

    Seed(): weight(1) {}
    Seed(const vcg::Point3f &point): p(point), weight(1) {}
    Seed(const vcg::Point3f &point, const float w): 
      p(point), weight(w) {}
    
    inline float Dist(const vcg::Point3f &q) const { 
      return weight * vcg::Distance(p,q); 
    }
    inline float SquaredDist(const vcg::Point3f &q) const  { 
      return weight * weight *vcg::SquaredDistance(p,q); 
    }  
  };
  

  class VoronoiPartition: public std::vector<Seed> {
  public:
    vcg::GridStaticPtr< std::vector<Seed> > ug;
    
    void SetBox(const vcg::Box3f &b) { box = b; }
    //call this before starting queries.
    void Init();
    float Closest(const vcg::Point3f &p, 
		  unsigned int &target, float radius = 0);
    
    unsigned int Locate(const vcg::Point3f &p) {
      unsigned int target;
      Closest(p, target);
      return target;
    } 
    vcg::Box3f box;
  };

}
#endif

