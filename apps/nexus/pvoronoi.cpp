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
Revision 1.5  2004/09/28 10:26:21  ponchio
Rewrote.

Revision 1.4  2004/09/21 00:53:23  ponchio
Lotsa changes.

Revision 1.3  2004/08/27 00:39:28  ponchio
Rewrote.

Revision 1.2  2004/07/01 21:35:34  ponchio
int -> Key

Revision 1.1  2004/06/24 14:32:45  ponchio
Moved from wrap/nexus

Revision 1.2  2004/06/24 14:19:20  ponchio
Debugged

Revision 1.1  2004/06/23 17:17:46  ponchio
Created


****************************************************************************/

#pragma warning(disable:4786 4804 4244 4018 4267 4311)
#include <stdio.h>
#include <iostream>
#include "pvoronoi.h"

using namespace std;
using namespace vcg;
using namespace nxs;

bool Seed::Dist(const Point3f &point, float &mindist, 
		   Point3f &res) {
  float newdist = Distance(p, point) * weight; 
  if(newdist < mindist) {
    mindist = newdist;
    res = p;
    return true;
  } else 
    return false;
}

void VoronoiPartition::Init() {
  assert(size() > 0);
  for(iterator i = begin(); i != end(); i++)
    box.Add((*i).p);

  ug.SetBBox(box);
  ug.Set(*(vector<Seed> *)this);
}

float VoronoiPartition::Closest(const vcg::Point3f &p, 
				unsigned int &target, float radius) {
  Point3f res;
  float mindist = 1e20;
  target = 0xffffffff;
  
  Seed *nsp = ug.GetClosest(p, mindist, res);
  if(nsp) 
    target = nsp-&*begin();

  return mindist;
}

Point3f VoronoiPartition::FindBorder(vcg::Point3f &p, float radius) {
  Point3f a = p;
  unsigned int atarget = Locate(a);
  Point3f &seed = operator[](atarget).p;

  if((a - seed).Norm() < radius/100) return p; //Bad luck.

  Point3f dir = (a - seed).Normalize();
  Point3f b = seed + dir*radius*1.1;
  unsigned int btarget = Locate(b);

  if(atarget == btarget) {
    //probably nothing on the side we are looking gor;
    return p;
  }
  Point3f m;
  for(unsigned int i = 0; i < 10; i++) {
    m = (a + b)/2;
    unsigned int mtarget = Locate(m);
    if(mtarget == atarget) a = m;
    else if(mtarget == btarget) b = m;
    else break; //something in the middle
  }
  return m;
}
