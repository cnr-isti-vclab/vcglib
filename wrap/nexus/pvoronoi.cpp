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

#pragma warning(disable:4786 4804 4244 4018 4267 4311)
#include "pvoronoi.h"
using namespace std;
using namespace vcg;
using namespace nxs;

bool Seed::Dist(const Point3f &point, float &mindist, 
		   Point3f &res) {
  float newdist = Distance(p, point); 
  if(newdist < mindist) {
    mindist = newdist;
    res = p;
    return true;
  } else 
    return false;
}

int VoronoiPartition::Add(vcg::Point3f p, float weight) {
  Seed ns(p,weight);	
  all_seeds.push_back(ns); 
  seedBuf.push_back(p);
  
  if(seedBuf.size() >= MAX_BUF) {
    for(unsigned int i = 0; i < seedBuf.size(); ++i)
      ug_seeds.push_back(seedBuf[i]);
    seedBuf.clear();
    ug.Set(ug_seeds);
  }
  return size();
}

float VoronoiPartition::Closest(const vcg::Point3f &p, int &target, float radius) {
  Point3f res;
  float mindist = 1e20;
  target = -1;
  
  if(ug_seeds.size()) {
    Seed *nsp = ug.GetClosest(p, mindist, res);
    if(nsp) {
      target = nsp-&*ug_seeds.begin();
    }
  }	
  
  for(unsigned int i=0;i<seedBuf.size();++i) {
    float dist = seedBuf[i].Dist(p);
    if(mindist > dist) {		
      target=ug_seeds.size()+i;
      mindist=dist;
    }
  }
  
  //assert(target >=0 );
  //assert (target < size()+seedBuf.size());
  return mindist;
}


void VoronoiPartition::iterator::operator++() {
  ++seed;
}
const VoronoiPartition::Key VoronoiPartition::iterator::operator*() {
  return seed;
}
bool VoronoiPartition::iterator::operator==(const VoronoiPartition::iterator &key) {
  return key.seed == seed;
}
bool VoronoiPartition::iterator::operator!=(const VoronoiPartition::iterator &key) {
  return key.seed != seed;
}

VoronoiPartition::iterator VoronoiPartition::begin() {
  iterator i;
  i.seed = 0;
  return i;
}
VoronoiPartition::iterator VoronoiPartition::end() {
  iterator i;
  i.seed = size();
  return i;
}
int VoronoiPartition::size() {
  return all_seeds.size();
}
void VoronoiPartition::clear() {
  all_seeds.clear();
  ug_seeds.clear();
  seedBuf.clear();
}
unsigned int VoronoiPartition::count(Key key) {
  return key > 0 && key < size();
}

VoronoiPartition::Key VoronoiPartition::Locate(const vcg::Point3f &p) {
  int target;
  Closest(p, target);
  assert(target != -1);
  return target;
}

float VoronoiPartition::Priority(const vcg::Point3f &p, Key key) {
  Seed &seed = all_seeds[key];
  return seed.Dist(p);
}
	
