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
Revision 1.4  2005/02/21 17:55:48  ponchio
debug debug debug

Revision 1.3  2005/01/21 17:09:13  ponchio
Porting and debug.

Revision 1.2  2004/12/04 13:24:27  ponchio
Fixed a couple of memory leak...

Revision 1.1  2004/11/30 22:50:30  ponchio
Level 0.


****************************************************************************/

#include <stdio.h>

#include "vpartition.h"
#include <ANN/ANN.h>

using namespace std;
using namespace vcg;
using namespace nxs;

VPartition::~VPartition() {
  if(bd) delete bd;
}

void VPartition::Init() {  
  if(bd) delete bd;
  buffer.resize(size() * 3);
  for(unsigned int i = 0; i < size(); i++) {
    for(int k = 0; k < 3; k++)
      buffer[i*3+k] = operator[](i)[k];
  }
  points.resize(size());
  for(unsigned int i = 0; i < size(); i++) {
    points[i] = &buffer[i*3];
  }
  bd = new ANNkd_tree(&*points.begin(), size(), 3);
}

void VPartition::Closest(const vcg::Point3f &p, unsigned int nsize, 
			       vector<int> &nears, 
			       vector<float> &dist) {
  double point[3];
  point[0] = p[0];
  point[1] = p[1];
  point[2] = p[2];
  if(nsize > size()) nsize = size();

  nears.resize(nsize);
  dist.resize(nsize);
  vector<double> dists;
  dists.resize(nsize);
  bd->annkSearch(&point[0], nsize, &*nears.begin(), &*dists.begin());
  for(unsigned int i = 0; i < nsize; i++)
    dist[i] = (float)dists[i];
}

void VPartition::Closest(const vcg::Point3f &p, 
			       int &target, float &dist) {
  double point[3];
  point[0] = p[0];
  point[1] = p[1];
  point[2] = p[2];
  double dists;
  bd->annkSearch(&point[0], 1, &target, &dists);
  assert(target >= 0);
  assert(target < size());

  dist = (float)dists;
}

void VPartition::Closest(const vcg::Point3f &p, 
			 vector<int> &targets,
			 vector<double> &dists,
			 float max_distance) {

  double point[3];  point[0] = p[0];  point[1] = p[1];  point[2] = p[2];

  int seeds = 6;
  while(1) {
    if(seeds > size()) seeds = size();
    targets.resize(seeds);
    dists.resize(seeds);
    bd->annkSearch(&point[0], seeds, &(targets[0]), &(dists[0]));
    for(int i = 0; i < seeds; i++) {
      if(dists[i] > max_distance) {
	targets.resize(i);
	dists.resize(i);
	break;
      }
    }
    if(targets.size() < seeds) break;
    if(seeds == size()) break;
    seeds *= 2;
  }
}

void VPartition::Closest(const vcg::Point3f &p, unsigned int nsize,
			 int *targets, 
			 double *dists) {
  double point[3];
  point[0] = p[0];
  point[1] = p[1];
  point[2] = p[2];
  bd->annkSearch(&point[0], nsize, targets, dists);
}

int VPartition::Locate(const vcg::Point3f &p) {

  double point[3];
  point[0] = p[0];
  point[1] = p[1];
  point[2] = p[2];

  int target = -1;
  double dists;
  bd->annkSearch(&point[0], 1, &target, &dists);

  return target;
} 

float VPartition::Radius(unsigned int seed) {
  assert(size() > 1);
  int nears[2];
  double dists[2];

  double point[3];
  Point3f &p = operator[](seed);
  point[0] = p[0];
  point[1] = p[1];
  point[2] = p[2];

  bd->annkSearch(&point[0], 2, nears, dists);

  if(dists[1] == 0) return 0.0f;
  assert(nears[0] == seed);
  assert(dists[0] == 0);
  
  return (float)sqrt(dists[1]);
}
