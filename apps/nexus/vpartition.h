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

#ifndef NXS_VPARTITION_H
#define NXS_VPARTITION_H

#include <vector>
#include <string>

#include <vcg/space/point3.h>


//TODO provide a Sort function, to sort spatially the seeds.

class ANNkd_tree;
class ANNbd_tree;
class ANNbruteForce;

namespace nxs {

class VPartition: public std::vector<vcg::Point3f> {
  public:
  VPartition(): bd(NULL) {}
  void Init();
  int Locate(const vcg::Point3f &p);

  //looks for the min distance point from seed.
  float Radius(unsigned int seed);
  void Closest(const vcg::Point3f &p, unsigned int nsize, 
	       std::vector<int> &nears, 
	       std::vector<float> &dist);
  void Closest(const vcg::Point3f &p, 
	       int &target, float &dist);
  

  
  ANNkd_tree *bd;
  std::vector<double> buffer;
  std::vector<double *> points;
};

} //namespace nxs
#endif
