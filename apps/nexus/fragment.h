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
Revision 1.6  2005/02/08 12:43:03  ponchio
Added copyright


****************************************************************************/

#ifndef NXS_FRAGMENT_H
#define NXS_FRAGMENT_H

#include <vector>

#include <ptypes/pstreams.h>
#include "nexus.h"

namespace nxs {

class VoronoiPartition;

struct BigLink {
  unsigned int start_vert;
  unsigned int end_patch;
  unsigned int end_vert;
  bool operator<(const BigLink &l) const {
    if(end_patch == l.end_patch) {
      if(start_vert == l.start_vert) {
	return end_vert < l.end_vert;
      } else
	return start_vert < l.start_vert;
    } else
      return end_patch < l.end_patch;
  }
};

class NxsPatch {
 public:
  //this fields is the patch number in the infragment
  //and the seeds id in the outfragment
  unsigned int patch;
  vcg::Sphere3f sphere;
  ANCone3f cone;
  
  std::vector<vcg::Point3f> vert; 
  std::vector<unsigned short> face;
  
  std::vector<Link> bord;
  
  void Write(pt::outstm *out);
  void Read(pt::instm *in);
};
 
class Fragment {
 public:
  unsigned int id;

  float error;

  std::vector<vcg::Point3f> seeds;
  std::vector<unsigned int> seeds_id;

  std::vector<NxsPatch> pieces;
  
  bool Write(pt::outstm *out);
  bool Read(pt::instm *in);

  //returns the index of the seed
  unsigned int Locate(const vcg::Point3f &p);
};

 void Join(Fragment &in, 
	   std::vector<vcg::Point3f> &newvert,
	   std::vector<unsigned int> &newface,
	   std::vector<BigLink> &newbord);

 void Split(Fragment &out, 
	    std::vector<vcg::Point3f> &newvert,
	    std::vector<unsigned int> &newface,
	    std::vector<BigLink> &newbord);
}

#endif
