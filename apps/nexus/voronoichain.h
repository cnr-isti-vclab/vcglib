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
Revision 1.3  2004/09/21 00:53:23  ponchio
Lotsa changes.

Revision 1.2  2004/09/16 14:25:16  ponchio
Backup. (lot of changes).

Revision 1.1  2004/08/26 18:03:48  ponchio
First draft.


****************************************************************************/

#ifndef NXS_VORONOI_CHAIN_H
#define NXS_VORONOI_CHAIN_H

#include <vector>
#include <map>

#include <vcg/space/box3.h>

#include "pchain.h"
#include "pvoronoi.h"
#include "crude.h"
#include "nexus.h"
#include "vfile.h"

namespace nxs {

class VoronoiChain: public PChain {
 public:
  unsigned int mean_size; //mean number of faces per patch
  unsigned int min_size;  //minimum number of faces per patch (soft)
  unsigned int max_size;  //max number of faces per patch (hard);

  VoronoiChain(unsigned int mean_s,
	       unsigned int min_s,
	       unsigned int max_s = 32000):
    mean_size(mean_s), min_size(min_s), max_size(max_s) {}
  virtual ~VoronoiChain() {}

  void Init(Crude &crude, float scaling, int steps);
  virtual unsigned int Locate(unsigned int level, const vcg::Point3f &p);
  void RemapFaces(Crude &crude, VFile<unsigned int> &face_remap,
		  std::vector<unsigned int> &patch_faces, float scaling, 
		  int steps);

  void BuildLevel(Nexus &nexus, unsigned int offset, float scaling, int steps);

  std::vector<VoronoiPartition> levels;

  //  coarse partition -> associated patches
  std::map<unsigned int, std::set<unsigned int> > newfragments;
  std::map<unsigned int, std::set<unsigned int> > oldfragments;
  // private:

  float radius;
  vcg::Box3f box;
};
}

#endif
