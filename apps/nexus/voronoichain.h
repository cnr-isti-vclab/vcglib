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
  VoronoiChain():scaling(0.5), patch_size(1000), patch_threshold(300) {}
  virtual ~VoronoiChain() {}
  void Initialize(unsigned int psize, unsigned int pthreshold);
  void Init(Crude &crude);

  virtual unsigned int Locate(unsigned int level, const vcg::Point3f &p);
  void RemapFaces(Crude &crude, VFile<unsigned int> &face_remap,
		  std::vector<unsigned int> &patch_faces);

  void BuildLevel(Nexus &nexus, unsigned int offset);

  std::vector<VoronoiPartition> levels;

  //  coarse partition -> associated patches
  std::map<unsigned int, std::set<unsigned int> > newfragments;
  std::map<unsigned int, std::set<unsigned int> > oldfragments;
  // private:
  float scaling;
  unsigned int patch_size;
  unsigned int patch_threshold;

  float radius;
  vcg::Box3f box;
};
}

#endif
