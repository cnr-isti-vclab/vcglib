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

#ifndef NXS_PCHAIN_H
#define NXS_PCHAIN_H

#include <vcg/space/point3.h>

namespace {

template <class Partition> class PChain {
 public:
  typedef typename Partition::Key Key;

  unsigned int Levels() {
    return levels.size();
  }

  Key Locate(unsigned int level, const vcg::Point3f &p) {
    assert(level < levels.size());
    return levels[level].Locate(p);
  }

  float Priority(unsigned int level, const vcg::Point3f &p, Key key) {
    assert(level < levels.size());
    return levels[level].Priority(level, p, key);
  }
 private:
  std::vector<Partition> levels;
};

}//namespace

#endif
