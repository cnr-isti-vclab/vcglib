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
Revision 1.1  2004/06/23 17:17:46  ponchio
Created


****************************************************************************/

#ifndef NXS_PCHAIN_H
#define NXS_PCHAIN_H

#include <stdio.h>
#include <vcg/space/point3.h>

/** Partition must be a class with a Key type, with 
    Levels, Locate, Priority, Save(FILE *), Load(FILE *)
    as in pvoronoi.h */
namespace {

template <class Partition> class PChain {
 public:
  typedef typename Partition::Key Key;
  std::vector<Partition> levels;

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
  bool Save(const std::string &file) {
    FILE *fp = fopen(file.c_str(), "wb+");
    if(!fp) return false;
    int n = Levels();
    fwrite(&n, sizeof(int), 1, fp);
    for(int i = 0; i < n; i++)
      levels[i].Save(fp);
    fclose(fp);
    return true;
  }
  bool Load(const std::string &file) {
    levels.clear(); 
    FILE *fp = fopen(file.c_str(), "rb");
    if(!fp) return false;
    int n;
    fread(&n, sizeof(int), 1, fp);
    levels.resize(n);
    for(int i = 0; i < n; i++) 
      levels[i].Load(fp);

    fclose(fp);
    return true;
  }
};

}//namespace

#endif
