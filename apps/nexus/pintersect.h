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
Revision 1.2  2004/07/01 21:46:36  ponchio
Added header.


****************************************************************************/
#ifndef PARTITION_INTERSECT_H
#define PARTITION_INTERSECT_H

#include <map>
#include <set>
#include <string>
#include <vcg/space/point3.h>

namespace nxs {

  struct IPair {
    unsigned int fine;
    unsigned int coarse;
    bool operator<(const IPair &p) const {
      if(fine == p.fine) return coarse < p.coarse;
      return fine < p.fine;
    }
  };

  struct ICell {
    unsigned int index;
    unsigned int count;
  };



template <class Partition> class PIntersect {
 public:
  Partition *fine;
  Partition *coarse;
  
  //  typedef typename Partition::Key Key;
  
  unsigned int last_cell;
  std::map<IPair, ICell>  pairs;
  std::set<unsigned int> cells;
  
  PIntersect(Partition *f = NULL, Partition *c = NULL):
    fine(f), coarse(c), last_cell(0) {}

  void Init(Partition *f, Partition *c, unsigned int of_cell) {
    fine = f;
    coarse = c;
    last_cell = of_cell;
  }

  unsigned int Locate(const vcg::Point3f &p) {
    IPair pp;
    pp.fine = fine->Locate(p);
    pp.coarse = coarse->Locate(p);
    if(pairs.count(pp)) { //already there
      ICell &cell = pairs[pp];
      cell.count++;
      return cell.index;
    } else {
      ICell &cell = pairs[pp];
      cell.index = last_cell++;
      cell.count = 1;
      cells.insert(cell.index);
      return cell.index;
    }

  }

  unsigned int Relocate(const vcg::Point3f &p) {
    //TODO!!!!!!
  }
  void Prune(unsigned int threshold) {
    //TODO!!!!!!
  }

  unsigned int Count(unsigned int cell) {
    return cells.count(cell);
  }

  bool Save(const std::string &file) {
    FILE *fp = fopen((file + ".chi").c_str(), "wb+");
    if(!fp)
      return false;
    
    unsigned int n = pairs.size();
    fwrite(&n, sizeof(unsigned int), 1, fp);
	
    std::map<IPair, ICell>::iterator i;
    for(i = pairs.begin(); i != pairs.end(); i++) {
      IPair pair = (*i).first;
      ICell cell = (*i).second;
      fwrite(&pair, sizeof(IPair), 1, fp);
      fwrite(&cell, sizeof(ICell), 1, fp);
    }
    fclose(fp);
    return true;
  }
  
  bool Load(const std::string &file) {
    pairs.clear();
    cells.clear();
    FILE *fp = fopen((file + ".chi").c_str(), "rb");
    if(!fp)
      return false;
    unsigned int n;
    fread(&n, sizeof(unsigned int), 1, fp);
    
    IPair pair;
    ICell cell;
    for(unsigned int i = 0; i < n; i++) {
      fread(&pair, sizeof(IPair), 1, fp);
      fread(&cell, sizeof(ICell), 1, fp);
      pairs[pair] = cell;
      cells.insert(cell.index);
    }
    fclose(fp);
    return true;
  }
};
 
}//namespace
#endif
