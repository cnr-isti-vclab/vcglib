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
Revision 1.1  2004/07/02 13:00:02  ponchio
Created.


****************************************************************************/
#ifndef NXS_BORDER_H
#define NXS_BORDER_H

namespace nxs {

struct Link {
  Link(): start_vert(0xffff), end_vert(0xffff), end_patch(0xffffffff) {}
  Link(unsigned short sv, unsigned short ev, unsigned int ep):
    start_vert(sv), end_vert(ev), end_patch(ep) {}
  unsigned short start_vert;
  unsigned short end_vert;
  unsigned int end_patch;
  bool IsNull() { return end_patch == 0xffffffff; }

  bool operator==(const Link &l) {
    return end_patch == l.end_patch && 
      end_vert == l.end_vert &&
      start_vert == l.start_vert;
  }
  bool operator<(const Link &l) {
    if(end_patch == l.end_patch) {
      if(start_vert == l.start_vert) {
	return end_vert < l.end_vert;
      } else
	return start_vert < l.start_vert;
    } else
      return end_patch < l.end_patch;
  }
};

class Border {
 public:
  Border(Link *l = NULL, unsigned short s = 0): start(l), size(s) {}
  unsigned int Size() { return size; }
  Link &operator[](unsigned int i) { return start[i]; }

  //TODO implement an iterator! 
 private:
  Link *start;
  unsigned short size;
};

}


#endif
