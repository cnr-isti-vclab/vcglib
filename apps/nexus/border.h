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
#ifndef NXS_BORDER_H
#define NXS_BORDER_H

namespace nxs {

struct Link {
  Link(): start_vertex(0xffff), end_vertex(0xffff), end_patch(0xffffffff) {}

  unsigned short start_vertex;
  unsigned short end_vertex;
  unsigned int end_patch;
  bool IsNull() { return start_vertex == 0xffff; }
};

class Border {
 public:
  unsigned int Size() { return size; }
  Link &operator[](unsigned int i) { return start[i]; }
 private:
  unsigned short size;
  Link *start;
};

}


#endif
