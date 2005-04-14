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
Revision 1.5  2005/02/17 15:39:44  ponchio
Reorderes statistics a bit.

Revision 1.4  2005/02/08 12:43:03  ponchio
Added copyright


****************************************************************************/

#ifndef NXS_PRELOAD_H
#define NXS_PRELOAD_H

#include <assert.h>

#include <vector>

#include <ptypes/pasync.h>

#include "extraction.h"

namespace nxs {

class NexusMt;

class Preload: public pt::thread{
 public:

  NexusMt *mt;
  
  pt::mutex lock;
  pt::trigger trigger;
  
  std::vector<Item> queue;

  unsigned int disk;     //kbytes readed from disk
  unsigned int disk_tri; //number of triangles readed from disk
  unsigned int total_disk;
  Preload(): thread(false), trigger(false, false) {}
  ~Preload() {
    waitfor();
  }
  
  void execute();
  
  void post(std::vector<Item> &patches) {
    trigger.reset();
    lock.enter();

    queue.reserve(patches.size());
    for(int i = patches.size() -1; i >= 0; i--)
      queue.push_back(patches[i]);

    trigger.post();
    lock.leave();
  }

  void cleanup() {}
};

}
#endif
