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

#include "preload.h"
#include "nexusmt.h"
#include <iostream>

using namespace std;
using namespace nxs;

 void Preload::execute() {

   total_disk = 0;
   disk = 0;
   
   assert(mt);
   while(!get_signaled()) {
     trigger.wait();
     lock.enter();
     while(!queue.size()) {
       trigger.reset();
       lock.leave();
       trigger.wait();
       lock.enter();
     }
     //TODO check we are not loading too much memory!
     assert(queue.size());
     Item &item = queue.back();
     if(item.error == 0 || mt->CanAdd(item)) {
       //we cannot flush since we are not in the openGL thread
       //and flushing includes VBO buffer flushing also.
       Entry &entry = (*mt)[item.id];
       if(!entry.patch)
	 disk += entry.disk_size;

       Patch &patch = mt->GetPatch(item.id, item.error, false);

       //test... make sure memory is in ram (if not on vbo that is.
       if(!entry.vbo_array)
	 total_disk += patch.Face(0)[0];

       queue.pop_back();
     } else
       queue.clear();
     lock.leave();
   }
 }
  
