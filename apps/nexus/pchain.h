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
Revision 1.3  2004/07/20 14:06:02  ponchio
Changed filename saving.

Revision 1.2  2004/06/25 16:47:13  ponchio
Various debug

Revision 1.1  2004/06/24 14:32:45  ponchio
Moved from wrap/nexus

Revision 1.2  2004/06/24 14:19:20  ponchio
Debugged

Revision 1.1  2004/06/23 17:17:46  ponchio
Created


****************************************************************************/

#ifndef NXS_PCHAIN_H
#define NXS_PCHAIN_H

#include <vcg/space/point3.h>
#include "crude.h"

namespace nxs {

class PChain {
 public:

  //  virtual void Init(Crude &crude, unsigned int max_level) = 0;

  /* Return an unique uint for couple of patch level, level+1 */
  virtual unsigned int Locate(unsigned int level, const vcg::Point3f &p) = 0;


  virtual void RemapFaces(Crude &crude, VFile<unsigned int> &face_remap,
			  std::vector<unsigned int> &patch_faces) = 0;

  /* virtual unsigned int Levels() = 0;
  
  virtual unsigned int Locate(unsigned int level, 
			      const vcg::Point3f &p) = 0;

  virtual float Priority(unsigned int level, 
  const vcg::Point3f &p, unsigned int key) = 0;*/
};

}//namespace

#endif
