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

#ifndef NXS_DECIMATE_H
#define NXS_DECIMATE_H

#include <vector>
#include "border.h"
#include "fragment.h"
#include <vcg/space/point3.h>
namespace nxs {

  enum Decimation { QUADRIC, CLUSTER };
  class BigLink;
  float Decimate(Decimation mode,
		 unsigned int target_faces, 
		 std::vector<vcg::Point3f> &newvert, 
		 std::vector<unsigned int> &newface,
		 std::vector<BigLink> &newbord);

}

#endif
