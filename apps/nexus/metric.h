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
Revision 1.5  2005/02/19 16:22:45  ponchio
Minor changes (visited and Cell)

Revision 1.4  2005/02/08 12:43:03  ponchio
Added copyright


****************************************************************************/

#ifndef NXS_METRIC_H
#define NXS_METRIC_H

#include <wrap/gui/frustum.h>
#include <vcg/space/sphere3.h>
#include "nexus.h"

namespace nxs {
  
  enum MetricKind { FRUSTUM, FLAT, DELTA };

  class Metric {
  public:
    vcg::Frustumf frustum;    
    bool culling;

    Metric(): culling(true) {}
    virtual void GetView() { frustum.GetView(); }
    virtual float GetError(Entry &entry, bool &visible) = 0;

  };
  
  class FlatMetric: public Metric {
  public:
    float GetError(Entry &entry, bool &visible) { 
      visible = true;
      return entry.error; 
    }
  };

  class FrustumMetric: public Metric {
  public:
    float GetError(Entry &entry, bool &visible) {
      visible = true;
      vcg::Sphere3f &sph = entry.sphere;
      float dist = (sph.Center() - frustum.ViewPoint()).Norm() - sph.Radius();

      if(dist < 0) return 1e20f;

      float error = entry.error/frustum.Resolution(dist);
      if(culling) {
	float remote = frustum.Remoteness(sph.Center(), sph.Radius());      
	if(remote > 0) {
	  visible = false;
	  error /= remote;
	} else if(entry.cone.Backface(sph, frustum.ViewPoint())) {
	  visible = false;
	}
      }
      return error;
    }
  };
}

#endif
