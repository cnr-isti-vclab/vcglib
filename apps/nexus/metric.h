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

#ifndef NXS_METRIC_H
#define NXS_METRIC_H

#include <wrap/gui/frustum.h>
#include <vcg/space/sphere3.h>
#include "nexus.h"

namespace nxs {
  
  enum MetricKind { FRUSTUM, FLAT, DELTA };

  class Metric {
  public:
    virtual void GetView() {}
    virtual float GetError(Entry &entry) = 0;
  };
  
  class FlatMetric: public Metric {
  public:
    float GetError(Entry &entry) { return entry.error; }
  };

  class FrustumMetric: public Metric {
  public:
    vcg::Frustumf frustum;

    virtual void GetView() { frustum.GetView(); }
    float GetError(Entry &entry) {
      vcg::Sphere3f &sphere = entry.sphere;
      float dist = (sphere.Center() - frustum.ViewPoint()).Norm() - sphere.Radius();
      //float dist = Distance(sphere, frustum.ViewPoint());
      if(dist < 0) 
	      return 1e20f;

      float remote = frustum.Remoteness(sphere.Center(), sphere.Radius());      
      if(remote > 0)
        return (entry.error/remote)/frustum.Resolution(dist);

      //if(frustum.IsOutside(sphere.Center(), sphere.Radius()))
	      //return -1;
     
      return entry.error/frustum.Resolution(dist);
    }
  };
}

#endif
