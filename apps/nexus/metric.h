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
