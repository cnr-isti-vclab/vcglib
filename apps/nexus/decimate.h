#ifndef NXS_DECIMATE_H
#define NXS_DECIMATE_H

#include <vector>
#include "border.h"
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
