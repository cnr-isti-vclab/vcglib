#ifndef NXS_ALGO_H
#define NXS_ALGO_H

#include <vector>

namespace nxs {
  
  class Nexus;
  class Patch;

  void ComputeNormals(Nexus &nexus);
  void ComputeTriStrip(unsigned short nfaces, unsigned short *faces, 
		    std::vector<unsigned short> &strip);
  void Reorder(unsigned int signature, nxs::Patch &patch);
}

#endif
