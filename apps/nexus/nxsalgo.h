#ifndef NXS_ALGO_H
#define NXS_ALGO_H

#include <vector>

namespace nxs {
  
  class Nexus;

  void ComputeNormals(Nexus &nexus);
  void ComputeTriStrip(unsigned short nfaces, unsigned short *faces, 
		    std::vector<unsigned short> &strip);
}

#endif
