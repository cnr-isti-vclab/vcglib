#ifndef NXS_TRISTRIP_H
#define NXS_TRISTRIP_H

#include <vector>

namespace nxs {

  void ComputeTriStrip(unsigned short nfaces, unsigned short *faces, 
		       std::vector<unsigned short> &strip);
}

#endif
