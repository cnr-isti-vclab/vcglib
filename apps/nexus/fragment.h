#ifndef NXS_FRAGMENT_H
#define NXS_FRAGMENT_H

#include <vector>

#include <ptypes/pstreams.h>
#include "nexus.h"

namespace nxs {

class VoronoiPartition;

class NxsPatch {
 public:
  unsigned int patch;
  std::vector<vcg::Point3f> vert; 
  std::vector<unsigned short> face;
  std::vector<Link> bord;

  void write(pt::outstm *out);
  void read(pt::instm *in);
};

class Fragment {
 public:
  unsigned int id;

  float error;
  Nexus::Update update;
  std::vector<NxsPatch> pieces;
  
  void write(pt::outstm *out);
  void read(pt::instm *in);
};

 void join(Fragment &in, 
	   std::vector<vcg::Point3f> &newvert,
	   std::vector<unsigned int> &newface,
	   std::vector<Link> &newbord);

 void split(Fragment &out, 
	    std::vector<vcg::Point3f> &newvert,
	    std::vector<unsigned int> &newface,
	    std::vector<Link> &newbord, 
	    VoronoiPartition &part);
}

#endif
