#ifndef NXS_FRAGMENT_H
#define NXS_FRAGMENT_H

#include <vector>

#include <ptypes/pstreams.h>
#include "nexus.h"

namespace nxs {

class VoronoiPartition;

struct BigLink {
  unsigned int start_vert;
  unsigned int end_patch;
  unsigned int end_vert;
  bool operator<(const BigLink &l) const {
    if(end_patch == l.end_patch) {
      if(start_vert == l.start_vert) {
	return end_vert < l.end_vert;
      } else
	return start_vert < l.start_vert;
    } else
      return end_patch < l.end_patch;
  }
};

class NxsPatch {
 public:
  //this fields is the patch number in the infragment
  //and the seeds id in the outfragment
  unsigned int patch;
  std::vector<vcg::Point3f> vert; 
  std::vector<unsigned short> face;

  std::vector<Link> bord;

  void Write(pt::outstm *out);
  void Read(pt::instm *in);
};

class Fragment {
 public:
  unsigned int id;

  float error;

  std::vector<vcg::Point3f> seeds;
  std::vector<unsigned int> seeds_id;

  std::vector<NxsPatch> pieces;
  
  bool Write(pt::outstm *out);
  bool Read(pt::instm *in);

  //returns the index of the seed
  unsigned int Locate(const vcg::Point3f &p);
};

 void Join(Fragment &in, 
	   std::vector<vcg::Point3f> &newvert,
	   std::vector<unsigned int> &newface,
	   std::vector<BigLink> &newbord);

 void Split(Fragment &out, 
	    std::vector<vcg::Point3f> &newvert,
	    std::vector<unsigned int> &newface,
	    std::vector<BigLink> &newbord);
}

#endif
