#ifndef NXS_NEXUS_H
#define NXS_NEXUS_H

#include <string>
#include <vector>
#include <set>
#include <vcg/space/point3.h>
#include <vcg/space/sphere3.h>
#include "vfile.h"
#include "patch.h"
#include "border.h"

namespace nxs {



class Nexus {
 public:

  struct Entry {
    Entry(): patch_start(0xffffffff), border_start(0xffffffff),
	 patch_size(0), border_size(0), border_used(0),
      nvert(0), nface(0), sphere(vcg::Sphere3f()) {}

    unsigned int patch_start;  //granularita' Chunk
    unsigned int border_start; //granuralita' Link
    unsigned short patch_size;  //in cuhnks
    unsigned short border_size; //in Links
    unsigned short border_used; //in Links

    //Data used for extraction
    unsigned short nvert;
    unsigned short nface;
    vcg::Sphere3f sphere;
    float error;
    unsigned short ram;
    unsigned short agp;
  };

  struct Update {
    std::vector<unsigned int> erased;
    std::vector<unsigned int> created;
  };


  Nexus();
  virtual ~Nexus();
  bool Create(const std::string &filename, Signature signature);
  virtual bool Load(const std::string &filename);
  virtual void Close();

  Patch GetPatch(unsigned int patch);
  Border GetBorder(unsigned int patch);

  bool IsCompressed() { return signature & NXS_COMPRESSED; }

  //MOVE to nexus_build.cpp

  unsigned int AddPatch(unsigned int nvert, unsigned int nface, 
			unsigned int nbord);

  //  unsigned int Join(std::vector<unsigned int> &patches);
  void Join(const std::set<unsigned int> &patches,
	    std::vector<vcg::Point3f> &vert,
	    std::vector<unsigned int> &faces,
	    std::vector<Link> &links);

  void Unify(float threshold = 0.0f);

  //TODO implement theese
  void CompactBorder(unsigned int patch);
  void CompactBorders();
  void CompactPatches();
  
  
  Signature signature;
  
  unsigned int totvert;
  unsigned int totface;
  vcg::Sphere3f sphere;
    
  std::vector<Entry> index;

  VFile<Chunk> patches;
  VFile<Link> borders; 
  
  std::vector<Update> history;
 private:
  FILE *index_file;
};

}

#endif
