#ifndef NXS_NEXUS_H
#define NXS_NEXUS_H

#include <string>
#include <vector>
#include <vcg/space/sphere3.h>
#include "vfile.h"
#include "patch.h"
#include "border.h"

namespace nxs {



class Nexus {
 public:

  class Entry {
  public:
    Entry(): patch_start(0xffffffff), border_start(0xffffffff),
      patch_size(0), border_size(0), 
      nvert(0), nface(0), sphere(vcg::Sphere3f()) {}
    unsigned int patch_start;  //granularita' Chunk
    unsigned int border_start; //granuralita' Link
    unsigned short patch_size;  //in cuhnks
    unsigned short border_size; //in Links

    unsigned short nvert;
    unsigned short nface;
    vcg::Sphere3f sphere;
  };

  Nexus();
  ~Nexus();
  bool Create(const std::string &filename);
  bool Load(const std::string &filename);
  void Close();

  Patch GetPatch(unsigned int patch);
  Border GetBorder(unsigned int patch);

  unsigned int AddPatch(unsigned int nvert, unsigned int nface, 
			unsigned int nbord);

  //  unsigned int Join(std::vector<unsigned int> &patches);
  void Join(std::vector<unsigned int> &patches,
	    std::vector<Point3f &vert,
	    std::vector<unsigned int> &faces,
	    std::vector<Link> &links);

  //TODO implement theese
  void CompactBorder(unsigned int patch);
  void CompactBorders();
  void CompactPatches();
  
  unsigned int totvert;
  unsigned int totface;
  unsigned int totchunks; //number of chunks.
  unsigned int totlinks;
  vcg::Sphere3f sphere;
    
  std::vector<Entry> index;

  VFile<Chunk> patches;
  VFile<Link> borders; 
 private:
  FILE *index_file;
};

}

#endif
