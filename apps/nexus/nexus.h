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
    Entry(): patch_offset(0xffffffff), border_offset(0xffffffff),
      patch_size(0), border_size(0), sphere(vcg::Sphere3f()) {}
    unsigned int patch_offset;  //granularita' Chunk
    unsigned int border_offset; //granuralita' Link
    unsigned short patch_size;  //in cuhnks
    unsigned short border_size; //in Links
    vcg::Sphere3f sphere;
  };

  Nexus();
  ~Nexus();
  bool Create(const std::string &filename);
  bool Load(const std::string &filename);
  void Close();

  Patch GetPatch(unsigned int patch);
  void GetBorder(unsigned int border, Border &border);

  //  unsigned int addPatch(Patch *builder);
  void AddBorder(unsigned int patch, std::vector<Link> &links);
  
  unsigned int totvert;
  unsigned int totface;
  unsigned int totchunks; //number of chunks.
  unsigned int totlinks;
  vcg::Sphere3f sphere;
    
  std::vector<Entry> index;

  FILE *index_file;
  VFile<Chunk> patches;
  VFile<Link> borders; 
};

}

#endif
