#ifndef NXS_NEXUS_BASE_H
#define NXS_NEXUS_BASE_H

#include <string>
#include <vector>

#include "pserver.h"

namespace nxs {

struct PatchInfo {
  unsigned short nvert;
  unsigned short nface;
  
  vcg::Sphere3f sphere;
  float error;
};


class NexusBase {
 public:

  //TODO optimize to be vector with offset.
  struct Update {
    std::vector<unsigned int> erased;
    std::vector<unsigned int> created;
  };

  NexusBase(): index_file(NULL) {}

  //  bool Create(const std::string &filename, Signature signature,
  //	      unsigned int chunk_size = 1024);
  //  bool Load(const std::string &filename, bool readonly = false);
  //  void Close();

  bool IsCompressed()    { return (signature & NXS_COMPRESSED) != 0; }
  bool HasStrips()       { return (signature & NXS_STRIP) != 0; }
  bool HasColors()       { return (signature & NXS_COLORS) != 0; }
  bool HasNormalsShort() { return (signature & NXS_NORMALS_SHORT) != 0; }
  bool HasNormalsFloat() { return (signature & NXS_NORMALS_FLOAT) != 0; }

  
  //BE CAREFUL: this 2 members get replicated into patchserver
  //TODO fix this nasty thing it is dangerous as it is.
  Signature signature;
  unsigned int chunk_size;
  
  unsigned int totvert;
  unsigned int totface;
  vcg::Sphere3f sphere;
    
  std::vector<PatchInfo> index;
  std::vector<Update> history;

  bool readonly;

 protected:
  FILE *index_file;
};

}

#endif
