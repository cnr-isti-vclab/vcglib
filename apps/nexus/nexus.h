#ifndef NXS_NEXUS_H
#define NXS_NEXUS_H

#include <string>
#include <vector>

#include "nexusbase.h"
#include "lrupserver.h"
#include "borderserver.h"

namespace nxs {

class Nexus: public NexusBase {
 public:

  Nexus() {}
  ~Nexus();
  
  bool Create(const std::string &filename, Signature signature,
	      unsigned int chunk_size = 1024);
  bool Load(const std::string &filename, bool readonly = false);
  void Close();

  unsigned int AddPatch(unsigned int nv, unsigned int nf, unsigned int nb);
  Patch &GetPatch(unsigned int patch, bool flush = true);
  Border GetBorder(unsigned int patch, bool flush = true);

  void AddBorder(unsigned int patch, Link &link);

  bool IsCompressed()    { return (signature & NXS_COMPRESSED) != 0; }
  bool HasStrips()       { return (signature & NXS_STRIP) != 0; }
  bool HasColors()       { return (signature & NXS_COLORS) != 0; }
  bool HasNormalsShort() { return (signature & NXS_NORMALS_SHORT) != 0; }
  bool HasNormalsFloat() { return (signature & NXS_NORMALS_FLOAT) != 0; }

  void MaxRamBuffer(unsigned int ram_size);

  //move to nxsalgo!
  void Unify(float threshold = 0.0f);


  /* Nexus data */
  
  //BE CAREFUL: this 2 members get replicated into patchserver
  //TODO fix this nasty thing it is dangerous as it is.
  //  Signature signature;
  //  unsigned int chunk_size;
  
  //  unsigned int totvert;
  //  unsigned int totface;
  //  vcg::Sphere3f sphere;
    
  //  std::vector<PatchInfo> index;

  LruPServer patches;
  BorderServer borders; 
  
  //  std::vector<Update> history;

  //  bool readonly;

  // private:
  //  FILE *index_file;
};

}

#endif
