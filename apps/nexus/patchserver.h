#ifndef NXS_PATCH_SERVER_H
#define NXS_PATCH_SERVER_H

#include "patch.h"
#include "file.h"

#include <vector>

namespace nxs {

struct PatchEntry { 
  Patch *patch;
  unsigned int patch_start;  //granularita' Chunk
  unsigned short patch_size;  //in chunks
  unsigned short ram_used;  // in chunks (used when compressed)
  unsigned int lru_pos;
};

class PatchServer: public File {
 public:
  struct PTime {
    unsigned int patch;
    unsigned int frame;

    PTime(unsigned int p = 0xffffffff, unsigned int f = 0xffffffff):
	 patch(p), frame(f) {}

    bool operator<(const PTime &p) const { return frame < p.frame; }
  };


  Signature signature;
  unsigned int chunk_size;
  unsigned int ram_size;
  unsigned int ram_used;
  unsigned int frame;
    

  bool Create(const std::string &filename, Signature signature, 
	      unsigned int chunk_size, unsigned int ram_size = 128000);
  bool Load(const std::string &filename, Signature sig, 
	    unsigned int chunk_size, bool readonly, 
	    unsigned int ram_size = 128000);

  void Close();

  bool ReadEntries(FILE *fp);
  bool WriteEntries(FILE *fp);

  void AddPatch(unsigned short nvert, unsigned short nface);
  Patch &GetPatch(unsigned int patch, unsigned short nvert, unsigned short nface,
		  bool flush = true);

  void Flush();
  void FlushAll();
  void Flush(unsigned int patch);

  std::vector<PatchEntry> patches;
  std::vector<PTime> lru;
};

}

#endif
