#ifndef NXS_PATCH_SERVER_H
#define NXS_PATCH_SERVER_H

#include "patch.h"
#include "file.h"

#include <vector>

namespace nxs {

struct PatchEntry { 
  Patch *patch;
  unsigned int patch_start;  //granularita' Chunk
  unsigned short ram_size;  //in chunks 
  unsigned short disk_size;  // in chunks (used when compressed)
  unsigned int lru_pos;
};

 struct VboBuffer {
   VboBuffer(unsigned int v = 0, unsigned int i = 0):
     vertex(v), index(i) {}
   unsigned int vertex;
   unsigned int index;
 };
 

class PatchServer: public File {
 public:
  struct PTime {
    unsigned int patch;
    unsigned int frame;

    PTime(unsigned int p = 0xffffffff, unsigned int f = 0xffffffff):
	 patch(p), frame(f) {}

    bool operator<(const PTime &p) const { return frame > p.frame; }
  };


  Signature signature;
  unsigned int chunk_size;

  unsigned int ram_size;
  unsigned int ram_used;
  unsigned int vbo_size;
  unsigned int vbo_used;
  unsigned int frame;

  //statistics:
  unsigned int ram_readed;
  unsigned int ram_flushed;
    

  bool Create(const std::string &filename, Signature signature, 
	      unsigned int chunk_size, unsigned int ram_size = 128000);
  bool Load(const std::string &filename, Signature sig, 
	    unsigned int chunk_size, bool readonly, 
	    unsigned int ram_size = 128000);

  void Close();

  bool ReadEntries(FILE *fp);
  bool WriteEntries(FILE *fp);

  void AddPatch(unsigned short nvert, unsigned short nface);
  Patch &GetPatch(unsigned int patch, 
		  unsigned short nvert, unsigned short nface,
		  bool flush = true);

  VboBuffer &GetVbo(unsigned int patch);

  void Flush(unsigned int patch);
  //return false if was not allocated.
  bool FlushVbo(unsigned int patch);
  void Flush();
  void FlushAll();
  
  void SetRamBufferSize(unsigned int ram_buffer);

  std::vector<PatchEntry> patches;
  std::vector<VboBuffer> vbos;
  std::vector<PTime> lru;
};

}

#endif
