#ifndef NXS_PATCH_SERVER_H
#define NXS_PATCH_SERVER_H

#include "patch.h"
#include "file.h"

#include <vector>

namespace nxs {

struct PatchEntry { 
  unsigned int patch_start;  //granularita' Chunk
  unsigned short ram_size;  //in chunks 
  unsigned short disk_size;  // in chunks (used when compressed)
  unsigned int lru_pos;
};

class PatchServer: public File {
 public:

  struct PTime {
    unsigned int npatch;
    unsigned int frame;
    
    Patch *patch;
    unsigned int vbo_array;
    unsigned int vbo_element;
    bool locked;

    PTime(unsigned int p = 0xffffffff, unsigned int f = 0xffffffff):
	 npatch(p), frame(f), patch(NULL), 
	 vbo_array(0), vbo_element(0) {}

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
    

  PatchServer(): chunk_size(1024), ram_size(128000000), vbo_size(32000000) {}
  bool Create(const std::string &filename, Signature signature, 
	      unsigned int chunk_size, unsigned int ram_size = 0);
  bool Load(const std::string &filename, Signature sig, 
	    unsigned int chunk_size, bool readonly, 
	    unsigned int ram_size = 0);

  void Close();

  bool ReadEntries(FILE *fp);
  bool WriteEntries(FILE *fp);

  void AddPatch(unsigned short nvert, unsigned short nface);
  Patch &GetPatch(unsigned int patch, 
		  unsigned short nvert, unsigned short nface,
		  bool flush = true);

  void GetVbo(unsigned int patch, unsigned int &element, unsigned int &array);

  void Flush(PTime &ptime);
  //return false if was not allocated.
  bool FlushVbo(PTime &ptime);
  void Flush();
  void FlushAll();
  
  void SetRamBufferSize(unsigned int ram_buffer);

  std::vector<PatchEntry> patches;
  std::vector<PTime> lru;
};

}

#endif
