#ifndef NXS_NEXUS_H
#define NXS_NEXUS_H

#include <string>
#include <vector>
#include <list>
#include <map>

#include <vcg/space/sphere3.h>

#include "patch.h"
#include "index_file.h"
#include "history.h"
#include "borderserver.h"

namespace nxs {

  /* Header fo nexus: 
  1Kb riservato per dati globali:
     Magic:        'n' 'x' 's' 0x00
     Signature:    unsigned int (maschera di bit)
     Chunk size:   unsigned int
     Index offset: unsigned int (offset to the index begin, 
                                 must be a multiple of chunk size)
     History offset: unsigned int: multiple of chunk_size

     Tot vert:     unsigned int
     Tot face:     unsigned int  
     Bound sphere: Sphere3f (4 float: Point3f center (x, y, z), (radius))

     11 * 4 = 44 bytes -> 4k per alignment purpoouses and reserving space. */

struct Entry { 
  unsigned int patch_start;  //granularita' Chunk
  unsigned short ram_size;  //in chunks 
  unsigned short disk_size;  // in chunks (used when compressed)

  unsigned short nvert;
  unsigned short nface;
  
  vcg::Sphere3f sphere;
  float error;

  Patch *patch;
  unsigned int vbo_array;
  unsigned int vbo_element;
};

class Nexus: public IndexFile<Entry> {
 public:
  //HEader data:
  Signature signature;
  unsigned int chunk_size;
  //unsigned int .IndexFile::offset;
  int64 history_offset;
  
  unsigned int totvert;
  unsigned int totface;
  vcg::Sphere3f sphere;
  
  History history;

  BorderServer borders;   
  

  Nexus() {}
  ~Nexus();
  
  bool Create(const std::string &filename, Signature signature,
	      unsigned int chunk_size = 1024);
  bool Load(const std::string &filename, bool readonly = false);
  void Close();
  void Flush(bool all = true);
  
  unsigned int AddPatch(unsigned int nv, unsigned int nf, unsigned int nb);
  Patch &GetPatch(unsigned int patch, bool flush = true);
  Border GetBorder(unsigned int patch, bool flush = true);

  unsigned int &MaxRam() { return ram_max; }
  //  void AddBorder(unsigned int patch, Link &link);

  //move to nxsalgo!
  void Unify(float threshold = 0.0f);

  bool IsCompressed()    { return (signature & NXS_COMPRESSED) != 0; }
  bool HasStrips()       { return (signature & NXS_STRIP) != 0; }
  bool HasColors()       { return (signature & NXS_COLORS) != 0; }
  bool HasNormalsShort() { return (signature & NXS_NORMALS_SHORT) != 0; }
  bool HasNormalsFloat() { return (signature & NXS_NORMALS_FLOAT) != 0; }

  unsigned int ram_max;
  unsigned int ram_used;		
 protected:

  std::list<unsigned int> pqueue;
  std::map<unsigned int, std::list<unsigned int>::iterator> index;

  Patch *LoadPatch(unsigned int id);
  virtual void FlushPatch(unsigned int id);

  bool LoadHeader();
  void SaveHeader();
};

}

#endif
