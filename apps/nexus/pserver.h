#ifndef NXS_PSERVER_H
#define NXS_PSERVER_H

#include <vector>
#include <list>
#include <map>
#include <iostream>

#include "patch.h"
#include "mfile.h"

namespace nxs { 

  /* HEader fo pserver
  1Kb riservato per dati globali:
     Magic:        'n' 'x' 's' 0x00
     Signature:    unsigned int (maschera di bit)
     Chunk size:   unsigned int
     Index offset: unsigned int (offset to the index begin, 
                                 must be a multiple of chunk size)
     Index size:   unsigned int (in number of entryies)
     History mode: unsigned int 0 means erased and created 
                                1 means ready for drawing
     History offset: unsigned int: multiple of chunk_size
     History size: unsigned int (in chars)

     Tot vert:     unsigned int
     Tot face:     unsigned int  
     Bound sphere: Sphere3f (4 float: Point3f center (x, y, z), (radius))*/



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
 
 class PServer: public IndexFile<Entry, 
public:

  std::list<unsigned int> pqueue;
  std::map<unsigned int, std::list<unsigned int>::iterator> index;
  
  Signature signature;
  unsigned int chunk_size;

  unsigned int ram_max;
  unsigned int ram_used;
  
  PServer(): chunk_size(1024), 
    ram_max(128000000), 
    ram_used(0) {}
  
  ~PServer() { PServer::Close(); }

  bool Create(const std::string &filename, Signature signature, 
	      unsigned int chunk_size = 1024);
  bool Load(const std::string &filename, Signature signature, 
	    bool readonly = true, unsigned int chunk_size = 1024);
  void Close();

  bool ReadEntries(FILE *fp);
  bool WriteEntries(FILE *fp);


  void AddPatch(unsigned short nvert, unsigned short nface);
  Patch *LoadPatch(unsigned int id);
  void FlushPatch(unsigned int id);
  void Flush();    

  bool IsLoaded(unsigned int patch);
  Entry &Lookup(unsigned int patch, std::vector<unsigned int> &flushed);
  void LoadVbo(unsigned int patch);
  void FlushVbo(unsigned int patch);

  bool IsCompressed()    { return (signature & NXS_COMPRESSED) != 0; }
  bool HasStrips()       { return (signature & NXS_STRIP) != 0; }
  bool HasColors()       { return (signature & NXS_COLORS) != 0; }
  bool HasNormalsShort() { return (signature & NXS_NORMALS_SHORT) != 0; }
  bool HasNormalsFloat() { return (signature & NXS_NORMALS_FLOAT) != 0; }

  void MaxRamBuffer(unsigned int ram_buffer);
};


}//namespace
#endif
