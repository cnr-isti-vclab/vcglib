#ifndef NXS_PSERVER_H
#define NXS_PSERVER_H

#include <vector>
#include <map>
#include <iostream>
#include "patch.h"
#include "mfile.h"

namespace nxs {

class PServer: public MFile {
 public:
  
  struct Entry { 
    unsigned int patch_start;  //granularita' Chunk
    unsigned short ram_size;  //in chunks 
    unsigned short disk_size;  // in chunks (used when compressed)
  };

  /*  struct Data {
    //    unsigned int npatch;
    Patch *patch;
    unsigned int vbo_array;
    unsigned int vbo_element;
    Data(): patch(NULL), vbo_array(0), vbo_element(0) {}
    };*/

  struct Item { //used by lru and pqueue.
    unsigned int patch;
    float priority;
    Item(unsigned int p, float f): patch(p), priority(f) {}
    bool operator<(const Item &i) const {
      return priority < i.priority;
    }
  };
  


  Signature signature;
  unsigned int chunk_size;

  unsigned int ram_max;
  unsigned int ram_used;

  //statistics:
  unsigned int ram_readed;
  unsigned int ram_flushed;
    
  //pt::rwlock ramlock; //read only thread safety...
  //pt::rwlock disklock; //read only thread safety...

  std::vector<Entry> entries;
  
  
  
  PServer(): chunk_size(1024), 
             ram_max(128000000), 
             ram_used(0) {}
  virtual ~PServer() {
    std::cerr << "Closing pserver" << std::endl;
    MFile::Close();
  }

  bool Create(const std::string &filename, Signature signature, 
	      unsigned int chunk_size, unsigned int ram_max = 128000000);
  bool Load(const std::string &filename, Signature sig, 
	    unsigned int chunk_size, bool readonly, 
	    unsigned int ram_max = 128000000);

  bool ReadEntries(FILE *fp);
  bool WriteEntries(FILE *fp);
  virtual void Close();

  void AddPatch(unsigned short nvert, unsigned short nface);
  Patch *LoadPatch(unsigned int id, unsigned short nv, unsigned short nf);
  void FlushPatch(unsigned int id, Patch *patch);
  
  virtual bool IsLoaded(unsigned int patch) = 0;
  //  virtual Patch &Lookup(unsigned int patch, 
  //		      unsigned short nv, unsigned short nf) = 0;
  virtual void Flush() = 0;

  void MaxRamBuffer(unsigned int ram_buffer);

  //void GetVbo(unsigned int patch, unsigned int &elem, unsigned int &array);  
  //bool FlushVbo(unsigned int patch, Data &data);
  //void MaxVboBuffer(unsigned int ram_buffer);
};


}//namespace
#endif
