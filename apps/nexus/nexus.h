#ifndef NXS_NEXUS_H
#define NXS_NEXUS_H

#include <string>
#include <vector>
#include <set>
#include <vcg/space/point3.h>

#include "patchserver.h"
#include "borderserver.h"

namespace nxs {


struct PatchInfo {
  unsigned short nvert;
  unsigned short nface;
  
  vcg::Sphere3f sphere;
  float error;
};



class Nexus {
 public:

  struct BorderEntry {
    unsigned int border_start; //granuralita' Link
    unsigned short border_size; //in Links
    unsigned short border_used; //in Links
  };
  
  struct PatchInfo {
    unsigned short nvert;
    unsigned short nface;
    
    vcg::Sphere3f sphere;
    float error;
  };


  //TODO optimize to be vector with offset.
  struct Update {
    std::vector<unsigned int> erased;
    std::vector<unsigned int> created;
  };


  Nexus();
  virtual ~Nexus();

  bool Create(const std::string &filename, Signature signature,
	      unsigned int chunk_size = 1024);
  virtual bool Load(const std::string &filename, bool readonly = false);
  virtual void Close();

  unsigned int AddPatch(unsigned int nv, unsigned int nf, unsigned int nb);
  Patch &GetPatch(unsigned int patch, bool flush = true);
  Border GetBorder(unsigned int patch, bool flush = true);

  bool IsCompressed()    { return (signature & NXS_COMPRESSED) != 0; }
  bool HasStrips()       { return (signature & NXS_STRIP) != 0; }
  bool HasColors()       { return (signature & NXS_COLORS) != 0; }
  bool HasNormalsShort() { return (signature & NXS_NORMALS_SHORT) != 0; }
  bool HasNormalsFloat() { return (signature & NXS_NORMALS_FLOAT) != 0; }

  void SetRamBufferSize(unsigned int ram_size);

  /*  void Join(const std::set<unsigned int> &patches,
	    std::vector<vcg::Point3f> &vert,
	    std::vector<unsigned int> &faces,
	    std::vector<Link> &links);*/

  //move to nxsalgo!
  void Unify(float threshold = 0.0f);


  /* Nexus data */
  
  //BE CAREFUL: this 2 members get replicated into patchserver
  //TODO fix this nasty thing it is dangerous as it is.
  Signature signature;
  unsigned int chunk_size;
  
  unsigned int totvert;
  unsigned int totface;
  vcg::Sphere3f sphere;
    
  std::vector<PatchInfo> index;

  PatchServer patches;
  BorderServer borders; 
  
  std::vector<Update> history;

  bool readonly;

 private:
  FILE *index_file;
};

}

#endif
