#ifndef NXS_REMAPPING_H
#define NXS_REMAPPING_H

#include <vector>
#include <string>
#include "nxstypes.h"
#include "vchain.h"
#include "nexus.h"
#include "vfile.h"

namespace nxs {
  
  struct BlockEntry {
    BlockEntry(int64 o = 0, unsigned int s = 0): offset(o), size(s) {}
    int64 offset;
    unsigned int size;
  };
  
  class BlockIndex: public std::vector<nxs::BlockEntry> {
  public:
    bool Save(const std::string &file);
    bool Load(const std::string &file);
  };

  void Remap(VChain &chain,
	     VFile<vcg::Point3f> &points,
	     VFile<unsigned int> &remap,
	     BlockIndex &index,
	     unsigned int target_size,
	     unsigned int min_size,
	     unsigned int max_size,
	     float scaling,
	     int step);
  
  void BuildPartition(VPartition &part,
		      VFile<vcg::Point3f> &points,
		      unsigned int target_size,
		      unsigned int min_size,
		      unsigned int max_size,
		      int steps);

  void BuildLevel(VChain &chain,
		  Nexus &nexus,
		  unsigned int offset, 
		  float scaling,
		   unsigned int target_size,
		  unsigned int min_size,
		  unsigned int max_size,
		  int steps);

  //removes small or really big patches.
  bool Optimize(VPartition &part, 
		unsigned int target_cells,
		unsigned int target_size,
		unsigned int min_size,
		unsigned int max_size,
		std::vector<vcg::Point3f> &centroids,
		std::vector<unsigned int> &counts, 
		bool join);

  int GetBest(VPartition &part, unsigned int seed,
		       std::vector<bool> &mark,
		       std::vector<unsigned int> &counts);
}

#endif
