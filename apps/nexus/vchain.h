#ifndef NXS_VCHAIN_H
#define NXS_VCHAIN_H

#include <vector>
#include <map>
#include <set>
#include "vpartition.h"

namespace nxs {

class VChain: public std::vector<VPartition *> {
 public:
  ~VChain();
  bool Save(const std::string &file);
  bool Load(const std::string &file);

  std::map<unsigned int, std::set<unsigned int> > newfragments;
  std::map<unsigned int, std::set<unsigned int> > oldfragments;
};

}
#endif
