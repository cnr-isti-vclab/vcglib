#ifndef NXS_NEXUS_MT_H
#define NXS_NEXUS_MT_H


#include "nexus.h"

namespace nxs {

class NexusMt: public Nexus {
 private:

  class Frag:public std::vector<unsigned int> {};

  struct Node {  
    std::vector<Node *> in;
    std::vector<Node *> out;
    std::vector<Frag> frags;    
    float error;
    bool visited;        
  };

  std::vector<Node> nodes;
 public:
  void LoadHistory();
  void ClearHistory();

  void ExtractFixed(std::vector<unsigned int> &selected, float error);

};

}

#endif
