#ifndef NXS_EXTRACTION_H
#define NXS_EXTRACTION_H

#include "history.h"		
#include <vector>
#include <wrap/gui/frustum.h>

namespace nxs {

class Metric;
class NexusMt; 

class Extraction {
 public:
  struct HeapNode {
    History::Node *node;
    float error;
    unsigned int extr;
    unsigned int draw;
    unsigned int disk;
    HeapNode(History::Node *_node, float _error): node(_node), error(_error),
	 extr(0), draw(0), disk(0) {}
    bool operator<(const HeapNode &node) const {
      return error < node.error; }
  };

  Metric *metric;

  vcg::Frustumf frustum;

  float target_error;
  unsigned int extr_used, extr_max;
  unsigned int draw_used, draw_max;
  unsigned int disk_used, disk_max;

  std::vector<bool> visited;
  std::vector<HeapNode> heap;
  std::vector<unsigned int> selected;

  Extraction();
  ~Extraction();

  void Extract(NexusMt *mt);
  //  void Update(std::vector<unsigned int> &selected);


 protected:
  void Select();
  void Visit(History::Node *node);

  bool Expand(HeapNode &node);
  void Diff(HeapNode &node);

 private:
  NexusMt *mt;
};


}//namespace

#endif
