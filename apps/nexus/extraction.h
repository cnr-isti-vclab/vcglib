#ifndef NXS_EXTRACTION_H
#define NXS_EXTRACTION_H

#include <set>
#include <vector>

#include <wrap/gui/frustum.h>

#include "history.h"		

namespace nxs {

class Metric;
class NexusMt; 

class Extraction {
 public:
  typedef History::Node Node;
  typedef History::Link Link;

  struct Cost {
    unsigned int extr;
    unsigned int draw;
    unsigned int disk;
    Cost(): extr(0), draw(0), disk(0) {}
  };
  
  struct HeapNode {
    Node *node;
    float error;

    HeapNode(Node *_node, float _error): node(_node), error(_error) {}
    bool operator<(const HeapNode &node) const {
      return error < node.error; }
    bool operator>(const HeapNode &node) const {
      return error > node.error; }
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
  unsigned int draw_size; //first in selected should be drawn

  //nodes that i can expand to
  std::vector<HeapNode> front;
  //nodes that i can contract
  std::vector<HeapNode> back;

  unsigned int tot_budget;

  Extraction();
  ~Extraction();

  void Extract(NexusMt *mt);
  void Update(NexusMt *mt);


 protected:         

  void Select();
  void Visit(Node *node);

  bool Expand(HeapNode &node);
  void Diff(Node *node, Cost &cost);

  bool Refine(HeapNode &node);
  bool Coarse(HeapNode &node);

  void Init();
 private:
  NexusMt *mt;
  Node *root;
  Node *sink;

  bool Visited(Node *node) {
    return visited[node - root];
  }

  float GetRefineError(Node *node);
};


}//namespace

#endif
