/****************************************************************************
* VCGLib                                                            o o     *
* Visual and Computer Graphics Library                            o     o   *
*                                                                _   O  _   *
* Copyright(C) 2004                                                \/)\/    *
* Visual Computing Lab                                            /\/|      *
* ISTI - Italian National Research Council                           |      *
*                                                                    \      *
* All rights reserved.                                                      *
*                                                                           *
* This program is free software; you can redistribute it and/or modify      *   
* it under the terms of the GNU General Public License as published by      *
* the Free Software Foundation; either version 2 of the License, or         *
* (at your option) any later version.                                       *
*                                                                           *
* This program is distributed in the hope that it will be useful,           *
* but WITHOUT ANY WARRANTY; without even the implied warranty of            *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
* GNU General Public License (http://www.gnu.org/licenses/gpl.txt)          *
* for more details.                                                         *
*                                                                           *
****************************************************************************/
/****************************************************************************
  History

$Log: not supported by cvs2svn $
Revision 1.7  2005/02/10 09:18:20  ponchio
Statistics.

Revision 1.6  2005/02/08 12:43:03  ponchio
Added copyright


****************************************************************************/

#ifndef NXS_EXTRACTION_H
#define NXS_EXTRACTION_H

#include <set>
#include <vector>

#include <wrap/gui/frustum.h>

#include "history.h"		

namespace nxs {

class Metric;
class NexusMt; 

struct Item {
  float error;
  unsigned int id;
  Item(unsigned int i = 0, float e = 0): id(i), error(e) {}
  bool operator<(const Item &item) const {
    return error < item.error;
  }
};

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

  //TODO make a pointer so no need to doubl check visibility...!
  vcg::Frustumf frustum;

  float target_error;
  float max_error;
  unsigned int extr_used, extr_max;
  unsigned int draw_used, draw_max;
  unsigned int disk_used, disk_max;

  std::vector<bool> visited;
  std::vector<bool> visible;
  std::map<unsigned int, float> errors;

  std::vector<Item> selected;
  unsigned int draw_size; //first in selected should be drawn

  std::vector<HeapNode> heap; //no realtime extraxtion

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

  bool Refine(HeapNode node);
  bool Coarse(HeapNode node);

  void Init();
 private:
  NexusMt *mt;
  Node *root;
  Node *sink;

  bool Visited(Node *node)            { return visited[node - root]; }
  void SetVisited(Node *node, bool v) { visited[node - root] = v; }
  bool Visible(Node *node)            { return visible[node - root]; }
  void SetVisible(Node *node, bool v) { visible[node - root] = v; }

  float GetRefineError(Node *node);
};


}//namespace

#endif
