#ifndef PARTITION_INTERSECT_H
#define PARTITION_INTERSECT_H

#include <map>
#include <vcg/space/point3.h>

namespace nxs {

template <class Partition> class PIntersect {
 public:
  Partition &fine;
  Partition &coarse;

  typedef typename Partition::Key Key;
  
  struct Pair {
    Key fine;
    Key coarse;
    bool operator<(const Pair &p) const {
      if(fine == p.fine) return coarse < p.coarse;
      return fine < p.fine;
    }
  };
  struct Cell {
    unsigned int index;
    unsigned int count;
  };

  unsigned int last_cell;
  std::map<Pair, Cell>  pairs;
  std::set<unsigned int> cells;
  
  PIntersect(Partition &f, Partition &c):
    fine(f), coarse(c), last_cell(0) {}

  unsigned int Locate(const vcg::Point3f &p) {
    Pair pp;
    pp.fine = fine.Locate(p);
    pp.coarse = coarse.Locate(p);
    if(pairs.count(pp)) { //already there
      Cell &cell = pairs[pp];
      cell.count++;
      return cell.index;
    } else {
      Cell &cell = pairs[pp];
      cell.index = last_cell++;
      cell.count = 1;
      cells.insert(cell.index);
      return cell.index;
    }

  }

  unsigned int Relocate(const vcg::Point3f &p) {
    //TODO!!!!!!
  }
  void Prune(unsigned int threshold) {
    //TODO!!!!!!
  }
  unsigned int Count(unsigned int cell) {
    return cells.count(cell);
  }
};

}

#endif
