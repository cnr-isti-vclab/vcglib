#ifndef NXS_BORDERSERVER_H
#define NXS_BORDERSERVER_H

#include "vfile.h"
#include "border.h"
#include <vector>

namespace nxs {

struct BorderEntry {
  unsigned int border_start; //granuralita' Link
  unsigned short border_size; //in Links
  unsigned short border_used; //in Links
};

class BorderServer: public VFile<Link> {
 public:
  void AddBorder(unsigned short nbord, unsigned int used = 0);
  Border GetBorder(unsigned int border, bool flush = true);
  //return true if you need to reread border as it changed location
  bool ResizeBorder(unsigned int border, unsigned int nbord);

  bool ReadEntries(FILE *fp);
  bool WriteEntries(FILE *fp);

  unsigned int BorderSize(unsigned int i) { 
    return borders[i].border_used; 
  }
  unsigned int BorderCapacity(unsigned int i) { 
    return borders[i].border_size; 
  }
  
  

  
  std::vector<BorderEntry> borders;

};

}

#endif
