#ifndef NXS_BORDERSERVER_H
#define NXS_BORDERSERVER_H

#include <list>
#include <map>

#include "index_file.h"
#include "border.h"


/*nell'header ci sta scritto solo:
  spazio riservato: 2 * sizeof(Link);
   magic: nxb0 (4 bytes)
   offset: (int64) */

namespace nxs {

  //This should be Border class!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
struct BorderEntry {
  unsigned int start; //granuralita' Link
  unsigned short size; //in Links             //TODO what are this? be clear!
  unsigned short used; //in Links
  Link *links;
};

class BorderServer: public IndexFile<BorderEntry> {
 public:
  BorderServer(): ram_max(1000000), ram_used(0) {}
  ~BorderServer() { Close(); }
  bool Create(const std::string &file);
  bool Load(const std::string &file, bool readonly = true);
  void Close();
  void Flush();

  void AddBorder(unsigned short nbord, unsigned int used = 0);
  Border GetBorder(unsigned int border, bool flush = true);
  //return true if you need to reread border as it changed location
  bool ResizeBorder(unsigned int border, unsigned int nbord);

  unsigned int BorderSize(unsigned int i) { 
    return operator[](i).used; 
  }
  unsigned int BorderCapacity(unsigned int i) { 
    return operator[](i).size; 
  }
 protected:
  unsigned int ram_max;
  unsigned int ram_used;		
  std::list<unsigned int> pqueue;
  std::map<unsigned int, std::list<unsigned int>::iterator> index;

  bool LoadHeader();
  void SaveHeader();
  
  void FlushBorder(unsigned int border);
  Link *GetRegion(unsigned int start, unsigned int size); //size in links.
};

}

#endif
