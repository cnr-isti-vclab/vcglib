#ifndef NXS_MFILE_H
#define NXS_MFILE_H


#include "file.h"
#include <string>
#include <vector>

namespace nxs {


class MFile {
 public:

  MFile() {}
  ~MFile() { Close(); }   
  
  //max is  so default is 1 G
  bool Create(const std::string &filename, 
	      unsigned int max_file_size = (1<<30));
  bool Load(const std::string &filename, bool readonly = false);
  void Close();
  
  unsigned long long Length() { return size; }
  void Redim(unsigned long long size);

  void SetPosition(unsigned long long pos);
  void ReadBuffer(void *data, unsigned int size);
  void WriteBuffer(void *data, unsigned int size);

  bool IsReadOnly() { return readonly; }

 protected:
  std::string filename;
  std::vector<File> files;
  unsigned int curr_pos;
  unsigned int curr_fp;
  unsigned long long size;
  unsigned int max_size;
  bool readonly;
 private:
  //all theese refer to the last in the fp.
  void AddFile();
  void RemoveFile();
  void RedimLast(unsigned int sz);
  unsigned int GetSize();
  std::string Name(unsigned int n);
};

}

#endif
