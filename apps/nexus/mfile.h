#ifndef NXS_MFILE_H
#define NXS_MFILE_H

#include <string>
#include <vector>

#include "file.h"
#include "nxstypes.h"

namespace nxs {

  /*#ifdef WIN32
typedef __int64 int64;
#else
typedef unsigned long long int64;
#endif*/

#define MFILE_MAX_SIZE (1<<30)

class MFile {
 public:

  MFile() {}
  ~MFile() { Close(); }   
  
  //max is  so default is 1 G
  bool Create(const std::string &filename, 
	      unsigned int max_file_size = MFILE_MAX_SIZE);
  bool Load(const std::string &filename, bool readonly = false);
  void Close();
  void Delete();
  
  int64 Length() { return size; }
  void Redim(int64 size);

  void SetPosition(int64 pos);
  void ReadBuffer(void *data, unsigned int size);
  void WriteBuffer(void *data, unsigned int size);

  bool IsReadOnly() { return readonly; }
  void SetReadOnly(bool rd) { readonly = rd; } //USE WITH CARE!!!!
 protected:
  std::string filename;
  std::vector<File *> files;
  unsigned int curr_pos;
  unsigned int curr_fp;
  int64 size;
  unsigned int max_size;
  bool readonly;
 private:
  //all theese refer to the last in the fp.
  bool AddFile();
  void RemoveFile();
  void RedimLast(unsigned int sz);
  unsigned int GetSize();
  std::string Name(unsigned int n);
};

}

#endif
