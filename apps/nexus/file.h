#ifndef NXS_FILE_H
#define NXS_FILE_H

//TODO move includes in cpp

#ifdef WIN32
#ifndef _WINDOWS_
#define _WINSOCKAPI_
#include <windows.h>
#endif
#else
#include <unistd.h>
#endif

#include <stdio.h>
#include <string>

namespace nxs {

class File {
 public:

  File(): fp(NULL) {}
  ~File() { Close(); }   

  File(const File &file) { 
    fp = file.fp;
    size = file.size;
    readonly = file.readonly;
    ((File &)file).fp = NULL;
  }  
  bool Create(const std::string &filename);
  bool Load(const std::string &filename, bool readonly = false);
  void Close();
  
  unsigned int Length() { return size; }
  void Redim(unsigned int elem);

  void SetPosition(unsigned int chunk);
  void ReadBuffer(void *data, unsigned int size);
  void WriteBuffer(void *data, unsigned int size);

  bool IsReadOnly() { return readonly; }

  static void Delete(const std::string &filename);

#ifdef WIN32
   HANDLE fp;
#else
   FILE *fp;  
#endif
   unsigned int size;
   bool readonly;
};

}

#endif
