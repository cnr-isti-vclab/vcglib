#include "file.h"
#include <assert.h>

using namespace std;
using namespace nxs;

bool File::Create(const string &filename) {
  size = 0;
#ifdef WIN32
  fp = CreateFile(filename.c_str(), GENERIC_READ | GENERIC_WRITE, 0, 
		  NULL, CREATE_ALWAYS, 0, NULL);
  if(fp == INVALID_HANDLE_VALUE) return false;
#else
  fp = fopen(filename.c_str(), "wb+");
  if(!fp) return false;        
#endif
    return true;
}

bool File::Load(const string &filename) {
#ifdef WIN32
  fp = CreateFile(filename.c_str(), GENERIC_READ | GENERIC_WRITE, 
		  0, NULL, OPEN_EXISTING, 0, NULL); 
  if(fp == INVALID_HANDLE_VALUE) return false;
#else
  fp = fopen(filename.c_str(), "rb+");
  if(!fp) return false;
#endif

#ifdef WIN32    
  size = GetFileSize(fp, NULL);
#else
  //TODO use stat()
  fseek(fp, 0, SEEK_END);
  size = ftell(fp);   
#endif
  return true;
}

void File::Close() {
  if(fp) {
#ifdef WIN32
    CloseHandle(fp);
#else
    fclose(fp);
#endif
    fp = NULL;
  }
}

void File::Resize(unsigned int elem) {
  assert(fp);
  if(elem > size) {
    
#ifdef WIN32
    if(INVALID_SET_FILE_POINTER == 
       SetFilePointer(fp, elem - 1, 0, FILE_BEGIN))
#else
    if(-1 == fseek(fp, elem - 1, SEEK_SET)) 
#endif
      assert(0 && "Could not resize");
    
    unsigned char a;
#ifdef WIN32
    DWORD tmp;
    WriteFile(fp, &a, 1, &tmp, NULL);
#else
    fwrite(&a, sizeof(unsigned char), 1, fp);
#endif
  } else {
    //TODO optimize: we do not need flush for buffers over elem.
#ifndef WIN32
    int fd = fileno(fp);
    ftruncate(fd, elem);
#else
    SetFilePointer(fp, elem, 0, FILE_BEGIN);
    SetEndOfFile(fp);
#endif
  }    
  size = elem;
}

void File::SetPosition(unsigned int pos) {
#ifdef WIN32
  SetFilePointer(fp, pos, 0, FILE_BEGIN);
#else
  fseek(fp, pos, SEEK_SET);
#endif
}

void File::ReadBuffer(void *data, unsigned int sz) {
#ifdef WIN32
  DWORD tmp;
  ReadFile(fp, data, sz, &tmp, NULL);
  if(tmp != sz)
    assert(0 && "Could not read");        
#else    
  if(sz != fread(data, 1, sz, fp)) 
    assert(0 && "Could not read");    
#endif
}

void File::WriteBuffer(void *data, unsigned int sz) {
#ifdef WIN32
  DWORD tmp;
  WriteFile(fp, data, sz, &tmp, NULL);
  assert(tmp == sz);
#else    
  if(sz != fwrite(data, 1, sz, fp)) 
    assert(0 && "Could not write");    
#endif
}
