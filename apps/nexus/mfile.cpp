#include "mfile.h"
#include <assert.h>
#include <iostream>

using namespace std;
using namespace nxs;

bool MFile::Create(const string &fname, unsigned int mxs) {
  filename = fname;
  files.clear();
  size = 0;
  readonly = false;
  assert(mxs <= (1<<30));
  max_size = mxs;
  return AddFile();  
}

bool MFile::Load(const string &fname, bool ronly) {
  filename = fname;
  files.clear();
  readonly = ronly;
  max_size = (1<<30);
  size = 0;
  
  while(1) {
    string name = Name(files.size());
    files.push_back(File());
    File &file = files.back();
    if(!file.Load(name, ronly)) {
      files.pop_back();
      break;
    }
    size += file.Length();
  }
  if(files.size() == 0) return false;
  if(files.size() == 1) { 
    assert(size <= max_size);
  } else {
    //SANITY TEST
    for(unsigned int i = 0; i < files.size() -2; i++) {
      if(files[i].Length() != files[i++].Length()) {
	//"Inconsistent file size for some file.\n";
	return false;
      }
      max_size = files[0].Length();
    }
  }
  return true;
}

void MFile::Close() {
  files.clear();
}

void MFile::Redim(int64 sz) {
  assert(!readonly);
  if(sz > size) {
    unsigned int totfile = (unsigned int)(sz/max_size);
    //TODO test rhis!!!!
    while(files.size() < totfile) {
      RedimLast(max_size);
      AddFile();
    }
    assert(size < sz);
    assert(sz - size < max_size);
    RedimLast(files.back().Length() + (unsigned int)(sz - size));
  } else {
    while(size - files.back().Length() > sz)
      RemoveFile();
    assert(sz <= size);
    RedimLast(files.back().Length() - (unsigned int)(size - sz));
  }    
}

void MFile::SetPosition(int64 pos) {
  assert(pos < size);
  curr_fp = (unsigned int)(pos/(int64)max_size);
  curr_pos = (unsigned int)(pos - (int64)max_size * (int64)curr_fp);
  assert(curr_fp < files.size());
  files[curr_fp].SetPosition(curr_pos);
}

void MFile::ReadBuffer(void *data, unsigned int sz) {
  while(sz + curr_pos > max_size) {
    unsigned int n = max_size - curr_pos;
    files[curr_fp].ReadBuffer(data, n);
    data = ((char *)data) + n;
    sz -= n;
    curr_fp++;
    curr_pos = 0;
    files[curr_fp].SetPosition(curr_pos);
  }
  files[curr_fp].ReadBuffer(data, sz);
}

void MFile::WriteBuffer(void *data, unsigned int sz) {
  assert(!readonly);
  while(sz + curr_pos > max_size) {
    unsigned int n = max_size - curr_pos;
    files[curr_fp].WriteBuffer(data, n);
    data = ((char *)data) + n;
    sz -= n;
    curr_fp++;
    curr_pos = 0;
    files[curr_fp].SetPosition(curr_pos);
  }
  files[curr_fp].WriteBuffer(data, sz);
}

 bool MFile::AddFile() {
   string name = Name(files.size());
   files.push_back(File());
   File &file = files.back();
   return file.Create(name);     
 }

 void MFile::RemoveFile() {
   assert(files.size());

   string name = Name(files.size()); 
   File &file = files.back();
   unsigned int last_size = file.Length();
   files.pop_back();
   size -= last_size;
#ifdef WIN32
   DeleteFile(name.c_str());
#else
   unlink(name.c_str());
#endif
 }

void MFile::RedimLast(unsigned int sz) {
  File &file = files.back();
  unsigned int last_size = file.Length();
  file.Redim(sz);
  size += sz - last_size;
}

std::string MFile::Name(unsigned int n) {
  char buffer[1024];
  if(n == 0) 
    sprintf(buffer, "%s", filename.c_str());
  else
    sprintf(buffer, "%s%d", filename.c_str(), n);
  return string(buffer);
}

