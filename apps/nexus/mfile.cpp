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

****************************************************************************/

#include "mfile.h"
#include <assert.h>
#include <iostream>

using namespace std;
using namespace nxs;

bool MFile::Create(const string &fname, unsigned int mxs) {
  Close();
  filename = fname;
  _size = 0;
  readonly = false;
  assert(mxs <= MFILE_MAX_SIZE);
  max_size = mxs;
  return AddFile();  
}

bool MFile::Load(const string &fname, bool ronly) {
  Close();
  filename = fname;
  readonly = ronly;
  max_size = MFILE_MAX_SIZE;
  _size = 0;
  
  while(1) {
    string name = Name(files.size());
    File *file = new File;
    files.push_back(file);
    if(!file->Load(name, ronly)) {
      files.pop_back();
      break;
    }
    _size += file->Length();
  }
  if(files.size() == 0) return false;
  if(files.size() == 1) { 
    assert(_size <= max_size);
  } else {
    //SANITY TEST
    for(unsigned int i = 0; i < files.size() -2; i++) {
      if(files[i]->Length() != files[i++]->Length()) {
	//"Inconsistent file size for some file.\n";
	return false;
      }
      max_size = files[0]->Length();
    }
  }
  return true;
}

void MFile::Close() {
  for(unsigned int i = 0; i < files.size(); i++)
    delete files[i];
  files.clear();
}

void MFile::Delete() {
  while(files.size())
    RemoveFile();
}

void MFile::Redim(int64 sz) {
  assert(!readonly);
  if(sz > _size) {
    unsigned int totfile = (unsigned int)(sz/max_size);
    //TODO test rhis!!!!
    while(files.size() <= totfile) {
      RedimLast(max_size);
      assert(_size == (int64)max_size * (int64)(files.size()));
      AddFile();
    }
    assert(_size <= sz);
    assert(sz - _size < max_size);
    assert(files.back()->Length() + (unsigned int)(sz - _size) < max_size);
    RedimLast(files.back()->Length() + (unsigned int)(sz - _size));
  } else {
    while(_size - files.back()->Length() > sz)
      RemoveFile();
    assert(sz <= _size);
    RedimLast(files.back()->Length() - (unsigned int)(_size - sz));
  }    
  assert(sz == _size);
}

void MFile::SetPosition(int64 pos) {
  assert(pos <= _size);
  curr_fp = (unsigned int)(pos/(int64)max_size);
  curr_pos = (unsigned int)(pos - (int64)max_size * (int64)curr_fp);  
  assert(curr_pos < max_size);
  assert(curr_fp < files.size());
  files[curr_fp]->SetPosition(curr_pos);
}

void MFile::ReadBuffer(void *data, unsigned int sz) {
  while(sz + curr_pos > max_size) {
    unsigned int n = max_size - curr_pos;
    files[curr_fp]->ReadBuffer(data, n);
    data = ((char *)data) + n;
    sz -= n;
    curr_fp++;
    assert(curr_fp < files.size());
    curr_pos = 0;
    files[curr_fp]->SetPosition(curr_pos);
  }
  files[curr_fp]->ReadBuffer(data, sz);
}

void MFile::WriteBuffer(void *data, unsigned int sz) {
  assert(!readonly);
  while(sz + curr_pos > max_size) {
    unsigned int n = max_size - curr_pos;
    files[curr_fp]->WriteBuffer(data, n);
    data = ((char *)data) + n;
    sz -= n;
    curr_fp++;
    assert(curr_fp < files.size());
    curr_pos = 0;
    files[curr_fp]->SetPosition(curr_pos);
  }
  files[curr_fp]->WriteBuffer(data, sz);
}

 bool MFile::AddFile() {
   string name = Name(files.size());
   File *file = new File;
   files.push_back(file);
   return file->Create(name);     
 }

 void MFile::RemoveFile() {
   assert(files.size());

   string name = Name(files.size()-1); 
   File *file = files.back();
   unsigned int last_size = file->Length();
   delete file;
   files.pop_back();
   _size -= last_size;
   cerr << "Removing file: " << name << endl;
#ifdef WIN32
   DeleteFile(name.c_str());
#else
   unlink(name.c_str());
#endif
 }

void MFile::RedimLast(unsigned int sz) {  
  assert(sz <= max_size);
  File &file = *files.back();
  unsigned int last_size = (int64)file.Length();
  file.Redim(sz);
  _size += sz - (int64)last_size;
}

std::string MFile::Name(unsigned int n) {
  char buffer[1024];
  if(n == 0) 
    sprintf(buffer, "%s", filename.c_str());
  else
    sprintf(buffer, "%s%d", filename.c_str(), n);
  return string(buffer);
}

