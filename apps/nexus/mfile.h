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
  
  int64 Length() { return _size; }
  void Redim(int64 size);

  void SetPosition(int64 pos);
  void ReadBuffer(void *data, unsigned int size);
  void WriteBuffer(void *data, unsigned int size);

  bool Opened() { return files.size() > 0; }
  bool IsReadOnly() { return readonly; }
  void SetReadOnly(bool rd) { readonly = rd; } //USE WITH CARE!!!!
 protected:
  std::string filename;
  std::vector<File *> files;
  unsigned int curr_pos;
  unsigned int curr_fp;
  int64 _size;
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
