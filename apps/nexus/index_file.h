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

#ifndef NXS_INDEX_FILE_H
#define NXS_INDEX_FILE_H

#include <assert.h>

#include <vector>

#include "mfile.h"

namespace nxs {

  /* WARINING when subclassing this class you must add a Close() 
     in the destructor! */

template <class T> 
class IndexFile: public MFile, public std::vector<T> {

 public:
  virtual ~IndexFile() {}

  bool Create(const std::string &filename, unsigned int header_size, 
	      unsigned int max_file_size = MFILE_MAX_SIZE) {
    clear();
    if(!MFile::Create(filename, max_file_size)) return false;
    MFile::Redim(header_size);
    offset = header_size;
    return true;
  }

  bool Load(const std::string &filename, bool readonly = false) {
    clear();
    if(!MFile::Load(filename, readonly)) return false;
    SetPosition(0);
    LoadHeader();

    SetPosition(offset);
    unsigned int tot;
    ReadBuffer(&tot, sizeof(unsigned int));
    resize(tot);
    ReadBuffer(&*begin(), size() * sizeof(T));
    return true;
  }

  void Close() {
    if(IsReadOnly()) return;
    if(files.size() == 0) return; //not loaded, not created or closed

    MFile::Redim(offset + size() * sizeof(T));
    SetPosition(offset);
    unsigned int tot = size();
    WriteBuffer(&tot, sizeof(unsigned int));
    WriteBuffer(&*begin(), size() * sizeof(T));
    SetPosition(0);
    SaveHeader();
    MFile::Close();
  }

  int64 Length() { //count the header but not the index...
    return offset;
  }

  void Redim(int64 size) {
    MFile::Redim(size);
    offset = size;
    assert(MFile::Length() == offset);
  }

 protected:
  int64 offset;
  
  //MUST set offset to its correct value
  virtual bool LoadHeader() = 0;
  //MUST save offset somewhere
  virtual void SaveHeader() = 0;
};

} //namespace

#endif
