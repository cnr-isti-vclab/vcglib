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
Revision 1.1  2004/06/24 14:32:45  ponchio
Moved from wrap/nexus

Revision 1.1  2004/06/24 14:18:58  ponchio
Created


****************************************************************************/

#ifndef NXS_MFHASH_H
#define NXS_MFHASH_H

#include "vfile.h"
#include <stdio.h>

namespace nxs {



class MFHash {
 public:
  struct Bucket {
    unsigned int key;
    unsigned int value;
    Bucket(): key(0xffffffff) {}
    Bucket(unsigned int k, unsigned int v): key(k), value(v) {}
    bool Empty() { return key == 0xffffffff; }
  };
  
  MFHash() {}
  bool Create(const std::string &file, unsigned int reserved = 32);
  bool Load(const std::string &file, unsigned int used = 0xffffffff);

  void Resize(unsigned int n);
  void Insert(unsigned int key, unsigned int value, bool rehash = true);
  template <class C> void GetValues(unsigned int key, C &container) {
    container.clear();
    unsigned int hash_size = buffer.Size();
    unsigned int j = key % hash_size;
    while(!buffer[j].Empty()) {
      if(buffer[j].key == key) {
	container.push_back(buffer[j].value);
      }
      j++;
      if(j >= hash_size) j = 0;
    } 
  }
  unsigned int Count(unsigned int key);
  void Clear();
  unsigned int Size();
  void Close();
 private:
  VFile<Bucket> buffer;
  unsigned int space;
    };
}//namespace
#endif
