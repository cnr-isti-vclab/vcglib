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
Revision 1.3  2004/07/02 17:40:30  ponchio
Debug.

Revision 1.2  2004/07/01 21:34:04  ponchio
Rehash bug.

Revision 1.1  2004/06/24 14:32:45  ponchio
Moved from wrap/nexus

Revision 1.1  2004/06/24 14:18:58  ponchio
Created


****************************************************************************/

#include "mfhash.h"
#include <iostream>

using namespace std;
using namespace nxs;

bool MFHash::Create(const string &file, unsigned int reserved) {
  if(!buffer.Create(file)) return false;
  
  buffer.Resize(reserved);
  Bucket empty;
  for(unsigned int i = 0; i < buffer.Size(); i++)
    buffer[i] = empty;
  space = reserved;
  return true;
}

bool MFHash::Load(const string &file, unsigned int used) {
  if(!buffer.Load(file)) return false;
  if(used != 0xffffffff) {
    space = buffer.Size() - used;
  } else {
    space = 0;
    for(unsigned int i = 0; i < buffer.Size(); i++)
      if(buffer[i].Empty()) space++;
  }
  return true;
}

void MFHash::Resize(unsigned int n) {
  assert(buffer.Size() - space <= n);
  //lets dump actual content
  FILE *fp = tmpfile();
  unsigned int count = 0;
  for(unsigned int i = 0; i < buffer.Size(); i++) {
    Bucket &bucket = buffer[i];
    if(!bucket.Empty()) {
      fwrite(&bucket, sizeof(Bucket), 1, fp);
      ++count;
    }
  }

  buffer.Resize(n);
  Clear();

  rewind(fp);
  Bucket bucket;
  for(unsigned int i = 0; i < count; i++) {
    fread(&bucket, sizeof(Bucket), 1, fp);
    Insert(bucket.key, bucket.value, false);
  }
  fclose(fp);
}

void MFHash::Insert(unsigned int key, unsigned int value, bool rehash) {
  if(buffer.Size() < 5)
    Resize(5);
  assert(space > 0);
  unsigned int hash_size = buffer.Size();
  unsigned int j = key % hash_size;
  Bucket bucket = buffer[j];
  while(!bucket.Empty()) {
    if(bucket.key == key && bucket.value == value)  //already here
      return;
    j++;
    if(j >= hash_size) j = 0;
    bucket = buffer[j];
  } 
  buffer[j] = Bucket(key, value);
  space--;

  if(rehash) {
    float ratio = space / (float)buffer.Size();
    if(ratio < 0.4) { //need to resize 
      Resize(buffer.Size() * 2 + 3);
    }
  }
}

unsigned int MFHash::Count(unsigned int key) {
  unsigned int count = 0;
  unsigned int hash_size = buffer.Size();
  unsigned int j = key % hash_size;

  while(!buffer[j].Empty()) {
    if(buffer[j].key == key)
      count++;
    j++;
    if(j >= hash_size) j = 0;
  } 
  return count;
}

void MFHash::Clear() {
  Bucket empty;
  for(unsigned int i = 0; i < buffer.Size(); i++)
    buffer[i] = empty;
  space = buffer.Size();
}

void MFHash::Close() {
  buffer.Close();
}

unsigned int MFHash::Size() {
  return buffer.Size() - space;
}
