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
Revision 1.29  2005/02/19 10:45:04  ponchio
Patch generalized and small fixes.

Revision 1.28  2005/02/08 12:43:03  ponchio
Added copyright


****************************************************************************/

#include <assert.h>

#include <iostream>
#include <set>

#include "nexus.h"

using namespace std;
using namespace vcg;
using namespace nxs;

Nexus::~Nexus() {
  Close();
}

bool Nexus::Create(const string &file, Signature &sig, unsigned int c_size) {
  signature = sig;
  totvert = 0;
  totface = 0;
  sphere = Sphere3f();
  chunk_size = c_size;
  unsigned int header_size = 256; //a bit more than 64 needed
  if(chunk_size > header_size) header_size = chunk_size;
  
  history.Clear();
  ram_used = 0;
  ram_max = 50 * (1<<20) / chunk_size;

  if(!IndexFile<Entry>::Create(file + ".nxp", header_size)) {
    cerr << "Could not create file: " << file << ".nxp" << endl;
    return false;
  }

  //Important: chunk_size must be 1 so that i can use Region in VFile.
  if(!borders.Create(file + ".nxb")) {
    cerr << "Could not create file: " << file << ".nxb" << endl;
    return false;
  }
  return true;
}


bool Nexus::Load(const string &file, bool rdonly) {
  if(!IndexFile<Entry>::Load(file + ".nxp", rdonly)) return false;
  ram_used = 0;
  ram_max = 50 * (1<<20) / chunk_size;

  history.Clear();
  SetPosition(history_offset);
  unsigned int history_size;
  ReadBuffer(&history_size, sizeof(unsigned int));

  char *buffer = new char[history_size];
  ReadBuffer(buffer, history_size);

  if(!history.Load(history_size, buffer)) {
    cerr << "Error loading history\n";
    return false;
  }

  borders.Load(file + ".nxb", rdonly);
  //TODO on nxsbuilder assure borders are loaded
  return true;
}

void Nexus::Close() { 
  if(!Opened()) return;

  Flush();

  if(!IsReadOnly()) {
    //set history_offset
    history_offset = 0;
    if(size()) {
      //we need to discover where is the last patch
      for(unsigned int i = 0; i < size(); i++) {
	Entry &e = operator[](i);
	if(e.patch_start + e.disk_size > history_offset)
	  history_offset = e.patch_start + e.disk_size;
      }
      //      history_offset = (back().patch_start + back().disk_size);
    }
    history_offset *= chunk_size;
    
    unsigned int history_size;
    char *mem = history.Save(history_size);
    Redim(history_offset + history_size + sizeof(unsigned int));
    SetPosition(history_offset);
    WriteBuffer(&history_size, sizeof(unsigned int));
    WriteBuffer(mem, history_size);
    delete []mem;
  }
  borders.Close();
  IndexFile<Entry>::Close();
}

void Nexus::SaveHeader() {
  unsigned int magic = 0x3053584e; // nxs0
  WriteBuffer(&magic, sizeof(unsigned int));
  unsigned int version = 1;
  WriteBuffer(&version, sizeof(unsigned int));

  WriteBuffer(&signature, sizeof(Signature));
  WriteBuffer(&chunk_size, sizeof(unsigned int));
  WriteBuffer(&offset, sizeof(int64));
  WriteBuffer(&history_offset, sizeof(int64));
  WriteBuffer(&totvert, sizeof(unsigned int));
  WriteBuffer(&totface, sizeof(unsigned int));
  WriteBuffer(&sphere, sizeof(Sphere3f));
}

bool Nexus::LoadHeader() {
  unsigned int magic;
  ReadBuffer(&magic, sizeof(unsigned int));
  if(magic != 0x3053584e) {
    cerr << "Invalid magic. Not a nxs file\n";
    return false;
  }
  //Current version is 1
  unsigned int version;
  ReadBuffer(&version, sizeof(unsigned int));
  if(version != NXS_CURRENT_VERSION) {
    cerr << "Old version. Sorry.\n";
    return false;
  }
  ReadBuffer(&signature, sizeof(Signature));
  ReadBuffer(&chunk_size, sizeof(unsigned int));
  ReadBuffer(&offset, sizeof(int64));
  ReadBuffer(&history_offset, sizeof(int64));
  ReadBuffer(&totvert, sizeof(unsigned int));
  ReadBuffer(&totface, sizeof(unsigned int));
  ReadBuffer(&sphere, sizeof(Sphere3f));
  return true;
}

void Nexus::Flush(bool all) {
  if(all) {
    std::map<unsigned int, list<unsigned int>::iterator>::iterator i;
    for(i = index.begin(); i != index.end(); i++) {
      unsigned int patch = (*i).first;
      FlushPatch(patch);
    }
    pqueue.clear();
    index.clear();
  } else {
    while(ram_used > ram_max) {    
      unsigned int to_flush = pqueue.back();
      pqueue.pop_back();
      index.erase(to_flush);        
      FlushPatch(to_flush);
    }
  }
}



Patch &Nexus::GetPatch(unsigned int patch, bool flush) { 
  Entry &entry = operator[](patch);
  if(index.count(patch)) {
    assert(entry.patch);
    list<unsigned int>::iterator i = index[patch];
    pqueue.erase(i);
    pqueue.push_front(patch);
    index[patch] = pqueue.begin();
  } else {
    while(flush && ram_used > ram_max) {    
      unsigned int to_flush = pqueue.back();
      pqueue.pop_back();
      index.erase(to_flush);        
      FlushPatch(to_flush);
    }
    assert(!entry.patch);
    entry.patch = LoadPatch(patch);
    pqueue.push_front(patch);
    list<unsigned int>::iterator i = pqueue.begin();
    index[patch] = i;      
  }                        
  return *(entry.patch);
}

Border &Nexus::GetBorder(unsigned int patch, bool flush) {
  return borders.GetBorder(patch);
}

unsigned int Nexus::AddPatch(unsigned int nvert, unsigned int nface,
			     unsigned int nbord) {

  Entry entry;
  entry.patch_start = 0xffffffff;
  entry.ram_size = Patch::ChunkSize(signature, nvert, nface, chunk_size);
  entry.disk_size = 0xffff;
  entry.nvert = nvert;
  entry.nface = nface;
  entry.error = 0;
  //sphere undefined.
  entry.patch = NULL;
  entry.vbo_array = 0;
  entry.vbo_element = 0;
  
  push_back(entry);
  
  borders.AddBorder(nbord);

  totvert += nvert;
  totface += nface;
  return size() - 1;
}

Patch *Nexus::LoadPatch(unsigned int idx) {
  assert(idx < size());
  Entry &entry = operator[](idx);
  if(entry.patch) return entry.patch;
  
  char *ram = new char[entry.ram_size * chunk_size];
#ifndef NDEBUG
  if(!ram) {
    cerr << "COuld not allocate ram!\n";
    exit(0);
  }
#endif
  
  Patch *patch = new Patch(signature, ram, entry.nvert, entry.nface);
  
  if(entry.patch_start != 0xffffffff) { //was allocated.
    assert(entry.disk_size != 0xffff);
    
    MFile::SetPosition((int64)entry.patch_start * (int64)chunk_size);
    
    if(signature.compr == 0) { //not compressed
      MFile::ReadBuffer(ram, entry.disk_size * chunk_size);
    } else {
      unsigned char *disk = new unsigned char[entry.disk_size * chunk_size];
      MFile::ReadBuffer(disk, entry.disk_size * chunk_size);
      
      patch->Decompress(entry.ram_size * chunk_size, 
			disk, entry.disk_size * chunk_size);
      delete []disk;
    } 
  } else {
    //zero all bytes... so compressio gets better with padding.
    memset(ram, 0, entry.ram_size * chunk_size);
  }
  ram_used += entry.ram_size;  
  entry.patch = patch;  
  return patch;
}

void Nexus::FlushPatch(unsigned int id) {
  Entry &entry = operator[](id);    
  assert(entry.patch);

  if(!MFile::IsReadOnly()) { //write back patch
    if(signature.compr) {
      unsigned int compressed_size;
      char *compressed = entry.patch->Compress(entry.ram_size * chunk_size,
					       compressed_size);
      if(entry.disk_size == 0xffff) {//allocate space 
	assert(entry.patch_start == 0xffffffff);
	entry.disk_size = (unsigned int)((compressed_size-1)/chunk_size) + 1;
	entry.patch_start = (unsigned int)(Length()/chunk_size);
	Redim(Length() + entry.disk_size * chunk_size);
      } else {
	//cerr << "OOOOPSPPPS not supported!" << endl;
	exit(-1);
      }
      MFile::SetPosition((int64)entry.patch_start * (int64)chunk_size);
      MFile::WriteBuffer(compressed, entry.disk_size * chunk_size);
      delete []compressed;
    } else {
      if(entry.disk_size == 0xffff) {
	entry.disk_size = entry.ram_size;
	entry.patch_start = (unsigned int)(Length()/chunk_size);
	Redim(Length() + entry.disk_size * chunk_size);
      }
      MFile::SetPosition((int64)entry.patch_start * (int64)chunk_size);
      MFile::WriteBuffer(entry.patch->fstart, entry.disk_size * chunk_size);
    }
  }

  delete [](entry.patch->fstart);
  delete entry.patch;  
  entry.patch = NULL;    
  ram_used -= entry.ram_size;      
}
