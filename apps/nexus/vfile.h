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
Revision 1.11  2004/10/01 16:00:12  ponchio
Added include <assert.h>

Revision 1.10  2004/09/30 23:56:33  ponchio
Backup (added strips and normals)

Revision 1.9  2004/07/20 14:04:32  ponchio
Improved efficience in operator[]

Revision 1.8  2004/07/15 14:32:49  ponchio
Debug.

Revision 1.7  2004/07/05 15:49:39  ponchio
Windows (DevCpp, mingw) port.

Revision 1.6  2004/07/04 15:23:48  ponchio
Debug

Revision 1.5  2004/07/02 17:41:37  ponchio
Debug.

Revision 1.4  2004/07/02 13:02:39  ponchio
Added GetRegion, Read and Write

Revision 1.3  2004/07/01 21:36:54  ponchio
*** empty log message ***

Revision 1.2  2004/06/25 16:47:13  ponchio
Various debug

Revision 1.1  2004/06/24 14:32:45  ponchio
Moved from wrap/nexus

Revision 1.3  2004/06/22 15:32:09  ponchio
Tested

Revision 1.2  2004/06/22 10:27:16  ponchio
*** empty log message ***

Revision 1.1  2004/06/22 00:39:56  ponchio
Created


****************************************************************************/

#ifndef VFILE_H
#define VFILE_H

#include "file.h"

#include <assert.h>
#include <map>
#include <list>
#include <string>

/**Vector structure on file with simulated mmapping.
 * a priority queue of buffers is used 
 * TODO: port to over 4Gb usable space
 *       some mechanism to report errors?
 *       use an Iterator?
 */

namespace nxs {

template <class T> class VFile: public File {
 public:
  
  struct Buffer {
    unsigned int key;
    unsigned int size; //in number of elements
    T *data;
  };

 protected:
  unsigned int n_elements;
  std::list<Buffer> buffers;
  typedef typename std::list<Buffer>::iterator list_iterator;

  std::map<unsigned int, list_iterator> index;   //TODO move to hash_map 
  Buffer *last_buffer;

  unsigned int chunk_size; //default buffer size (expressed in number of T)
  unsigned int queue_size;

 public:
  class iterator {
  public:
    iterator(unsigned int p = 0, VFile *b = 0): n(p), buffer(b) {}
    T &operator*() { return (*buffer)[n]; }
    void operator++() { n++; }
    bool operator!=(const iterator &i) { return i.n != n; }
  private:
    unsigned int n;
    VFile *buffer;
  };
  
  VFile(): last_buffer(NULL) {}
  ~VFile() { Close(); }
  bool Create(const std::string &filename, 
	      unsigned int _chunk_size = 4096/sizeof(T), 
	      unsigned int _queue_size = 1000) {

    assert(_chunk_size > 0);
    n_elements = 0;
    last_buffer = NULL;
    chunk_size = _chunk_size;
    queue_size = _queue_size;

    return File::Create(filename);
  }

  bool Load(const std:: string &filename, 
	    unsigned int _chunk_size = 4096/sizeof(T), 
	    unsigned int _queue_size = 1000) {

    assert(_chunk_size > 0);
    last_buffer = NULL;
    chunk_size = _chunk_size;
    queue_size = _queue_size;

    if(!File::Load(filename)) return false;
    n_elements = size/sizeof(T);
    return true;
  }

  void Close() {
      Flush();
  }

  void Flush() {
    list_iterator i;
    for(i = buffers.begin(); i != buffers.end(); i++)       
      FlushBuffer(*i);
    buffers.clear();
    index.clear();
    last_buffer = NULL;
  }

  void FlushBuffer(Buffer buffer) {
    SetPosition(buffer.key * chunk_size * sizeof(T));
    WriteBuffer((char *)(buffer.data), buffer.size * sizeof(T));
    delete []buffer.data;
  }

  void Resize(unsigned int elem) {
    Flush();
    File::Resize(elem * sizeof(T));
    n_elements = elem;
  }

  /** Remember that T is a valid pointer only until next call of
   * getElement or setElement
   */
  T &operator[](unsigned int n) { 

    assert(n < n_elements);

    unsigned int chunk = n/chunk_size;
    unsigned int offset = n - chunk*chunk_size;
    assert(offset < chunk_size * sizeof(T));

    if(last_buffer && last_buffer->key == chunk) 
        return *(last_buffer->data + offset);
        
    if(index.count(chunk)) {
      last_buffer = &*index[chunk];
      return *((*(index[chunk])).data + offset);
    }
    
    if(buffers.size() > queue_size) {
      Buffer &buffer= buffers.back();
      assert(buffer.key != chunk);
      FlushBuffer(buffer);      
      index.erase(buffer.key);  
      buffers.pop_back();
    }
    
    Buffer buffer;
    buffer.key = chunk;
    buffer.size = chunk_size;

    if(buffer.size + chunk * chunk_size > n_elements)
      buffer.size = n_elements - chunk * chunk_size;

    buffer.data = new T[buffer.size];  

    buffers.push_front(buffer);   
    index[buffer.key] = buffers.begin();   
    last_buffer = &*buffers.begin();
    
    SetPosition(chunk * chunk_size * sizeof(T));
    ReadBuffer((char *)(buffer.data), buffer.size * sizeof(T));    

    return *(buffer.data + offset);
  }

  /** you can get a region instead of an element but:
      1)region must be Chunk aligned.
      2)you get impredictable results if regions overlap or mix with operator[]
  */
  T *GetRegion(unsigned int start, unsigned int size, bool flush = true) {
    assert(start + size <= n_elements);
    assert((size % chunk_size) == 0);
    assert((start % chunk_size) == 0);
    if(size == 0) return NULL;
    
    unsigned int chunk = start/chunk_size;

    if(index.count(chunk)) 
      return ((*(index[chunk])).data);
    
    while(flush && buffers.size() > queue_size) {
      Buffer &buffer= buffers.back();
      FlushBuffer(buffer);      
      index.erase(buffer.key);  
      buffers.pop_back();
    }
    
    Buffer buffer;
    buffer.key = chunk;
    buffer.size = size;
    buffer.data = new T[buffer.size];  

    buffers.push_front(buffer);    
    index[chunk] = buffers.begin();   

    SetPosition(chunk * chunk_size * sizeof(T));
    ReadBuffer((char *)(buffer.data), buffer.size * sizeof(T));
    return buffer.data;
  } 

  void PushBack(const T &t) {
    Resize(n_elements+1);
    operator[](n_elements-1) = t;
  }
  
  unsigned int Size() { return n_elements; }
  unsigned int ChunkSize() { return chunk_size; }
  unsigned int QueueSize() { return queue_size; }
  iterator Begin() { return iterator(0, this); }
  iterator End() { return iterator(Size(), this); }
};

}//namespace

#endif
