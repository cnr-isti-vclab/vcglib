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

#include <unistd.h>
#include <errno.h>
//#include <hash_map>
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

template <class T> class VFile {
 public:
  
  struct Buffer {
    unsigned int key;
    unsigned int size; //in number of elements
    T *data;
  };

 private:

  FILE *fp;  
  std::list<Buffer> buffers;
  typedef typename std::list<Buffer>::iterator list_iterator;

  std::map<unsigned int, list_iterator> index;   //TODO move to hash_map 

  unsigned int chunk_size; //default buffer size (expressed in number of T)
  unsigned int queue_size;
  unsigned int n_elements; //size of the vector

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
  
  VFile(): fp(NULL) {}
  ~VFile() { Close(); }                    
  bool Create(const std::string &filename, 
	      unsigned int _chunk_size = 4096/sizeof(T), 
	      unsigned int _queue_size = 1000) {

    assert(_chunk_size > 0);
    chunk_size = _chunk_size;
    queue_size = _queue_size;
    n_elements = 0;    

    fp = fopen(filename.c_str(), "wb+");
    if(!fp) return false;        

    return true;
  }

  bool Load(const std:: string &filename, 
	    unsigned int _chunk_size = 4096/sizeof(T), 
	    unsigned int _queue_size = 1000) {

    assert(_chunk_size > 0);
    chunk_size = _chunk_size;
    queue_size = _queue_size;

    fp = fopen(filename.c_str(), "rb+");
    if(!fp) return false;

    fseek(fp, 0, SEEK_END);
    n_elements = ftell(fp)/ sizeof(T);   
    return true;
  }

  void Close() {
    if(fp) {
      Flush();
      fclose(fp);
      fp = NULL;
    }
  }

  void Flush() {
    list_iterator i;
    for(i = buffers.begin(); i != buffers.end(); i++)       
      FlushBuffer(*i);
    buffers.clear();
    index.clear();
  }

  void FlushBuffer(Buffer buffer) {
    fseek(fp, buffer.key * chunk_size * sizeof(T), SEEK_SET);
    if(buffer.size != fwrite(buffer.data, sizeof(T), buffer.size, fp)) {
      assert(0 && "Could not write");
    }
    delete []buffer.data;
  }

  void Resize(unsigned int elem) {
    assert(fp);
    Flush();
    if(elem > n_elements) {
      if(-1 == fseek(fp, elem*sizeof(T) -1, SEEK_SET)) {
	assert(0 && "Could not resize");
      }
      unsigned char a;
      fwrite(&a, sizeof(unsigned char), 1, fp);
    } else {
      //TODO optimize: we do not need flush for buffers over elem.
      int fd = fileno(fp);
      ftruncate(fd, elem*sizeof(T));
    }    
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

    if(index.count(chunk)) 
      return *((*(index[chunk])).data + offset);
    
    if(buffers.size() > queue_size) {
      Buffer &buffer= buffers.back();
      assert(buffer.key != chunk);
      FlushBuffer(buffer);      
      index.erase(buffer.key);  
      buffers.pop_back();
    }
    
    Buffer buffer;
    buffer.key = chunk;
    buffer.data = new T[chunk_size];  
    buffer.size = chunk_size;
    if(buffer.size + chunk * chunk_size > n_elements)
      buffer.size = n_elements -chunk * chunk_size;

    buffers.push_front(buffer);   
    index[buffer.key] = buffers.begin();   

    if(fseek(fp, chunk * chunk_size * sizeof(T), SEEK_SET)) {
      assert(0 && "failed to fseek");
      return *(buffer.data);
    }
    
    if(buffer.size != fread(buffer.data, sizeof(T), buffer.size, fp)) {    
      if(feof(fp)) {
	assert(0 && "end of file");
      } else {     
	assert(0 && "failed reading!");
      }
      return (*buffer.data);
    }

    return *(buffer.data + offset);
  }
  
  /** you can get a region instead of an element but:
      1)region must be Chunk aligned.
      2)you get impredictable results if regions overlap or mix with operator[]
  */
  T *GetRegion(unsigned int start, unsigned int size) {
    assert(start + size <= n_elements);

    unsigned int chunk = start/chunk_size;
    unsigned int offset = start - chunk*chunk_size;
    assert(offset == 0);

    if(index.count(chunk)) 
      return ((*(index[chunk])).data);
    
    if(buffers.size() > queue_size) {
      Buffer &buffer= buffers.back();
      FlushBuffer(buffer);      
      index.erase(buffer.key);  
      buffers.pop_back();
    }
    
    Buffer buffer;
    buffer.key = chunk;
    buffer.data = new T[chunk_size * size];  
    buffer.size = chunk_size * size;
    if(buffer.size + chunk * chunk_size > n_elements)
      buffer.size = -chunk * chunk_size + n_elements;

    buffers.push_front(buffer);    
    index[chunk] = buffers.begin();   

    if(fseek(fp, chunk * chunk_size * sizeof(T), SEEK_SET)) {
      assert(0 && "failed to fseek");
      return buffer.data;
    }
    
    if(buffer.size != fread(buffer.data, sizeof(T), buffer.size, fp)) {    
      if(feof(fp)) {
	assert(0 && "end of file");
      } else {     
	assert(0 && "failed reading!");
      }
      return buffer.data;
    }
    return buffer.data;
  }

  /** Allocate copy and return (remember to delete it) memory.
      size is number of elements. **/
  T *Read(unsigned int begin, unsigned int size) {
    //to make its simple we flush and read.
    Flush();
    T *buf = new T[size];
    if(fseek(fp, begin * sizeof(T), SEEK_SET)) {
      assert(0 && "failed to fseek");
      return buf;
    }
    if(size != fread(buf, sizeof(T), size, fp)) {
      if(feof(fp)) {
	assert(0 && "end of file");
      } else {     
	assert(0 && "failed reading!");
      }
    }
    return buf;
  }
  
  void Write(T *buf, unsigned int begin, unsigned int size) {
    //to make its simple we flush and read.
    Flush();
   if(fseek(fp, begin * sizeof(T), SEEK_SET)) {
      assert(0 && "failed to fseek");
   }
   if(size != fwrite(buf, sizeof(T), size, fp)) {
     if(feof(fp)) {
       assert(0 && "end of file");
     } else {     
       assert(0 && "failed reading!");
     }
   }
  }

  unsigned int Size() { return n_elements; }
  unsigned int ChunkSize() { return chunk_size; }
  unsigned int QueueSize() { return queue_size; }
  iterator Begin() { return iterator(0, this); }
  iterator End() { return iterator(Size(), this); }

 protected:
  void SetPosition(unsigned int chunk) {
    fseek(fp, chunk * chunk_size * sizeof(T), SEEK_SET);
  }
};

}//namespace

#endif
