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
Revision 1.1  2004/06/22 00:39:56  ponchio
Created


****************************************************************************/

#ifndef VFILE_H
#define VFILE_H

#include <errno.h>
//#include <hash_map>
#include <map>
#include <list>
#include <string>
#include <iostream>

/**Vector structure on file with simulated mmapping.
 * a priority queue of buffers is used 
 * TODO: port to over 4Gb usable space
 *       some mechanism to report errors?
 *       use an Iterator?
 */

template <class T> class VFile {
public:

  struct Buffer {
    unsigned int key;
    unsigned int size;
    T *data;
  };

  typedef std::list<Buffer>::iterator iterator;

 private:

  FILE *fp;  
  std::map<unsigned int, iterator> index;   //TODO move to hash_map 
  std::list<Buffer> buffers;
  unsigned int chunk_size;  //default buffer size
  unsigned int chunk_bits; //log2(chunk_size);
  unsigned int queue_size;
  unsigned int n_elements; //size of the vector

 public:
  
  VFile(): fp(NULL) {}
  ~VFile() { if(fp) close(); }                    
  bool create(const std::string &filename, 
	      unsigned int _chunk_bits = 12, 
	      unsigned int _queue_size = 1000) {

    fp = fopen(filename.c_str(), "wb+");
    if(!fp)
      return false;        
    chunk_bits = _chunk_bits;
    chunk_size = 1<<_chunk_bits;
    queue_size = _queue_size;
    n_elements = 0;    
    return true;
  }

  bool load(const std:: string &filename, 
	    unsigned int _chunk_bits = 12, 
	    unsigned int _queue_size = 1000) {

    fp = fopen(filename.c_str(), "rb+");
    if(!fp) return false;

    //troviamone la lunghezza
    fseek(fp, -1, SEEK_END);
     chunk_bits = _chunk_bits;
    chunk_size = 1<<_chunk_bits;
    queue_size = _queue_size;

    n_elements = ftell(fp)/ sizeof(T);   
    assert(n_elements >= chunk_size);    
    return true;
  }

  void close() {
    flush();
    fclose(fp);
    fp = 0;
  }

  void flush() {
    iterator i;
    for(i = buffers.begin(); i != buffers.end(); i++)       
      flushBuffer(*i);
    buffers.clear();
    index.clear();
  }

  void flushBuffer(Buffer buffer) {
    fseek(fp, buffer.key * sizeof(T), SEEK_SET);
    if(buffer.size != fwrite(buffer.data, sizeof(T), buffer.size, fp)) {
      cerr << "Could not write!" << endl;
      exit(0);
    }
    delete []buffer.data;
  }

  void resize(unsigned int elem) {
    if(elem < n_elements) 
      return;
    if(n_elements < chunk_size) 
      n_elements = chunk_size;
    while(elem > n_elements && n_elements < 256000000)
      n_elements *= 2;
    while(elem > n_elements)
      n_elements += 256000000;
    if(-1 == fseek(fp, n_elements * sizeof(T), SEEK_SET)) {
      cerr << "Could not resize!" << endl;
      exit(-1);
    }
    fwrite(&elem, sizeof(unsigned int), 1, fp);
  }

  /** Remember that T is a valid pointer only until next call of
   * getElement or setElement
   */
  T *getElement(unsigned int n) { 
    
    if(n > n_elements) {
      cerr << "Overflow!" << endl;
      return NULL;
    }    
    unsigned int chunk = (n >> chunk_bits) << chunk_bits;
    unsigned int offset = n - chunk;    
        
    if(index.count(chunk)) 
      return (*(index[chunk])).data + offset;

    if(buffers.size() > queue_size) {
      Buffer &buffer= buffers.back();
      flushBuffer(buffer);      
      index.erase(buffer.key);  
      buffers.pop_back();
    }
    
    Buffer buffer;
    buffer.key = chunk;
    buffer.data = new T[chunk_size * chunks];      
    buffer.size = chunks * chunk_size;
    if(fseek(fp, chunk * sizeof(T), SEEK_SET)) {
      cerr << "failed to fseek" << endl;
      return NULL;
    }
    
    if(buffer.size != fread(buffer.data, sizeof(T), buffer.size, fp)) {    
      if(!ferror(fp)) {
        cerr << "end of file" << endl;
      } else {     
	cerr << "failed reading!: " << errno << endl;
      }
      return NULL;
    }

    buffers.push_front(buffer);    
    index[chunk] = buffers.begin();      
    return buffer.data + offset;
  }

  /**use this for directly writing on the file... 
   * be careful to flush (unless you never readed or flushed)
   */

  void setElement(unsigned int i, T &t) {
    *getElement(i) = t;
  }

  void setPosition(unsigned int chunk) {
    fseek(fp, chunk * sizeof(T), SEEK_SET);
  }
  
  unsigned int size() { return n_elements; }
  unsigned int chunkSize() { return chunk_size; }
  unsigned int queueSize() { return queue_size; }
};

#endif
