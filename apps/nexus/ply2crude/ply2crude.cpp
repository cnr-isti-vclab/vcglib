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
Revision 1.1  2004/06/23 00:08:05  ponchio
Created


****************************************************************************/

///TODO: allow other kinds of ply to be readed.

#include <string>
#include <iostream>
#include <wrap/ply/plylib.h>
#include "../crude.h"

using namespace std;
using namespace vcg;
using namespace vcg::ply;
using namespace nxs;

struct PlyVertex {
    float v[3];
};

struct PlyFace {
  unsigned int f[3];
  unsigned char flags;  
};

PropDescriptor plyprop1[3]= {
  {"vertex","x",T_FLOAT,T_FLOAT,offsetof(PlyVertex,v[0]),0,0,0,0,0},
  {"vertex","y",T_FLOAT,T_FLOAT,offsetof(PlyVertex,v[1]),0,0,0,0,0},
  {"vertex","z",T_FLOAT,T_FLOAT,offsetof(PlyVertex,v[2]),0,0,0,0,0}
};

PropDescriptor plyprop2[1]=	{
  {"face", "vertex_indices",T_INT,T_UINT,offsetof(PlyFace,f[0]),
   1,0,T_UCHAR,T_UCHAR,offsetof(PlyFace,flags) }
};

int main(int argc, char *argv[]) {
  if(argc <= 2) {
    cerr << "Usage: " << argv[0] << " <input1.ply> <...> <inputN.ply> <output>\n";
    return 0;
  }
  string output = argv[argc-1];
  //test last one is not a ply
  if(output.substr(output.size()-4, output.size()) == ".ply") {
    cerr << "Last argument is output (so not a .ply)\n";
    return -1;
  }
  Crude crude;
  if(!crude.Create(output, 0, 0)) {
    cerr << "Could not create crude output\n";
    return -1;
  }
  Box3f box;
  box.SetNull();
  for(int k = 1; k < argc-1; k++) {
    PlyFile pf;  
    //Opening ply file
    int val = pf.Open(argv[k], PlyFile::MODE_READ);
    if(val == -1) 	{	
      cerr << "Could not open file '" << argv[k] << "'\n";
      return false;
    }
    
    //testing for required vertex fields.
    if( pf.AddToRead(plyprop1[0])==-1 || 	
	pf.AddToRead(plyprop1[1])==-1 || 
	pf.AddToRead(plyprop1[2])==-1) 	{
      cerr << "Error Ply file has not one of the required elements :"
	   << "xyz coords\n";
      return false;
    }

    //testing for required face fields.
    if( pf.AddToRead(plyprop2[0])==-1 ) {
      cerr << "Error Ply file has not one of the required elements:"
	   << "faces\n";
      return false;
    }

    for(unsigned int i = 0; i < pf.elements.size(); i++) {
      if(!strcmp( pf.ElemName(i),"vertex")) {
	unsigned int n_vertices = pf.ElemNumber(i);
	unsigned int offset = crude.Vertices();
	crude.Resize(offset + n_vertices, crude.Faces());

	cerr << "Adding " << n_vertices << " n_vertices" << endl;
	pf.SetCurElement(i);   
	PlyVertex vertex;
	Point3f p;
	for(unsigned v = offset; v < offset + n_vertices; v++) {
	  pf.Read((void *) &vertex);
	  p[0] = vertex.v[0];
	  p[1] = vertex.v[1];
	  p[2] = vertex.v[2];
	  box.Add(p);
	  crude.SetVertex(v, vertex.v);
	}

      } else if( !strcmp( pf.ElemName(i),"face") ) {
	unsigned int n_faces = pf.ElemNumber(i);
	unsigned int offset = crude.Faces();
	crude.Resize(crude.Vertices(), offset + n_faces);

	cerr << "Adding " << n_faces << " n_faces" << endl;
	pf.SetCurElement(i);    
	PlyFace face;
	for(unsigned v = offset; v < offset + n_faces; v++) {
	  pf.Read((void *) &face);
	  assert(face.f[0] < crude.Vertices() && 
		 face.f[1] < crude.Vertices() && 
		 face.f[2] < crude.Vertices());

	  crude.SetFace(v, face.f);
	}
      }
    }
    pf.Destroy();
  }
  crude.GetBox() = box;
  
  return 0;
}
