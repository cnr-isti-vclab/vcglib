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
Revision 1.1  2004/07/02 13:03:01  ponchio
Created


****************************************************************************/
#include <iostream>
#include <vector>
#include <set> //DEBUG
#include <map>

#include "../nexus.h"
#include "../crude.h"
#include "../vert_remap.h"
#include "../vfile.h"

using namespace std;
using namespace vcg;
using namespace nxs;


struct RemapLink {
  unsigned int rel_vert;
  unsigned int patch;
  unsigned int abs_vert;
};

int main(int argc, char *argv[]) {
  if(argc <= 2) {
    cerr << "Usage: " << argv[0] << " <crude input> <nexus output>\n";
    return -1;
  }
  Crude crude;
  if(!crude.Load(argv[1])) {
    cerr << "Could not load crude:" << argv[1] << endl;
    return -1;
  }

  VertRemap vert_remap;
  if(!vert_remap.Load(argv[1])) {
    cerr << "Could not load vert_remap files: " << argv[1] << "\n";
    return -1;
  }

  VFile<unsigned int> face_remap;
  if(!face_remap.Load(argv[1] + string(".rmf"))) {
    cerr << "Could not load face_remap file: " << argv[1] << ".frm\n";
    return -1;
  }



  Nexus nexus;
  if(!nexus.Create(argv[2])) {
    return -1;
  }
  


  //lets count faces, 
  vector<unsigned int> patch_faces;

  for(unsigned int i = 0; i < face_remap.Size(); i++) {
    unsigned int patch = face_remap[i];
    assert(patch != 0xffffffff);
    if(patch >= patch_faces.size())
      patch_faces.resize(patch+1, 0);
    patch_faces[patch]++;
    nexus.totface++;
  }

  cerr << "Vertremap size: " << vert_remap.Size() << endl;

  
  //lets count vertices
  vector<unsigned int> patch_verts;
  patch_verts.resize(patch_faces.size(), 0);

  unsigned int unreferenced = 0;
  // counting vertices using border_size for number of vertices  
  for(unsigned int i = 0; i < vert_remap.all.Size(); i++) {
    unsigned int patch = vert_remap.all[i];
    if(patch == 0xffffffff) {
      unreferenced++;
      continue;
    }
    assert(patch < patch_verts.size());
    patch_verts[patch]++;
    nexus.totvert++;
  }

  unsigned int totbord = 0;
  VFile<MFHash::Bucket> &border = vert_remap.borders.buffer;
  cerr << "border: " << border.Size() << endl;
  for(unsigned int i = 0; i < border.Size(); i++) {
    MFHash::Bucket &bucket = border[i];
    if(bucket.key == 0xffffffff) continue;
    unsigned int patch = bucket.value;
    assert(patch < patch_verts.size());
    patch_verts[patch]++;
    totbord++;
    //    nexus.totvert++; (?)
  }

  if(unreferenced)
    cerr << "Warning: found " << unreferenced << " vertices.\n";

  cerr << "Triangles found: " << nexus.totface << endl;
  cerr << "Vertex found: " << nexus.totvert << endl;
  cerr << "Borders found: " << totbord <<  endl;


  nexus.index.resize(patch_verts.size());
  //now that we know various sizes, lets allocate space
  for(unsigned int i = 0; i < nexus.index.size(); i++) {
    Nexus::Entry &entry = nexus.index[i];

    if(patch_faces[i] == 0 || patch_verts[i] == 0) 
      cerr << "Warning! Empty patch.\n";

    entry.patch_offset = nexus.totchunks;
    entry.patch_size = Patch::ChunkSize(patch_verts[i], patch_faces[i]);
    
    nexus.totchunks += entry.patch_size;
    entry.border_offset = 0;
  }

  nexus.patches.Resize(nexus.totchunks);

  //now we sort the faces into the patches (but still using absolute indexing
  //instead of relative indexing
  for(unsigned int i = 0; i < crude.face.Size(); i++) {
    Crude::Face &face = crude.face[i];
    unsigned int npatch = face_remap[i];

    Nexus::Entry &entry = nexus.index[npatch];
    Patch patch = nexus.GetPatch(npatch);
    if(entry.border_offset == 0) { //first time we find patch
      entry.border_offset = 0xffffffff;
      patch.VertSize() = patch_verts[npatch];
      patch.FaceSize() = 0;
    }

    Crude::Face *faces = (Crude::Face *)patch.VertBegin();
    faces[patch.FaceSize()] = face;
    patch.FaceSize()++;
  }


  //finally for every patch we collect the vertices
  //and fill the patch.
  for(unsigned int i = 0; i < nexus.index.size(); i++) {
    Patch patch = nexus.GetPatch(i);
    assert(patch.FaceSize() == patch_faces[i]);

    Nexus::Entry &entry = nexus.index[i];

    //make a copy of faces (we need to write there :P)
    Crude::Face *faces = new Crude::Face[patch_faces[i]];
    memcpy(faces, (Crude::Face *)patch.VertBegin(), 
	   patch.FaceSize() * sizeof(Crude::Face));

    //collect all vertices we need.
    //TODO an hash_map would be faster?
    map<unsigned int, unsigned short> remap;
    for(unsigned int k = 0; k < patch.FaceSize(); k++) {
      Crude::Face &face = faces[k];
      
      for(int j = 0; j < 3; j++) {
        if(!remap.count(face[j])) {          
	  assert(remap.size() < patch.VertSize());
	  Point3f &v = crude.vert[face[j]];
          patch.VertBegin()[remap.size()] = v;
	  entry.sphere.Add(v);
          remap[face[j]] = remap.size();
        }
	patch.FaceBegin()[k*3 + j] = remap[face[j]];
      }
    }

    /*
    

    //costruiamo la mappa dei bordi:
    border_remap.resize(border_offset + info.nbord);
    border_remap_offset.push_back(border_offset);
    border_remap_size.push_back(0);

    vector<unsigned int> bord;
    hash_map<unsigned int, unsigned short>::iterator m;         
    for(m = remap.begin(); m != remap.end(); m++) {
      unsigned int abs_vertex = (*m).first;
      unsigned short rel_vertex = (*m).second;
       vert_remap.getValues(abs_vertex, bord);
       if(bord.size() == 0) {
         cerr << "EEEEEEEEEKKKKKKK" << endl;
         exit(0);
       }
       if(bord.size() ==1) continue; // non e' un bordo.
       BigLink link;
       link.abs_vert = abs_vertex;
       link.rel_vert = rel_vertex;
       
       
       for(unsigned int j = 0; j < bord.size(); j++) {
         if(bord[j] == i) continue;
          link.patch = bord[j];                   
          border_remap.setElement(border_offset++, link);
          border_remap_size.back()++;
          
       }
       }*/    
    




    delete []faces;
  }
  for(unsigned int i = 0; i < nexus.index.size(); i++) {
    nexus.sphere.Add(nexus.index[i].sphere);
  }
  nexus.Close();
  return 0;
}
