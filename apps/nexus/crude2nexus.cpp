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
Revision 1.3  2004/07/15 14:32:49  ponchio
Debug.

Revision 1.2  2004/07/05 15:49:39  ponchio
Windows (DevCpp, mingw) port.

Revision 1.1  2004/07/04 15:30:00  ponchio
Changed directory structure.

Revision 1.3  2004/07/04 14:28:05  ponchio
Finito e debuggato

Revision 1.2  2004/07/02 17:42:12  ponchio
Backup.

Revision 1.1  2004/07/02 13:03:01  ponchio
Created


****************************************************************************/
#include <iostream>
#include <vector>
#include <set> //DEBUG
#include <map>

#include "nexus.h"
#include "crude.h"
#include "vert_remap.h"
#include "vfile.h"

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

  
  VFile<RemapLink> border_remap;
  if(!border_remap.Create(argv[1] + string(".tmp"))) {
    cerr << "Could not create temporary border remap file\n";
    return -1;
  }

  Nexus nexus;
  if(!nexus.Create(argv[2])) {
    return -1;
  }
  


  //lets count faces, 
  vector<unsigned int> patch_faces;
  unsigned int totface = 0;
  unsigned int totvert = 0;

  for(unsigned int i = 0; i < face_remap.Size(); i++) {
    unsigned int patch = face_remap[i];
    assert(patch != 0xffffffff);
    if(patch >= patch_faces.size())
      patch_faces.resize(patch+1, 0);
    patch_faces[patch]++;
    totface++;
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
    totvert++;
  }

  unsigned int totbord = 0;
  VFile<MFHash::Bucket> &border = vert_remap.borders.buffer;
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
    cerr << "Warning: found " << unreferenced << " unreferenced vertices.\n";

  cerr << "Triangles found: " << totface << endl;
  cerr << "Vertex found: " << totvert << endl;
  cerr << "Borders found: " << totbord <<  endl;


  nexus.index.resize(patch_verts.size());

  unsigned int totchunks = 0;
  //now that we know various sizes, lets allocate space
  for(unsigned int i = 0; i < nexus.index.size(); i++) {
    Nexus::Entry &entry = nexus.index[i];

    if(patch_faces[i] == 0 || patch_verts[i] == 0) 
      cerr << "Warning! Empty patch.\n";

    entry.patch_start = totchunks;
    entry.patch_size = Patch::ChunkSize(patch_verts[i], patch_faces[i]);
    
    totchunks += entry.patch_size;
    entry.border_start = 0;
  }

  nexus.patches.Resize(totchunks);

  //now we sort the faces into the patches (but still using absolute indexing
  //instead of relative indexing
  for(unsigned int i = 0; i < crude.face.Size(); i++) {
    Crude::Face &face = crude.face[i];
    unsigned int npatch = face_remap[i];

    Nexus::Entry &entry = nexus.index[npatch];
    Patch patch = nexus.GetPatch(npatch);
    if(entry.border_start == 0) { //first time we find patch
      entry.border_start = 0xffffffff;
      entry.nvert = patch_verts[npatch];
      //      patch.VertSize() = patch_verts[npatch];
      entry.nface = 0;
      //      patch.FaceSize() = 0;
    }

    Crude::Face *faces = (Crude::Face *)patch.VertBegin();
    faces[entry.nface] = face;
    entry.nface++;
    //    faces[patch.FaceSize()] = face;
    //    patch.FaceSize()++;
  }


  //finally for every patch we collect the vertices
  //and fill the patch.
  //we need to remember start and size in border_remap;
  //  vector<unsigned int> border_start;
  //  vector<unsigned int> border_size;

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
    unsigned int count = 0;
    map<unsigned int, unsigned short> remap;
    for(unsigned int k = 0; k < patch.FaceSize(); k++) {
      Crude::Face &face = faces[k];
      
      for(int j = 0; j < 3; j++) {
        if(!remap.count(face[j])) {          
	  assert(count < patch.VertSize());
	  Point3f &v = crude.vert[face[j]];
          patch.VertBegin()[remap.size()] = v;
	  entry.sphere.Add(v);
          remap[face[j]] = count++;
        }
	patch.FaceBegin()[k*3 + j] = remap[face[j]];
      }
    }
    assert(count == remap.size());
    assert(entry.nvert == remap.size());

    //record start of border:
    entry.border_start = border_remap.Size();

    //TODO hash_set?
    set<unsigned int> border_patches;
    map<unsigned int, unsigned short>::iterator m;
    for(m = remap.begin(); m != remap.end(); m++) {
      RemapLink link;
      link.abs_vert = (*m).first;
      link.rel_vert = (*m).second;

      vert_remap.GetValues(link.abs_vert, border_patches);
      assert(border_patches.size() >= 1);
      if(border_patches.size() == 1) continue; //its not a border

      set<unsigned int>::iterator s;
      for(s = border_patches.begin(); s != border_patches.end(); s++) {
	if((*s) == i) continue; 
	link.patch = *s;
	border_remap.PushBack(link);
      }
    }
    //and number of borders:
    entry.border_size = border_remap.Size() - entry.border_start;
    delete []faces;
  }
  //we can now update bounding sphere.
  for(unsigned int i = 0; i < nexus.index.size(); i++) 
    nexus.sphere.Add(nexus.index[i].sphere);

  //and last convert RemapLinks into Links
  nexus.borders.Resize(border_remap.Size());

  for(unsigned int i = 0; i < nexus.index.size(); i++) {
    Nexus::Entry &local = nexus.index[i];

    // K is the main iterator (where we write to in nexus.borders)
    for(unsigned int k = local.border_start;
	k < local.border_start + local.border_size; k++) {
      
      RemapLink start_link = border_remap[k];
      assert(start_link.rel_vert < local.nvert);

      Nexus::Entry &remote = nexus.index[start_link.patch];

      bool found = false;
      for(unsigned int j = remote.border_start;
	  j < remote.border_start + remote.border_size; j++) {
	
	RemapLink end_link = border_remap[j];
	assert(end_link.rel_vert < remote.nvert);

	if(start_link.abs_vert == end_link.abs_vert &&
	   end_link.patch == i) { //found the match
	  assert(!found);
	  nexus.borders[k] = Link(start_link.rel_vert, 
				  end_link.rel_vert, start_link.patch);
	  found = true;
	}
      }
      assert(nexus.borders[k].start_vert < local.nvert);
      assert(found);
    }
  }
  nexus.borders.Flush();

  //Checking border consistency:
  for(unsigned int i = 0; i < nexus.index.size(); i++) {
    Border border = nexus.GetBorder(i);
    Nexus::Entry &entry = nexus.index[i];
    for(unsigned int k = 0; k < border.Size(); k++) {
      Link &link = border[k];
      if(link.start_vert >= entry.nvert) {
	cerr << "K: " << k << endl;
	cerr << "patch: " << i << " nvert: " << entry.nvert << " startv: " 
	     << link.start_vert << endl;
	cerr << "bstart: " << entry.border_start 
	     << "bsize: " << entry.border_size << endl;
      }
      assert(link.end_patch < nexus.index.size());
      assert(link.start_vert < entry.nvert);
      Nexus::Entry &remote = nexus.index[link.end_patch];
      assert(link.end_vert < remote.nvert);
    }
    
  }

  nexus.Close();
  return 0;
}
