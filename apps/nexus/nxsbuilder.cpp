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
Revision 1.3  2004/12/01 03:32:46  ponchio
Level 2 (debug).

Revision 1.2  2004/12/01 03:24:32  ponchio
Level 2.

Revision 1.1  2004/12/01 01:15:03  ponchio
Level 0.

Revision 1.26  2004/11/30 22:49:39  ponchio
Level 0.


****************************************************************************/

#ifdef WIN32
#include <wrap/system/getopt.h>
#else
#include <unistd.h>
#endif

#include <iostream>

#include "nxstypes.h"
#include "crude.h"
#include "remapping.h"
#include "decimate.h"
#include "watch.h"


using namespace std;
using namespace vcg;
using namespace nxs;

/*struct PIndex {
  int64 offset;
  unsigned int lenght;
};

class FaceIndex: public std::vector<PIndex> {};*/


void usage() { 
  cerr << "Usage: voronoinxs <crude> <output> [options]\n"
    "  Options:\n"
    " -f N: use N faces per patch (default 1000, max 32000)\n"
    " -t N: mini faces per patch (default 200)\n"
    " -l N: number of levels\n"
    " -s F: scaling factor (0 < F < 1, default 0.5)\n"
    " -o N: number of optimization steps\n"
    " -d <method>: decimation method: quadric, cluster. (default quadric)\n"
    " -b N: ram buffer size (in bytes)\n"
    " -p N: which fase perform:\n"
    "       0) Remap faces\n"
    "       1) Sort faces\n"
    "       2) Build patches\n"
    "       3) Build borders\n"
    "       4) Build levels\n\n"
    " If you use the step option, all other parameters MUST stay the same\n\n";

}


void FirstStep(const string &crudefile, const string &output,
	       unsigned int patch_size, unsigned int patch_threshold,
	       float scaling, unsigned int optimization_steps) {
  Crude crude;

  if(!crude.Load(crudefile, true)) {
    cerr << "Could not open crude input: " << crudefile << endl;
    exit(0);
  }
  
  if(patch_size > crude.vert.Size()/2) {
    cerr << "Patch size too big: " << patch_size << " * 2 > " 
	 << crude.vert.Size() << endl;
    exit(0);
  }
  
  if(patch_threshold == 0xffffffff)
    patch_threshold = patch_size/4; 
  

  VChain vchain;

  VFile<unsigned int> face_remap;
  if(!face_remap.Create(output + ".rfm")) {
    cerr << "Could not create remap files: " << output << ".rmf\n";
    exit(0);
  }
  face_remap.Resize(crude.Faces());

  VFile<Point3f> baricenters;
  if(!baricenters.Create(output + ".bvr")) {
    cerr << "Could not create temporary baricenters file\n";
    exit(0);
  } 
  baricenters.Resize(crude.Faces());
  for(unsigned int i = 0; i < crude.Faces(); i++) {
    baricenters[i] = crude.GetBari(i);
  }

  BlockIndex face_index;

  Remap(vchain, baricenters, face_remap, face_index, 
	patch_size, patch_threshold, 65000, scaling,
	optimization_steps);

  if(!vchain.Save(output + ".vchain")) {
    cerr << "Could not save file: " << output << ".vchain\n";
    exit(0);
  }
  if(!face_index.Save(output + ".rfi")) {
    cerr << "Could not save file: " << output << ".rmi\n";
    exit(0);
  }

  baricenters.Delete();
}


void SecondStep(const string &crudefile, const string &output) {
  Crude crude;
  
  if(!crude.Load(crudefile, true)) {
    cerr << "Could not open crude input: " << crudefile << endl;
    exit(0);
  }
  VFile<unsigned int> face_remap;
  if(!face_remap.Load(output + ".rfm", true)) {
    cerr << "Could not load: " << output << ".rmf\n;";
    exit(0);
  }
  assert(face_remap.Size() == crude.Faces());

  VFile<Crude::Face> sorted;
  if(!sorted.Create(output + ".faces")) {
    cerr << "Could not create sorted faces\n";
    exit(0);
  }
  sorted.Resize(face_remap.Size());

  BlockIndex face_index;
  if(!face_index.Load(output + ".rfi")) {
    cerr << "Could not load index\n";
    exit(0);
  }
  //  cerr << "Face index size: " << face_index.size() << endl;

  //Sorting now.
  vector<unsigned int> done;
  done.resize(face_index.size(), 0);

  for(unsigned int i = 0; i < face_remap.Size(); i++) {
    unsigned int patch = face_remap[i];
    assert(patch < face_index.size());
    assert(patch >= 0);
    int64 offset = face_index[patch].offset + done[patch]++;
    sorted[offset] = crude.GetFace(i);
  }

#ifndef NDEBUG
  for(int i = 0; i < done.size(); i++)
    assert(done[i] == face_index[i].size);
#endif

  /*#ifndef NDEBUG
  //test:
  for(unsigned int i = 0; i < sorted.Size(); i++) {
  Crude::Face face = sorted[i];
  assert(face[0] < crude.Vertices());
  assert(face[1] < crude.Vertices());
  assert(face[2] < crude.Vertices());
  }
  #endif*/

  //once sorted
  crude.Close();
  sorted.Close();

  /*  TODO fix this (after debug!)
      WARNING what if multiple files?

  if(0 != unlink((crudefile + ".crf").c_str())) {
    cerr << "Could not remove " << crudefile << ".crf\n";
    exit(0);
  }
  if(0 != rename((output + ".faces").c_str(), (crudefile + ".crf").c_str())) {
    cerr << "Could not rename to: " << crudefile + ".crf\n";
    exit(0);
  }
  face_remap.Close();  */
  //TODO remove the file... (after finishing debug!)
  //  face_remap.Delete();
}

void ThirdStep(const string &crudefile, const string &output,
	       unsigned int chunk_size) {

  cerr << "Third step!\n";
  Crude crude;
  
  if(!crude.Load(crudefile, true)) {
    cerr << "Could not open crude input: " << crudefile << endl;
    exit(0);
  }
  
  VFile<Crude::Face> sorted;
  if(!sorted.Load(crudefile + ".faces", true)) {
    cerr << "Could not load sorted faces\n";
    exit(0);
  }
  BlockIndex face_index;
  if(!face_index.Load(output + ".rfi")) {
    cerr << "Could not load index\n";
    exit(0);
  }
  
  VFile<unsigned int> vert_remap;
  if(!vert_remap.Create(output + ".rvm")) {
    cerr << "Could not create: " << output << ".rvm\n";
    exit(0);
  }
  
  BlockIndex vert_index;

  Nexus nexus;
  //TODO here i really need no ram_buffer.....
  nexus.patches.SetRamBufferSize(0);
  if(!nexus.Create(output, NXS_FACES, chunk_size)) {
    cerr << "Could not create nexus output: " << output << endl;
    getchar();
    exit(0);
  }

  Report report(face_index.size());
  for(unsigned int patch = 0; patch < face_index.size(); patch++) {
    report.Step(patch);

    unsigned int vcount = 0;
    unsigned int fcount = 0;
    map<unsigned int, unsigned short> vremap;

    vector<Point3f> vertices;
    vector<unsigned short> faces;

    int64 &offset = face_index[patch].offset;
    unsigned int size = face_index[patch].size;
    for(unsigned int i = offset; i < offset + size; i++) {
      //TODO fix this after debug
      //      Crude::Face face = crude.GetFace(i);
      Crude::Face face = sorted[i];
      if(face[0] == face[1] || face[1] == face[2] || face[0] == face[2]) 
	continue; //degenerate
      for(int j = 0; j < 3; j++) {
	assert(face[j] < crude.Vertices());
	if(!vremap.count(face[j])) {          
	  Point3f &v = crude.vert[face[j]];
	  vertices.push_back(v);
	  vremap[face[j]] = vcount++;
	}
	faces.push_back(vremap[face[j]]);
	fcount++;
      }
    }
    assert(vcount == vertices.size());
    assert(fcount == faces.size());

    if(vcount > 65000) {
      cerr << "Too many vertices in this patch: " 
	   << vcount << endl;
      exit(0);
    }

    if(fcount > 65000) {
      cerr << "Too many faces in this patch: " 
	   << fcount << endl;
      exit(0);
    }
    unsigned int patch_idx = nexus.AddPatch(vertices.size(),
					    faces.size()/3,
					    0); //no borders!
    Patch &patch = nexus.GetPatch(patch_idx);
    memcpy(patch.FaceBegin(), &*faces.begin(), fcount * sizeof(short));
    memcpy(patch.VertBegin(), &*vertices.begin(), vcount * sizeof(Point3f));
    for(int i = 0; i < vertices.size(); i++)
      nexus.index[patch_idx].sphere.Add(vertices[i]);

    //saving vert_remap
    int64 vroffset = vert_remap.Size();
    vert_index.push_back(BlockEntry(vroffset, vcount));
    vert_remap.Resize(vroffset + vcount);

    map<unsigned int, unsigned short>::iterator r;
    for(r = vremap.begin(); r != vremap.end(); r++) {
      assert((*r).second < vcount);
      assert(vroffset + (*r).second < vert_remap.Size());
      vert_remap[vroffset + (*r).second] = (*r).first;
    }
  }

  //we can now update bounding sphere.
  for(unsigned int i = 0; i < nexus.index.size(); i++) 
    nexus.sphere.Add(nexus.index[i].sphere);
  
  Nexus::Update update;
  for(unsigned int i = 1; i < nexus.index.size(); i++) {
    update.created.push_back(i);
  }
  nexus.history.push_back(update);
  
  update.created.clear();
  update.created.push_back(0);
  for(unsigned int i = 1; i < nexus.index.size(); i++) {
    update.erased.push_back(i);
  }
  nexus.history.push_back(update);

  if(!vert_index.Save(output + ".rvi")) {
    cerr << "Could not save: " << output << ".rvi\n";
    exit(0);
  }
}

void FourthStep(const string &crudefile, const string &output, 
		unsigned int ram_buffer) {
  Nexus nexus;

  nexus.patches.SetRamBufferSize(ram_buffer);
  if(!nexus.Load(output)) {
    cerr << "Could not load nexus " << output << endl;
    exit(0);
  }
  //TODO Clear borders in case of failure!

  VFile<unsigned int> vert_remap;
  if(!vert_remap.Load(crudefile + ".rvm", true)) {
    cerr << "Could not load: " << crudefile << ".rvm\n";
    exit(0);
  }

  BlockIndex vert_index;
  if(!vert_index.Load(output + ".rvi")) {
    cerr << "Could not load index\n";
    exit(0);
  }
  Report report(nexus.index.size());

  for(int start = 0; start < nexus.index.size(); start++) {
    report.Step(start);
    Nexus::PatchInfo &s_entry = nexus.index[start];

    vector<Link> links;
    
    map<unsigned int, unsigned short> vremap;
    for(unsigned int i = 0; i < vert_index[start].size; i++) {
      unsigned int global = vert_remap[vert_index[start].offset + i];
      vremap[global] = i;
    }

    for(int end = 0; end < nexus.index.size(); end++) {
      if(start == end) continue;

      Nexus::PatchInfo &e_entry = nexus.index[end];
      float dist = Distance(s_entry.sphere, e_entry.sphere);
      if(dist > 0.1) continue;
      
      for(unsigned int i = 0; i < vert_index[end].size; i++) {
	unsigned int global = vert_remap[vert_index[end].offset + i];
	if(vremap.count(global)) {
	  Link link;
	  link.start_vert = vremap[global];
	  link.end_vert = i;
	  link.end_patch = end;
	  links.push_back(link);
	}
      }
    }
    //TODO Horribili visu (interfaccia di cacca!)
    nexus.borders.ResizeBorder(start, 3 * links.size());
    nexus.borders.borders[start].border_used = links.size();
    Border border = nexus.GetBorder(start);
    memcpy(&(border[0]), &*links.begin(), links.size() * sizeof(Link));
  }
}

void FifthStep() {
}

int main(int argc, char *argv[]) {
  
  /* Parameters: */
  unsigned int patch_size = 1000;  //step 0
  unsigned int patch_threshold = 0xffffffff; //Step 0
  float scaling = 0.5; //step 0, 4
  unsigned int optimization_steps = 5; //step 0, 4

  Decimation decimation = CLUSTER; //step 4
  unsigned int max_level = 0xffffffff; //step 4

  unsigned int ram_buffer = 128000000; //step 2, 3, 4
  unsigned int chunk_size = 1024;      //step 2, 3, 4
  int step = -1; //means all of them.
  
  int option;
  while((option = getopt(argc, argv, "f:t:l:s:d:o:b:c:p:")) != EOF) {
    switch(option) {
    case 'f': patch_size = atoi(optarg);
      if(patch_size == 0 || patch_size > 32000) {
	cerr << "Invalid number of faces per patch: " << optarg << endl;
	return -1;
      }
      break;
    case 't': patch_threshold = atoi(optarg);
      if(patch_threshold == 0 || patch_threshold > patch_size) {
	cerr << "Invalid patch threshold: " << optarg << endl;
	return -1;
      }
      break;
    case 'l': max_level = atoi(optarg);
      if(max_level == 0) {
	cerr << "Invalid number of levels: " << optarg << endl;
	return -1;
      }
      break;
    case 's': scaling = (float)atof(optarg);
      if(scaling <= 0 || scaling >= 1) {
	cerr << "Invalid scaling: " << optarg << endl;
	cerr << "Must be 0 < scaling < 1" << endl;
      }
      break;
    case 'd': 
      if(!strcmp("quadric", optarg)) 
	decimation = QUADRIC;
      else if(!strcmp("cluster", optarg)) 
	decimation = CLUSTER;
      else {
	cerr << "Unknown decimation method: " << optarg << endl;
	return -1;
      }
      break;
    case 'o': optimization_steps = atoi(optarg); break;
    case 'p': step = atoi(optarg); break;
    case 'b': ram_buffer = atoi(optarg); 
      if(ram_buffer == 0) {
	cerr << "Invalid ram buffer: " << optarg << endl;
	return -1;
      }
      break;
    case 'c': chunk_size = atoi(optarg);
      if(chunk_size == 0) {
	cerr << "Invalid chunk size: " << optarg << endl;
	return -1;
      }
      break;
    default: cerr << "Unknown option: " << (char)option << endl;
      return -1;
    }
  }

  if(optind != argc -2) {
    usage();
    return -1;
  }
  string crudefile = argv[optind];
  string output = argv[optind+1];

  if(step < 0 || step == 0)
    FirstStep(crudefile, output, patch_size, patch_threshold,
	      scaling, optimization_steps);
  if(step < 0 || step == 1)
    SecondStep(crudefile, output);

  if(step < 0 || step == 2)
    ThirdStep(crudefile, output, chunk_size);

  if(step < 0 || step == 3)
    FourthStep(crudefile, output, ram_buffer);
  return 0;
}

