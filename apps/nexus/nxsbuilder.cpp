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
Revision 1.14  2005/01/21 17:09:13  ponchio
Porting and debug.

Revision 1.13  2005/01/17 17:35:48  ponchio
Small changes and adding realtime extr.

Revision 1.12  2005/01/14 15:25:29  ponchio
Revolution.

Revision 1.11  2004/12/13 00:44:48  ponchio
Lotsa changes...

Revision 1.10  2004/12/09 22:33:28  ponchio
Different splitting optimization.

Revision 1.9  2004/12/04 13:22:55  ponchio
*** empty log message ***

Revision 1.8  2004/12/03 21:19:00  ponchio
Fixed a couple of memory leak...

Revision 1.7  2004/12/03 01:20:56  ponchio
Debug

Revision 1.6  2004/12/02 20:22:42  ponchio
Level 5;

Revision 1.5  2004/12/02 19:10:18  ponchio
Bounding sphere fix.

Revision 1.4  2004/12/01 16:00:35  ponchio
Level 3

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

#ifdef WIN32
#include <hash_map>
#endif

#include <iostream>


#include "nxstypes.h"
#include "crude.h"
#include "remapping.h"
#include "decimate.h"
#include "fragment.h"
#include "nxsalgo.h"
#include "nxsdispatcher.h"
#include "watch.h"


using namespace std;
using namespace vcg;
using namespace nxs;

void BuildFragment(Nexus &nexus, VPartition &part,
		   set<unsigned int> &patches, 
		   Fragment &fragment);

void SaveFragment(Nexus &nexus, VChain &chain,
		  Fragment &fragin,
		  Fragment &fragout);

void ReverseHistory(vector<History::Update> &history);


unsigned int current_level;
vector<unsigned int> patch_levels;

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
    int64 offset = face_index[patch].offset + done[patch]++;
    sorted[offset] = crude.GetFace(i);
  }

#ifndef NDEBUG
  for(int i = 0; i < done.size(); i++)
    assert(done[i] == face_index[i].size);
#endif

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
  nexus.MaxRam() = 0;
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

    //TODO deal with this case adding another patch at the end... its not that difficult!
    
    //This can happen on degenerate cases when we have a lot of detached faces....
    if(vcount > 65000 && fcount > 65000) {
      cerr << "Too many vertices or faces in patch: " << patch << " v: " << vcount 
        << " f: " << fcount << endl;
      exit(0);
    }
  
    unsigned int patch_idx = nexus.AddPatch(vertices.size(),
					    faces.size()/3,
					    0); //no borders!
    Patch &patch = nexus.GetPatch(patch_idx);
    memcpy(patch.FaceBegin(), &*faces.begin(), fcount * sizeof(short));
    memcpy(patch.VertBegin(), &*vertices.begin(), vcount * sizeof(Point3f));

    Sphere3f &sphere = nexus[patch_idx].sphere;
    for(int i = 0; i < vertices.size(); i++)
      sphere.Add(vertices[i]);
    sphere.Radius() *= 1.01;

#ifndef NDEBUG
    for(int i = 0; i < vertices.size(); i++) {
      assert(sphere.IsIn(vertices[i]));
    }
#endif
   

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
    if(vcount < 100) {
      cerr << "Small patch: " << vcount << "\n";
    }
  }  

  //we can now update bounding sphere.
  for(unsigned int i = 0; i < nexus.size(); i++) 
    nexus.sphere.Add(nexus[i].sphere);
  
  History::Update update;
  for(unsigned int i = 1; i < nexus.size(); i++) {
    update.created.push_back(i);
  }
  nexus.history.updates.push_back(update);
  
  update.created.clear();
  update.created.push_back(0);
  for(unsigned int i = 1; i < nexus.size(); i++) {
    update.erased.push_back(i);
  }
  nexus.history.updates.push_back(update);

  if(!vert_index.Save(output + ".rvi")) {
    cerr << "Could not save: " << output << ".rvi\n";
    exit(0);
  }
  nexus.Close();
}

void FourthStep(const string &crudefile, const string &output, 
		unsigned int ram_buffer) {
  Nexus nexus;  
  if(!nexus.Load(output)) {
    cerr << "Could not load nexus " << output << endl;
    exit(0);
  }
  nexus.MaxRam() = ram_buffer / nexus.chunk_size;  
  //TODO Clear borders in case of failure!

  VFile<unsigned int> vert_remap;
  if(!vert_remap.Load(output + ".rvm", true)) {
    cerr << "Could not load: " << crudefile << ".rvm\n";
    exit(0);
  }  

  BlockIndex vert_index;
  if(!vert_index.Load(output + ".rvi")) {
    cerr << "Could not load index\n";
    exit(0);
  }

  /*float max_radius = 0;
  VPartition grid;
  for(int start = 0; start < nexus.size(); start++) {    
    Entry &entry = nexus[start];
    Sphere3f &sphere = entry.sphere;
    if(sphere.Radius() > max_radius) max_radius = sphere.Radius();
    grid.push_back(sphere.Center());
  }
  grid.Init();
  max_radius *= max_radius * 4;
  vector<int> nears;
  vector<float> dists;   */
    

  Report report(nexus.size());

  for(int start = 0; start < nexus.size(); start++) {    
    report.Step(start);
    Entry &s_entry = nexus[start];
    Sphere3f &sphere = s_entry.sphere;

    vector<Link> links;   
#ifdef WIN32
    hash_map<unsigned int, unsigned short> vremap;
#else
    map<unsigned int, unsigned short> vremap;
#endif
    for(unsigned int i = 0; i < vert_index[start].size; i++) {
      unsigned int global = vert_remap[vert_index[start].offset + i];      
      vremap[global] = i;
    }    

    /*unsigned int n_nears = 10;    
    while(1) {
      if(n_nears > grid.size()) n_nears = grid.size();  
      grid.Closest(sphere.Center(), n_nears, nears, dists);
      if(dists.back() > max_radius) break;
      if(n_nears == grid.size()) break;
      n_nears *= 2;      
    }
      
    for(int n = 0; n < nears.size(); n++) {
      unsigned int end = nears[n]; */
    for(int end = 0; end < nexus.size(); end++) {
      if(start == end) continue;      

      Entry &e_entry = nexus[end];
      float dist = Distance(s_entry.sphere, e_entry.sphere);
      
      if(dist > s_entry.sphere.Radius() + e_entry.sphere.Radius()) {
        continue;
      }
      
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

    Border &border = nexus.GetBorder(start);
    nexus.borders.ResizeBorder(start, 3 * links.size());
    border.used = links.size();        
    memcpy(&(border[0]), &*links.begin(), links.size() * sizeof(Link));
  }
}

void FifthStep(const string &crudefile, const string &output, 
	       unsigned int ram_buffer, 
	       unsigned int optimization_steps, 
	       unsigned int patch_size, 
	       unsigned int patch_threshold,
	       Decimation decimation, 
	       float scaling, 
	       unsigned int max_level) {

  Nexus nexus;  
  if(!nexus.Load(output)) {
    cerr << "Could not load nexus " << output << endl;
    exit(0);
  }
  nexus.MaxRam() = ram_buffer / nexus.chunk_size;

  VChain vchain;
  if(!vchain.Load(output + ".vchain")) {
    cerr << "Could not load : " << output << ".vchain\n";
    exit(0);
  }
  nexus.history.Clear();
  History::Update update;
  for(unsigned int i = 0; i < nexus.size(); i++) {
    update.created.push_back(i);
    patch_levels.push_back(0);
  }
  nexus.history.updates.push_back(update); 
  Unify(nexus, 0.0f);
  //  nexus.Unify();
  nexus.Flush();


  Dispatcher dispatcher(&nexus, &vchain);
  dispatcher.mode = decimation;
  dispatcher.scaling = scaling;
  if(!dispatcher.Init("servers.txt")) {
    cerr << "Could not parse server file: " << "servers.txt"
	 << " proceding locally\n";
  }
  
  unsigned int oldoffset = 0;
  for(unsigned int level = 1; level < max_level; level++) {
    current_level = level;
    cerr << "Level: " << level << endl;

    unsigned int newoffset = nexus.size();
    BuildLevel(vchain, nexus, oldoffset, scaling, 
	       patch_size, patch_threshold, 65000,
	       optimization_steps);
    
    Report report(vchain.oldfragments.size());

    unsigned int fcount = 0;
    map<unsigned int, set<unsigned int> >::iterator fragment;
    for(fragment = vchain.oldfragments.begin(); 
	fragment != vchain.oldfragments.end(); fragment++) {
      report.Step(fcount++);

      Fragment *fragin = new Fragment;
      BuildFragment(nexus, *vchain[level+1], 
		    (*fragment).second, *fragin);

      dispatcher.SendFragment(fragin);


      /*
      //this can be executed on a remote host

      //TODO move this part to remote....
      vector<Point3f> newvert;
      vector<unsigned int> newface;
      vector<BigLink> newbord;
      Join(fragin, newvert, newface, newbord);

      float error = Decimate(decimation,
			     (unsigned int)((newface.size()/3) * scaling), 
			     newvert, newface, newbord);
      Fragment fragout;
      fragout.error = error;
      fragout.id = fragin.id;
      fragout.seeds = fragin.seeds;
      fragout.seeds_id = fragin.seeds_id;
      Split(fragout, newvert, newface, newbord);//, vchain.levels[level+1]); 

      
      SaveFragment(nexus, vchain, fragin, fragout);
      */
      dispatcher.processmsgs();
    }
    //TODO porcata!!!!
    while(dispatcher.frags.size()) {
      //      cerr << "frags: " << dispatcher.frags.size() << endl;
      dispatcher.processmsgs();
    }

    report.Finish();

    if(vchain.oldfragments.size() == 1) break;

    vchain.oldfragments = vchain.newfragments;
    oldoffset = newoffset;
  }

  //last level clean history:
  update.created.clear();
  update.erased.clear();
  map<unsigned int, set<unsigned int> >::iterator fragment;
  for(fragment = vchain.newfragments.begin(); 
      fragment != vchain.newfragments.end(); fragment++) {
    set<unsigned int> &fcells = (*fragment).second;
    set<unsigned int>::iterator s;
    for(s = fcells.begin(); s != fcells.end(); s++)
      update.erased.push_back(*s);
  }
  nexus.history.updates.push_back(update);
  ReverseHistory(nexus.history.updates);

  //  TestBorders(nexus);
  nexus.Close();

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

  if(step < 0 || step == 4) 
    FifthStep(crudefile, output, 
	      ram_buffer, 
	      optimization_steps, 
	      patch_size, patch_threshold,
	      decimation, 
	      scaling, max_level);
  return 0;
}

void BuildFragment(Nexus &nexus, VPartition &part,
		   set<unsigned int> &patches, 
		   Fragment &fragment) {

  set<unsigned int>::iterator f;
  for(f = patches.begin(); f != patches.end(); f++) {
    fragment.pieces.push_back(NxsPatch());
    NxsPatch &nxs = fragment.pieces.back();
    nxs.patch = *f;

    Patch &patch = nexus.GetPatch(*f);
    Border &border = nexus.GetBorder(*f);

    for(unsigned int k = 0; k < patch.nf; k++) {
      assert(patch.Face(k)[0] != patch.Face(k)[1]);
      assert(patch.Face(k)[0] != patch.Face(k)[2]);
      assert(patch.Face(k)[1] != patch.Face(k)[2]);
    }


    nxs.vert.resize(patch.nv);
    nxs.face.resize(patch.nf * 3);
    memcpy(&*nxs.vert.begin(), patch.VertBegin(), patch.nv * sizeof(Point3f));
    memcpy(&*nxs.face.begin(), patch.FaceBegin(), patch.nf * 3*sizeof(short));
    for(unsigned int i = 0; i < border.Size(); i++) {
      Link &link = border[i];
      if(!link.IsNull() && 
	       patch_levels[link.end_patch] == current_level-1) {
	      assert(link.end_patch != *f);
	      nxs.bord.push_back(link);
      }
    }
  }

  set<unsigned int> seeds;
  vector<int> nears;
  vector<float> dists;
  int nnears = 10;
  if(part.size() < 10) nnears = part.size();
  for(f = patches.begin(); f != patches.end(); f++) {
    Point3f &center = nexus[*f].sphere.Center();
    part.Closest(center, nnears, nears, dists);
    for(int i = 0; i < nnears; i++) 
      seeds.insert(nears[i]);
  }
  for(f = seeds.begin(); f != seeds.end(); f++) {
    Point3f &p = part[*f];
    fragment.seeds.push_back(p);
    fragment.seeds_id.push_back(*f);
  }
}

void SaveFragment(Nexus &nexus, VChain &chain,
		  Fragment &fragin,
		  Fragment &fragout) {

  set<unsigned int> orig_patches;

  History::Update update;  
  for(unsigned int i = 0; i < fragin.pieces.size(); i++) {
    NxsPatch &patch = fragin.pieces[i];
    update.erased.push_back(patch.patch);
    orig_patches.insert(patch.patch);
  }
  
  vector<unsigned int> patch_remap;
  patch_remap.resize(fragout.pieces.size());
  
  for(unsigned int i = 0; i < fragout.pieces.size(); i++) {
    NxsPatch &patch = fragout.pieces[i];
    //TODO detect best parameter below for bordersize...
    unsigned int bordsize = 6 * patch.bord.size();
    if(bordsize > 65500) bordsize = 65499;
    unsigned int patch_idx = nexus.AddPatch(patch.vert.size(),
					    patch.face.size()/3,
					    bordsize);
    patch_levels.push_back(current_level);
    Entry &entry = nexus[patch_idx];
    entry.error = fragout.error;

    patch_remap[i] = patch_idx;
    chain.newfragments[patch.patch].insert(patch_idx);
    update.created.push_back(patch_idx);
  }
  
  //here i put for every patch all its new links.
  map<unsigned int, set<Link> > newlinks;

  for(unsigned int i = 0; i < fragout.pieces.size(); i++) {
    NxsPatch &outpatch = fragout.pieces[i];
    unsigned int patch_idx = patch_remap[i];
    
    for(unsigned int k = 0; k < outpatch.face.size(); k += 3) {
      assert(k+2 < outpatch.face.size());
      assert(outpatch.face[k]   != outpatch.face[k+1]);
      assert(outpatch.face[k]   != outpatch.face[k+2]);
      assert(outpatch.face[k+1] != outpatch.face[k+2]);
    }
    
    Patch &patch = nexus.GetPatch(patch_idx);
    memcpy(patch.FaceBegin(), &outpatch.face[0], 
	          outpatch.face.size() * sizeof(unsigned short));
    memcpy(patch.VertBegin(), &outpatch.vert[0], 
	          outpatch.vert.size() * sizeof(Point3f));
    
    Entry &entry = nexus[patch_idx];
    for(unsigned int v = 0; v < outpatch.vert.size(); v++) {
      entry.sphere.Add(outpatch.vert[v]);
      nexus.sphere.Add(outpatch.vert[v]);
    } 
    entry.sphere.Radius() *= 1.01;

    //remap internal borders
    for(unsigned int k = 0; k < outpatch.bord.size(); k++) {
      Link &link = outpatch.bord[k];
      if(link.end_patch >= (1<<31)) { //internal
	      link.end_patch = patch_remap[link.end_patch - (1<<31)];
	      assert(link.end_patch != patch_remap[i]);
	      newlinks[patch_idx].insert(link);
      } 
    } 
    //TODO not efficient!
    //processing external borders
    for(unsigned int l = 0; l < outpatch.bord.size(); l++) {
      Link &link = outpatch.bord[l];
      if(link.end_patch >= (1<<31)) continue;

      unsigned short &start_vert = link.start_vert;
      unsigned int &start_patch = patch_idx;

      //Adding vertical connection
      newlinks[patch_idx].insert(link);

      Link up(link.end_vert, link.start_vert, patch_idx);
      newlinks[link.end_patch].insert(up);
      
      assert(link.end_patch != patch_idx);
      Border &cborder = nexus.GetBorder(link.end_patch);
      for(unsigned int k = 0; k < cborder.Size(); k++) {
        //cerr << "Cborder.size: " << cborder.Size() << endl;
        //cerr << "K: " << k << endl;
	      Link &clink = cborder[k];
	      assert(!clink.IsNull());

	      if(clink.start_vert != link.end_vert) continue;
	      if(patch_levels[clink.end_patch] < current_level-1) continue;
	      //boy i want only the external links!
	      if(orig_patches.count(clink.end_patch)) continue;

    	  unsigned short &end_vert = clink.end_vert;
	      unsigned int &end_patch = clink.end_patch;
	    
	      //TODO FIX THIS!!!!
	      assert(clink.end_patch != start_patch);
	      assert(clink.end_patch != link.end_patch);

	      Link newlink;

	      newlink.start_vert = start_vert;
	      newlink.end_vert = end_vert;
	      newlink.end_patch = end_patch;
	
	      newlinks[patch_idx].insert(newlink);

	      newlink.start_vert = end_vert;
	      newlink.end_vert = start_vert;
	      newlink.end_patch = start_patch;
	
	      newlinks[end_patch].insert(newlink);
      }
    }
  }
  map<unsigned int, set<Link> >::iterator n;
  for(n = newlinks.begin(); n != newlinks.end(); n++) {
    set<Link>::iterator l;
    unsigned int patch = (*n).first;
    set<Link> &links = (*n).second;
    Border &border = nexus.GetBorder(patch);    
    unsigned int bstart = border.Size();
    nexus.borders.ResizeBorder(patch, border.Size() + links.size());                
    for(l = links.begin(); l != links.end(); l++) {
      Link link = *l;
      border[bstart++] = link;
    }
  }
  nexus.history.updates.push_back(update);
}

void ReverseHistory(vector<History::Update> &history) {
  vector<History::Update> revert = history;
  history.clear();
  for(int i = revert.size()-1; i >= 0; i--)
    history.push_back(revert[i]);    
  //std::reverse(history.begin(), history.end());
  vector<History::Update>::iterator i;
  for(i = history.begin(); i != history.end(); i++)
    swap((*i).erased, (*i).created);
}
