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
Revision 1.19  2004/10/22 10:37:32  ponchio
Split is now in fragment.

Revision 1.18  2004/10/21 13:40:16  ponchio
Debugging.

Revision 1.17  2004/10/21 12:22:21  ponchio
Small changes.

Revision 1.16  2004/10/19 01:23:02  ponchio
Daily backup (fragment...)

Revision 1.15  2004/10/15 11:41:03  ponchio
Tests and small changes.

Revision 1.14  2004/10/10 17:19:42  ponchio
Added compression and debugged.

Revision 1.13  2004/10/09 17:32:25  ponchio
Ram buffer option added (last)

Revision 1.12  2004/10/09 17:29:04  ponchio
Ram buffer option added (again)

Revision 1.11  2004/10/09 17:22:57  ponchio
Ram buffer option added.

Revision 1.10  2004/10/09 16:51:36  ponchio
Windows porting small changes.

Revision 1.9  2004/10/08 15:12:04  ponchio
Working version (maybe)

Revision 1.8  2004/10/06 16:40:47  ponchio
Fixed degenerate faces.

Revision 1.7  2004/10/04 16:49:54  ponchio
Daily backup. Preparing for compression.

Revision 1.6  2004/10/01 16:54:57  ponchio
Daily backup.

Revision 1.5  2004/09/30 00:27:42  ponchio
Lot of changes. Backup.

Revision 1.4  2004/09/21 00:53:23  ponchio
Lotsa changes.

Revision 1.3  2004/09/17 15:25:09  ponchio
First working (hopefully) release.

Revision 1.2  2004/09/16 14:25:16  ponchio
Backup. (lot of changes).

Revision 1.1  2004/08/26 18:03:48  ponchio
First draft.


****************************************************************************/

#ifdef WIN32
#include <wrap/system/getopt.h>
#else
#include <unistd.h>
#endif

#include <iostream>
using namespace std;

#include "crude.h"
#include "nexus.h"

#include "voronoichain.h"
#include "vert_remap.h"

#include "decimate.h"
#include "fragment.h"
#include "nxsbuild.h"
#include "watch.h"
using namespace vcg;
using namespace nxs;


void BuildFragment(Nexus &nexus, VoronoiPartition &part,
		   set<unsigned int> &patches, 
		   Nexus::Update &update,
		   Fragment &fragment);

void SaveFragment(Nexus &nexus, VoronoiChain &chain,
		  unsigned int level,
		  Nexus::Update &update,
		  Fragment &fragout);

void ReverseHistory(vector<Nexus::Update> &history);

void TestBorders(Nexus &nexus);
void TestPatches(Nexus &nexus);

int main(int argc, char *argv[]) {

  Decimation decimation = CLUSTER;
  unsigned int patch_size = 1000;
  unsigned int patch_threshold = 0xffffffff;
  unsigned int optimization_steps = 5;
  bool stop_after_remap = false;
  unsigned int max_level = 0xffffffff;
  float scaling = 0.5;
  unsigned int ram_buffer = 128000000;
  unsigned int chunk_size = 1024;

  int option;
  while((option = getopt(argc, argv, "f:t:l:s:d:ro:b:c:")) != EOF) {
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
    case 'r': stop_after_remap = true; break;
    case 'o': optimization_steps = atoi(optarg); break;
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

  //Test that there are still 2 arguments
  if(optind != argc - 2) {
    cerr << "Usage: " << argv[0] << " <crude> <output> [options]\n";
    cerr << "  Options:\n";
    cerr << " -f N: use N faces per patch (default 1000, max 32000)\n";
    cerr << " -t N: mini faces per patch (default 200)\n";
    cerr << " -l N: number of levels\n";
    cerr << " -s F: scaling factor (0 < F < 1, default 0.5)\n";
    cerr << " -o N: nomber of optimization steps\n";
    cerr << " -d <method>: decimation method: quadric, cluster. (default quadric)\n";
    cerr << " -b N: ram buffer size (in bytes)\n";
    cerr << " -r  : stop after remapping fase\n";
    return -1;
  }

  Crude crude;
  if(!crude.Load(argv[optind])) {
    cerr << "Could not open crude input: " << argv[optind] << endl;
    return -1;
  }
  
  if(patch_size > crude.vert.Size()/2) {
    cerr << "Patch size too big: " << patch_size << " * 2 > " << crude.vert.Size()
	 << endl;
    return -1;
  }

  string output = argv[optind+1];

  Nexus nexus;
  nexus.patches.SetRamBufferSize(ram_buffer);
  if(!nexus.Create(output, NXS_FACES, chunk_size)) {
    cerr << "Could not create nexus output: " << output << endl;
    return -1;
  }


  if(patch_threshold == 0xffffffff)
    patch_threshold = patch_size/4;
  VoronoiChain vchain(patch_size, patch_threshold);
  //  vchain.scaling = scaling;

  //Now building level 0 of the Nexus

  VFile<unsigned int> face_remap;
  if(!face_remap.Create(output + ".rmf")) {
    cerr << "Could not create remap files: " << output << ".rmf\n";
    return -1;
  }
  face_remap.Resize(crude.Faces());

  VertRemap vert_remap;
  if(!vert_remap.Create(output)) {
    cerr << "Could not create remap file: " << output << ".rmv and .rmb\n";
    return -1;
  }
  vert_remap.Resize(crude.Vertices());
  
  VFile<RemapLink> border_remap;
  if(!border_remap.Create(output + string(".tmp"))) {
    cerr << "Could not create temporary border remap file\n";
    return -1;
  }

  /* BUILDING FIRST LEVEL */

  //Remapping faces and vertices using level 0 and 1 of the chain
  cerr << "Remapping faces.\n";
  vector<unsigned int> patch_faces;
  
  vchain.RemapFaces(crude, face_remap, patch_faces, 
		    scaling, optimization_steps);

  cerr << "Remapping vertices.\n";
  vector<unsigned int> patch_verts;
  patch_verts.resize(patch_faces.size(), 0);
  RemapVertices(crude, vert_remap, face_remap, patch_verts);

  if(stop_after_remap) return 0;

  cerr << "Allocating space\n";
  //allocate chunks for patches and copy faces (absoklute indexing) into them.
  NexusAllocate(crude, nexus, face_remap, patch_faces, patch_verts);

  cerr << "Filling first level\n";
  //insert vertices and remap faces, prepare borders
  NexusFill(crude, nexus, vert_remap, border_remap);

  //  NexusFixBorder(nexus, border_remap);

  //filling history
  Nexus::Update update;
  for(unsigned int i = 0; i < nexus.index.size(); i++) 
    update.created.push_back(i);
  nexus.history.push_back(update); 
 
  //unify vertices otherwise you may get cracks.
  nexus.Unify();
  nexus.patches.FlushAll();

  /* BUILDING OTHER LEVELS */

  Report report;

  unsigned int oldoffset = 0;
  for(unsigned int level = 1; level < max_level; level++) {
    cerr << "Level: " << level << endl;

    unsigned int newoffset = nexus.index.size();
    vchain.BuildLevel(nexus, oldoffset, scaling, optimization_steps);
    
    report.Init(vchain.oldfragments.size(), 1);
    unsigned int fcount = 0;
    vector<Nexus::Update> level_history;
    map<unsigned int, set<unsigned int> >::iterator fragment;
    for(fragment = vchain.oldfragments.begin(); 
	fragment != vchain.oldfragments.end(); fragment++) {
      report.Step(fcount++);

      

      Fragment fragin;
      BuildFragment(nexus, vchain.levels[level+1], (*fragment).second, 
		    update, fragin);




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
      Split(fragout, newvert, newface, newbord, vchain.levels[level+1]); 





      SaveFragment(nexus, vchain, level, update, fragout);
      level_history.push_back(update);
    }
    report.Finish();

    for(unsigned int i = 0; i < level_history.size(); i++)
      nexus.history.push_back(level_history[i]);

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
  nexus.history.push_back(update);
  ReverseHistory(nexus.history);

  //  TestBorders(nexus);
  nexus.Close();

  //TODO remove vert_remap, face_remap, border_remap
  return 0;
}

void BuildFragment(Nexus &nexus, VoronoiPartition &part,
		   set<unsigned int> &patches, 
		   Nexus::Update &update,
		   Fragment &fragment) {
  set<unsigned int>::iterator f;
  for(f = patches.begin(); f != patches.end(); f++) {
    fragment.pieces.push_back(NxsPatch());
    NxsPatch &nxs = fragment.pieces.back();
    nxs.patch = *f;

    Patch &patch = nexus.GetPatch(*f);
    Border border = nexus.GetBorder(*f);

    nxs.vert.resize(patch.nv);
    nxs.face.resize(patch.nf * 3);
    memcpy(&*nxs.vert.begin(), patch.VertBegin(), patch.nv * sizeof(Point3f));
    memcpy(&*nxs.face.begin(), patch.FaceBegin(), patch.nf * 3*sizeof(short));
    for(unsigned int i = 0; i < border.Size(); i++) {
      Link &link = border[i];
      if(!link.IsNull())
	nxs.bord.push_back(link);
    }
  }

  update.created.clear();
  update.erased.clear();
  
  set<unsigned int> &fcells = patches;
  set<unsigned int>::iterator s;
  for(s = fcells.begin(); s != fcells.end(); s++) {
    update.erased.push_back(*s);
  }

  //copy all seeds! //TODO copy only closest ones
  for(unsigned int i = 0; i < part.size(); i++) {
    fragment.seeds.push_back(part[i]);
    fragment.seeds_id.push_back(i);
  }
}


void SaveFragment(Nexus &nexus, VoronoiChain &chain,
		  unsigned int level,
		  Nexus::Update &update, 
		  Fragment &fragout) {
  
  vector<unsigned int> patch_remap;
  patch_remap.resize(fragout.pieces.size());

  for(unsigned int i = 0; i < fragout.pieces.size(); i++) {
    NxsPatch &patch = fragout.pieces[i];
    //TODO detect best parameter below for bordersize...
    unsigned int patch_idx = nexus.AddPatch(patch.vert.size(),
					    patch.face.size()/3,
					    6 * patch.bord.size());
    Nexus::PatchInfo &entry = nexus.index[patch_idx];
    entry.error = fragout.error;

    patch_remap[i] = patch_idx;
    chain.newfragments[patch.patch].insert(patch_idx);
    update.created.push_back(patch_idx);
  }

  for(unsigned int i = 0; i < fragout.pieces.size(); i++) {
    NxsPatch &outpatch = fragout.pieces[i];
    unsigned int patch_idx = patch_remap[i];
    
    Patch &patch = nexus.GetPatch(patch_idx);
    memcpy(patch.FaceBegin(), &outpatch.face[0], 
	   outpatch.face.size() * sizeof(unsigned short));
    memcpy(patch.VertBegin(), &outpatch.vert[0], 
	   outpatch.vert.size() * sizeof(Point3f));
    
    Nexus::PatchInfo &entry = nexus.index[patch_idx];
    for(unsigned int v = 0; v < outpatch.vert.size(); v++) {
      entry.sphere.Add(outpatch.vert[v]);
      nexus.sphere.Add(outpatch.vert[v]);
    } 
    //remap borders
    for(unsigned int i = 0; i < outpatch.bord.size(); i++) {
      Link &link = outpatch.bord[i];
      if(link.end_patch >= (1<<31)) //internal
	link.end_patch = patch_remap[link.end_patch - (1<<31)];
      else { //if external add the reverse border
	Border rborder = nexus.GetBorder(link.end_patch);

	unsigned int pos = rborder.Size();
	if(nexus.borders.ResizeBorder(link.end_patch, pos+1)) {
	  rborder = nexus.GetBorder(link.end_patch);
	}
      
	assert(rborder.Size() < rborder.Available());
	assert(rborder.Available() > pos);

	Link newlink; 
	newlink.start_vert = link.end_vert;
	newlink.end_vert = link.start_vert;
	newlink.end_patch = patch_idx;
	rborder[pos] = newlink;
      }
    }
    Border border = nexus.GetBorder(patch_idx);
    assert(border.Available() >= outpatch.bord.size());
    if(nexus.borders.ResizeBorder(patch_idx, outpatch.bord.size())) {
      border = nexus.GetBorder(patch_idx);
    }
    memcpy(&(border[0]), &(outpatch.bord[0]), 
	   outpatch.bord.size() * sizeof(Link)); 
  }
}

void ReverseHistory(vector<Nexus::Update> &history) {
  reverse(history.begin(), history.end());
  vector<Nexus::Update>::iterator i;
  for(i = history.begin(); i != history.end(); i++)
    swap((*i).erased, (*i).created);
}

void TestPatches(Nexus &nexus) {
  cerr << "TESTING PATCHES!!!!" << endl;
  for(unsigned int p = 0; p  < nexus.index.size(); p++) {
    Patch &patch = nexus.GetPatch(p);
    for(unsigned int i = 0; i < patch.nf; i++)
      for(int k = 0; k < 3; k++)
        if(patch.Face(i)[k] >= patch.nv) {
          cerr << "Totface: " << patch.nf << endl;
          cerr << "Totvert: " << patch.nv << endl;
          cerr << "Face: " << i << endl;
          cerr << "Val:  " << patch.Face(i)[k] << endl;
          exit(0);
        }
  }
}
void TestBorders(Nexus &nexus) {
  //check border correctnes
  nexus.borders.Flush();
  for(unsigned int i = 0; i < nexus.index.size(); i++) {
    Border border = nexus.GetBorder(i);
    for(unsigned int k = 0; k < border.Size(); k++) {
      Link &link = border[k];
      if(link.IsNull()) continue;
      if(link.end_patch >= nexus.index.size()) {
	cerr << "Patch: " << i << endl;
	cerr << "Bsize: " << border.Size() << endl;
	cerr << "Bava: " << border.Available() << endl;
	cerr << "K: " << k << endl;
	cerr << "end: " << link.end_patch << endl;
      }
      assert(link.end_patch < nexus.index.size());
    }
  }
}
