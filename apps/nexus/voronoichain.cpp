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
Revision 1.1  2004/08/26 18:03:47  ponchio
First draft.


****************************************************************************/

#include <iostream>

#include "voronoichain.h"

using namespace std;
using namespace vcg;
using namespace nxs;

void VoronoiChain::Initialize(unsigned int psize, unsigned int pthreshold) {
  patch_size = psize;
  patch_threshold = pthreshold;
}

void VoronoiChain::Init(Crude &crude) {
  box = crude.GetBox();
  radius = VoronoiPartition::OptimalRadius(crude, patch_size);
  float radius0 = radius;
  float radius1 = radius * 1.4;
  //  float radius1 = VoronoiPartition::OptimalRadius(crude, patch_size*2);

  levels.push_back(VoronoiPartition());
  levels.push_back(VoronoiPartition());
  VoronoiPartition &part0 = levels[0];
  VoronoiPartition &part1 = levels[1];
  part0.Init(crude.GetBox());
  part1.Init(crude.GetBox());

  VFile<Point3f>::iterator iter;
  for(iter = crude.vert.Begin(); iter != crude.vert.End(); ++iter) { 
    Point3f &v = *iter;
    unsigned int target_patch;

    float dist = part0.Closest(v, target_patch);
    if(dist >= radius0 || dist == -1)
      part0.Add(v, radius0);

    dist = part1.Closest(v, target_patch);
    if(dist >= radius1 || dist == -1)
      part1.Add(v, radius1);
  }

  //here goes some optimization pass.
}

unsigned int VoronoiChain::Locate(unsigned int level, 
				  const vcg::Point3f &p) {
  return levels[level].Locate(p);
  /*  assert(levels.size() > level+1);
  unsigned int fine = levels[level].Locate(p);
  unsigned int coarse = levels[level+1].Locate(p);
  return fine + coarse * levels[level].size();*/
}

void VoronoiChain::RemapFaces(Crude &crude, 
			      VFile<unsigned int> &face_remap,
			      vector<unsigned int> &patch_faces) {
  
  Init(crude);

  //TODO: improve quality of patches and implement threshold.
  typedef  map<pair<unsigned int, unsigned int>, unsigned int> FragIndex;

  //  map<pair<unsigned int, unsigned int>, unsigned int> patches;
  FragIndex patches;

  unsigned int totpatches = 0;

  Point3f bari;
  for(unsigned int i = 0; i < crude.Faces(); i++) {
    bari = crude.GetBari(i);
    //    unsigned int patch = Locate(0, bari);
    unsigned int fine = Locate(0, bari);
    unsigned int coarse = Locate(1, bari);

    unsigned int patch;
    
    if(!patches.count(make_pair(coarse, fine))) {
      patch = totpatches;
      patches[make_pair(coarse, fine)] = totpatches++;
    } else
      patch = patches[make_pair(coarse, fine)];

    face_remap[i] = patch;
    if(patch_faces.size() <= patch) 
      patch_faces.resize(patch+1, 0);
    patch_faces[patch]++;
  }



  //prune faces (now only 0 faces);
  unsigned int tot_patches = 0;
  vector<int> patch_remap;
  for(unsigned int i = 0; i < patch_faces.size(); i++) {
    //if below threshold (and can join faces)
    if(patch_faces[i] == 0)
      patch_remap.push_back(-1);
    else
      patch_remap.push_back(tot_patches++);
  }


  //building fragments
  FragIndex::iterator f;
  for(f = patches.begin(); f != patches.end(); f++) {
    unsigned int coarse = (*f).first.first;
    unsigned int fine = (*f).first.second;
    unsigned int patch = (*f).second;
    oldfragments[coarse].insert(patch_remap[patch]);
  }  

  //remapping faces
  for(unsigned int i = 0; i < face_remap.Size(); i++) {
    unsigned int patch = face_remap[i];
    assert(patch != 0xffffffff);
    assert(patch_remap[patch] != -1);
    face_remap[i] = patch_remap[patch];
  }

  //remapping patch_faces
  for(unsigned int i = 0; i < patch_faces.size(); i++) {
    assert(patch_remap[i] <= (int)i);
    if(patch_remap[i] != -1) {
      assert(patch_faces[i] > 0);
      patch_faces[patch_remap[i]] = patch_faces[i];
    }
  }

  patch_faces.resize(tot_patches);
}

void VoronoiChain::BuildLevel(Nexus &nexus, unsigned int offset) {
  unsigned int target_faces =  (int)(patch_size * 
				     pow(scaling, (float)levels.size()));

  float rad = radius * pow(sqrt(1/scaling), (float)levels.size());

  levels.push_back(VoronoiPartition());
  VoronoiPartition &part = levels[levels.size()-1];
  part.Init(box);

  for(unsigned int idx = offset; idx < nexus.index.size(); idx++) {
    Patch patch = nexus.GetPatch(idx);
    for(unsigned int i = 0; i < patch.VertSize(); i++) {
      Point3f &v = patch.Vert(i);
      unsigned int target_patch;

      float dist = part.Closest(v, target_patch);
      if(dist >= rad || dist == -1)
	part.Add(v, rad);
    }
  }
  cerr << "radius: " << rad << " ... cells: " << part.size() << endl;
  newfragments.clear();
  //TODO add some optimization
}
