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
Revision 1.3  2004/09/17 15:25:09  ponchio
First working (hopefully) release.

Revision 1.2  2004/09/16 14:25:16  ponchio
Backup. (lot of changes).

Revision 1.1  2004/08/26 18:03:47  ponchio
First draft.


****************************************************************************/

#include <iostream>

#include "voronoichain.h"

using namespace std;
using namespace vcg;
using namespace nxs;


void VoronoiChain::Init(Crude &crude, float scaling, int steps) {
  box = crude.GetBox();
  box.Offset(box.max - box.min);
  //how many cells do i need?
  cerr << "scaling: " << scaling << endl;
  unsigned int f_cells = crude.Faces() / mean_size;
  unsigned int c_cells = (unsigned int)(scaling * f_cells);

  cerr << "mean size: " << mean_size << endl;
  cerr << "f cells: " << f_cells << endl;
  cerr << "c_cells: " << c_cells << endl;

  levels.push_back(VoronoiPartition());
  levels.push_back(VoronoiPartition());
  VoronoiPartition &fine = levels[0];
  VoronoiPartition &coarse = levels[1];

  fine.Init(box);
  coarse.Init(box);

  srand(0);
  
  vector<Seed> fine_seeds;
  vector<Seed> coarse_seeds;

  float fine_vmean = mean_size/2;
  float coarse_vmean = (mean_size/scaling)/2;
  for(unsigned int i = 0; i < crude.Vertices(); i++) {
    int f = (int)(fine_vmean*rand()/(RAND_MAX + 1.0));
    int c = (int)(coarse_vmean *rand()/(RAND_MAX + 1.0));
    if(f == 1) {
      Point3f &point = crude.GetVertex(i);
      fine_seeds.push_back(Seed(point, 1));
    }
    if(c == 1) {
      Point3f &point = crude.GetVertex(i);
      coarse_seeds.push_back(Seed(point, 1));
    }
  }
  cerr << "fine_seeds.size: " << fine_seeds.size() << endl;
  cerr << "coarse_seeds.size: " << coarse_seeds.size() << endl;
  fine.all_seeds = fine_seeds;
  coarse.all_seeds =  coarse_seeds;
  fine.reload();
  coarse.reload();

  //here goes some optimization pass.
  //Fine optimization.
  vector<Point3f> fcentroids;
  vector<unsigned int> fcount;
  for(unsigned int i = 0; i < steps; i++) {
    cerr << "Optimization step 0: " << i << "/" << steps << endl;
    fcentroids.clear();
    fcount.clear();
    fcentroids.resize(fine.size(), Point3f(0, 0, 0));
    fcount.resize(fine.size(), 0);
    
    for(unsigned int v = 0; v < crude.Vertices(); v++) {
      unsigned int ftarget;
      float dist = fine.Closest(crude.vert[v], ftarget);
      assert(ftarget != -1);
      fcentroids[ftarget] += crude.vert[v];
      fcount[ftarget]++;
    }
    for(unsigned int v = 0; v < fine.size(); v++) {
      assert(fcount[v] != 0);
       
      fine[v].p = fcentroids[v]/fcount[v];
      //0.3 is related to the fact is doubled the box size.
      fine[v].weight = pow(fcount[v]/fine_vmean, 0.3f);
      //      fine.bbox.Add(fine[v].p);
    }
    //    fine.Init(fine.bbox);
    fine.reload();
  }    

  //Coarse optimization
  vector< map<unsigned int, Point3f> > ccentroids;
  vector< map<unsigned int, unsigned int> > ccount;

  for(unsigned int i = 0; i < steps; i++) {
    cerr << "Optimization step 1: " << i << "/" << steps << endl;
    ccentroids.clear();
    ccount.clear();
    ccentroids.resize(coarse.size());
    ccount.resize(coarse.size());

    for(unsigned int v = 0; v < crude.Vertices(); v++) {
      unsigned int ftarget;
      float dist = fine.Closest(crude.vert[v], ftarget);
      assert(ftarget != -1);

      unsigned int ctarget;
      dist = coarse.Closest(crude.vert[v], ctarget);
      assert(ctarget != -1);

      map<unsigned int, Point3f> &centroids = ccentroids[ctarget];
      map<unsigned int, unsigned int> &count = ccount[ctarget];

      if(!centroids.count(ftarget))
	centroids[ftarget]= Point3f(0, 0, 0);
	
      if(!count.count(ftarget))
	count[ftarget] = 0;

      centroids[ftarget] += crude.vert[v];
      count[ftarget]++;
    }
    
    for(unsigned int v = 0; v < coarse.size(); v++) {

      map<unsigned int, Point3f> &centroids = ccentroids[v];
      map<unsigned int, unsigned int> &count = ccount[v];
      

      coarse[v].p = Point3f(0, 0, 0);
      float weight = 0;
      unsigned int tot_size =0;
      map<unsigned int, Point3f>::iterator k;
      for(k = centroids.begin();k != centroids.end(); k++) {
	unsigned int size = count[(*k).first];
	tot_size += size;
	//coarse[v].p += (*k).second / (size * size);
	//weight += 1/(float)size;
	coarse[v].p += (*k).second / size;
	weight += 1;
	//	coarse[v].p += (*k).second;
	//	weight += size;
      }
      assert(weight > 0);
      coarse[v].p /= weight;
      //TODO find a solution
      //      coarse[v].weight = pow(tot_size/coarse_vmean, 0.25f);

    }
    coarse.reload();
  }
}

unsigned int VoronoiChain::Locate(unsigned int level, 
				  const vcg::Point3f &p) {
  return levels[level].Locate(p);
  /*  assert(levels.size() > level+1);
  unsigned int fine = levels[level].Locate(p);
  unsigned int coarse = levels[level+1].Locate(p);
  return fine + coarse * levels[level].size();*/
}

void VoronoiChain::RemapFaces(Crude &crude, 			      VFile<unsigned int> &face_remap,
			      vector<unsigned int> &patch_faces,
			      float scaling, int steps) {
  
  Init(crude, scaling, steps);

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

    //BEWARE unkomment this!
    face_remap[i] = patch;
    //face_remap[i] = fine;
    //    face_remap[i] = coarse;
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

void VoronoiChain::BuildLevel(Nexus &nexus, unsigned int offset, 
			      float scaling, int steps) {
  //  unsigned int target_faces =  (int)(mean_size * 
  //				     pow(0.5f, (float)levels.size()));

  unsigned int totface = 0;
  for(unsigned int idx = offset; idx < nexus.index.size(); idx++)
    totface += nexus.index[idx].nface;
  
  levels.push_back(VoronoiPartition());
  VoronoiPartition &coarse = levels[levels.size()-1];
  VoronoiPartition &fine = levels[levels.size()-2];

  coarse.Init(box);
  
  fine.reload();

  srand(0);
  float coarse_vmean = (totface/2)/(fine.size() * scaling);
  
  cerr << "initing random seeds\n";

  vector<Seed> coarse_seeds;
  
  for(unsigned int idx = offset; idx < nexus.index.size(); idx++) {
    Patch patch = nexus.GetPatch(idx);
    for(unsigned int i = 0; i < patch.nv; i++) {
      int c = (int)(coarse_vmean*rand()/(RAND_MAX + 1.0));
      if(c == 1) {
	Point3f &v = patch.Vert(i);
	coarse_seeds.push_back(v);
      }
    }
  }
  if(coarse_seeds.size() == 0)
    coarse_seeds.push_back(Point3f(0, 0, 0));
  coarse.all_seeds =  coarse_seeds;
  coarse.reload();

//Coarse optimization
  vector< map<unsigned int, Point3f> > ccentroids;
  vector< map<unsigned int, unsigned int> > ccount;

  for(unsigned int step = 0; step < steps; step++) {
    cerr << "Optimization step " << levels.size()-1 << ":" 
	 << step << "/" << steps << endl;
    ccentroids.clear();
    ccount.clear();
    ccentroids.resize(coarse.size());
    ccount.resize(coarse.size());

    for(unsigned int idx = offset; idx < nexus.index.size(); idx++) {
      Patch patch = nexus.GetPatch(idx);

      for(unsigned int i = 0; i < patch.nv; i++) {
	Point3f &v = patch.Vert(i);
	unsigned int ftarget;
	float dist = fine.Closest(Point3f(1,1,1), ftarget);

	dist = fine.Closest(v, ftarget);
	assert(ftarget != -1);

	unsigned int ctarget;
	dist = coarse.Closest(v, ctarget);
	assert(ctarget != -1);

	map<unsigned int, Point3f> &centroids = ccentroids[ctarget];
	map<unsigned int, unsigned int> &count = ccount[ctarget];

	if(!centroids.count(ftarget))
	  centroids[ftarget]= Point3f(0, 0, 0);
	
	if(!count.count(ftarget))
	  count[ftarget] = 0;
	
	centroids[ftarget] += v;
	count[ftarget]++;
      }
    }

    cerr << "recentring" << endl;
    for(unsigned int v = 0; v < coarse.size(); v++) {
      
      map<unsigned int, Point3f> &centroids = ccentroids[v];
      map<unsigned int, unsigned int> &count = ccount[v];
      
      coarse[v].p = Point3f(0, 0, 0);
      float weight = 0;
      unsigned int tot_size =0;
      map<unsigned int, Point3f>::iterator k;
      for(k = centroids.begin();k != centroids.end(); k++) {
	unsigned int size = count[(*k).first];
	tot_size += size;
	coarse[v].p += (*k).second / size;
	weight += 1;
      }
      assert(weight > 0);
      coarse[v].p /= weight;
      //TODO find a solution!
      //      coarse[v].weight = pow(tot_size/coarse_vmean, 0.25f);
    }
    coarse.reload();
  }
  newfragments.clear();
  //TODO add some optimization
}
