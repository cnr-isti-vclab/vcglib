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
Revision 1.14  2004/10/19 04:23:29  ponchio
*** empty log message ***

Revision 1.13  2004/10/15 16:45:27  ponchio
Vbo added.

Revision 1.12  2004/10/15 11:41:03  ponchio
Tests and small changes.

Revision 1.11  2004/10/10 17:19:42  ponchio
Added compression and debugged.

Revision 1.10  2004/10/09 14:46:47  ponchio
Windows porting small changes.

Revision 1.9  2004/10/08 15:12:04  ponchio
Working version (maybe)

Revision 1.8  2004/10/04 16:49:54  ponchio
Daily backup. Preparing for compression.

Revision 1.7  2004/10/01 16:54:57  ponchio
Daily backup.

Revision 1.6  2004/09/30 00:27:42  ponchio
Lot of changes. Backup.

Revision 1.5  2004/09/28 10:26:07  ponchio
Voronoi partition changes.

Revision 1.4  2004/09/21 00:53:23  ponchio
Lotsa changes.

Revision 1.3  2004/09/17 15:25:09  ponchio
First working (hopefully) release.

Revision 1.2  2004/09/16 14:25:16  ponchio
Backup. (lot of changes).

Revision 1.1  2004/08/26 18:03:47  ponchio
First draft.


****************************************************************************/

#include <iostream>

#include "voronoichain.h"
#include "watch.h"

using namespace std;
using namespace vcg;
using namespace nxs;


void VoronoiChain::Init(Crude &crude, float scaling, int steps) {
  unsigned int f_cells = crude.Faces() / mean_size;
  unsigned int c_cells = (unsigned int)(scaling * f_cells);

  //cerr << "mean size: " << mean_size << endl;
  //cerr << "f cells: " << f_cells << endl;
  //cerr << "c_cells: " << c_cells << endl;

  levels.push_back(VoronoiPartition());
  levels.push_back(VoronoiPartition());
  VoronoiPartition &fine = levels[0];
  VoronoiPartition &coarse = levels[1];
  fine.SetBox(crude.GetBox());
  coarse.SetBox(crude.GetBox());

  srand(0);
  
  float fine_vmean = mean_size/2.0f;
  float coarse_vmean = (mean_size/scaling)/2;
  for(unsigned int i = 0; i < crude.Vertices(); i++) {
    int f = (int)(fine_vmean*rand()/(RAND_MAX + 1.0));
    int c = (int)(coarse_vmean *rand()/(RAND_MAX + 1.0));
    if(f == 1) {
      Point3f &point = crude.GetVertex(i);
      fine.push_back(Seed(point, 1));
    }
    if(c == 1) {
      Point3f &point = crude.GetVertex(i);
      coarse.push_back(Seed(point, 1));
    }
  }
  //TODO! Check for duplicates (use the closest :P)
  //cerr << "fine_seeds.size: " << fine.size() << endl;
  //cerr << "coarse_seeds.size: " << coarse.size() << endl;
  fine.Init();
  coarse.Init();

//here goes some optimization pass.
   //Fine optimization.
  Report report;
  vector<Point3f> fcentroids;
  vector<unsigned int> fcount;
  for(int i = 0; i < steps; i++) {
    cerr << "Optimization step: " << i+1 << "/" << steps << endl;
    fcentroids.clear();
    fcount.clear();
    fcentroids.resize(fine.size(), Point3f(0, 0, 0));
    fcount.resize(fine.size(), 0);
    
    report.Init(crude.Vertices());
    for(unsigned int v = 0; v < crude.Vertices(); v++) {
      if(v & 0xffff) report.Step(v);
      unsigned int ftarget;
      float dist = fine.Closest(crude.vert[v], ftarget);
      assert(ftarget != -1);
      fcentroids[ftarget] += crude.vert[v];
      fcount[ftarget]++;
    }
    for(unsigned int v = 0; v < fine.size(); v++) {
      if(fcount[v] == 0) continue;
      fine[v].p = fcentroids[v]/(float)fcount[v];
      fine[v].weight = (float)pow(fcount[v]/(float)fine_vmean, 0.2f);
    }
    fine.Init();
  }     

  //remove small or zero patches.
  vector<Seed> seeds;
  for(unsigned int i = 0; i < fine.size(); i++) {
    if(fcount[i] > min_size)
      seeds.push_back(fine[i]);
  }
  fine.clear();
  for(unsigned int i = 0; i < seeds.size(); i++)
    fine.push_back(seeds[i]);

  if(fine.size() == 0) fine.push_back(Point3f(0,0,0));
  fine.Init();

  //here goes some optimization pass.
  //Coarse optimization.
  vector<Point3f> ccentroids;
  vector<unsigned int> ccount;
  vector<float> radius;
  for(int i = 0; i < steps; i++) {
    cerr << "Optimization step: " << i+1 << "/" << steps << endl;
    ccentroids.clear();
    ccount.clear();
    ccentroids.resize(coarse.size(), Point3f(0, 0, 0));
    ccount.resize(coarse.size(), 0);
    radius.resize(coarse.size(), 0);
    
    report.Init(crude.Vertices());
    for(unsigned int v = 0; v < crude.Vertices(); v++) {
      if(v & 0xffff) report.Step(v);
      unsigned int ctarget = 0xffffffff;
      float dist = coarse.Closest(crude.vert[v], ctarget);
      assert(ctarget != 0xffffffff);
      ccentroids[ctarget] += crude.vert[v];
      ccount[ctarget]++;
      if(dist > radius[ctarget]) radius[ctarget] = dist;
    }
    for(unsigned int v = 0; v < coarse.size(); v++) {
      if(ccount[v] == 0) continue;
      
      coarse[v].p = ccentroids[v]/(float)ccount[v];
      coarse[v].weight = (float)pow(ccount[v]/(float)coarse_vmean, 0.2f);

      if(radius[v] == 0) continue;

      //repel from fine seeds
      //      coarse[v].p = fine.FindBorder(coarse[v].p, radius[v]);
    }
    coarse.Init();
  }    


  //remove small or zero patches.
  seeds.clear();
  for(unsigned int i = 0; i < coarse.size(); i++) {
    if(ccount[i] > (int)min_size)
      seeds.push_back(coarse[i]);
  }
  coarse.clear();
  for(unsigned int i = 0; i < seeds.size(); i++)
    coarse.push_back(seeds[i]);

  if(coarse.size() == 0) coarse.push_back(Point3f(0,0,0));
  coarse.Init();


  //Coarse optimization
  /*  vector< map<unsigned int, Point3f> > ccentroids;
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
    coarse.Init();
    }*/
}

unsigned int VoronoiChain::Locate(unsigned int level, 
				  const vcg::Point3f &p) {
  return levels[level].Locate(p);
  /*  assert(levels.size() > level+1);
  unsigned int fine = levels[level].Locate(p);
  unsigned int coarse = levels[level+1].Locate(p);
  return fine + coarse * levels[level].size();*/
}

//TODO move this to nxsbuild
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

    face_remap[i] = patch;
    //face_remap[i] = fine;

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
  unsigned int totface = 0;
  unsigned int totvert = 0;
  for(unsigned int idx = offset; idx < nexus.index.size(); idx++) {
    totface += nexus.index[idx].nface;
    totvert += nexus.index[idx].nvert;
  }

  
  levels.push_back(VoronoiPartition());
  VoronoiPartition &coarse = levels[levels.size()-1];
  VoronoiPartition &fine = levels[levels.size()-2];
  coarse.SetBox(fine.box);
  fine.Init();
  
  unsigned int tot_coarse = (unsigned int)(fine.size() * scaling);
  
  //TODO this method for selecting the seeds is ugly!
  float ratio = tot_coarse/(float)(nexus.index.size() - offset);
  float cratio = 0;
  for(unsigned int idx = offset; idx < nexus.index.size(); idx++) {
    cratio += ratio;
    if(cratio > 1) {
      Patch patch = nexus.GetPatch(idx);
      Point3f &v = patch.Vert(0);
      coarse.push_back(v);
      cratio -= 1;
    }
  }

  if(coarse.size() == 0) {
    Patch patch = nexus.GetPatch(0);
    coarse.push_back(patch.Vert(0));
  }

  float coarse_vmean = totvert/(float)coarse.size();

  coarse.Init();
  cerr << "Coarse size: " << coarse.size() << endl;
  cerr << "Coarse mean: " << coarse_vmean << " mean_size: " << mean_size << endl;

  Report report;
//here goes some optimization pass.
  //Coarse optimization.
  vector<Point3f> ccentroids;
  vector<unsigned int> ccount;
  vector<float> radius;
  for(int i = 0; i < steps; i++) {
    cerr << "Optimization step: " << i+1 << "/" << steps << endl;
    ccentroids.clear();
    ccount.clear();
    ccentroids.resize(coarse.size(), Point3f(0, 0, 0));
    ccount.resize(coarse.size(), 0);
    radius.resize(coarse.size(), 0);

    report.Init(nexus.index.size());
    for(unsigned int idx = offset; idx < nexus.index.size(); idx++) {
      report.Step(idx);
      Patch patch = nexus.GetPatch(idx);
      for(unsigned int i = 0; i < patch.nv; i++) {

	unsigned int ctarget = 0xffffffff;
	float dist = coarse.Closest(patch.Vert(i), ctarget);
	assert(ctarget != 0xffffffff);
	ccentroids[ctarget] += patch.Vert(i);
	ccount[ctarget]++;
	if(dist > radius[ctarget]) radius[ctarget] = dist;
      }
    }
    for(unsigned int v = 0; v < coarse.size(); v++) {
      if(ccount[v] == 0) continue;
       
      coarse[v].p = ccentroids[v]/(float)ccount[v];
      //0.3 is related to the fact is doubled the box size.
      coarse[v].weight = (float)pow(ccount[v]/(float)coarse_vmean, 0.2f);

      if(radius[v] == 0) continue;
      //coarse[v].p = fine.FindBorder(coarse[v].p, radius[v]);
    }
    coarse.Init();
  }    
  vector<Seed> seeds;
  for(unsigned int i = 0; i < coarse.size(); i++) {
    if(ccount[i] > (int)min_size)
      seeds.push_back(coarse[i]);
  }

  coarse.clear();
  for(unsigned int i = 0; i < seeds.size(); i++)
    coarse.push_back(seeds[i]);

  if(coarse.size() == 0) coarse.push_back(Point3f(0,0,0));
  coarse.Init();

//Coarse optimization
/*  vector< map<unsigned int, Point3f> > ccentroids;
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
	float dist = fine.Closest(v, ftarget);

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
    coarse.Init();
    }*/
  newfragments.clear();
  //TODO add some optimization
}
