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
Revision 1.19  2004/11/03 16:31:38  ponchio
Trying to fix big patches.

Revision 1.18  2004/10/30 20:17:03  ponchio
Fixed big patches problem.

Revision 1.17  2004/10/29 16:33:29  ponchio
Trying to fix big patches.

Revision 1.16  2004/10/22 14:31:56  ponchio
Some controls added.

Revision 1.15  2004/10/21 12:22:21  ponchio
Small changes.

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


void print(Point3f p) {
  cerr << p[0] << " " << p[1] << " " << p[2] << endl;
}

//return first non zero distance point.
float getClosest(const Point3f &seed, VoronoiPartition &part) {
  vector<int> nears;
  vector<float> dists;
  float dist = 0;
  int count = 1;
  while(dist == 0) {
    if(count > part.size()) {
      cerr << "This should never happen!!!!\n";
      exit(0);
    }
    part.Closest(seed, count, nears, dists);
    for(int k = 0; k < count; k++) {
      int c = nears[k];
      assert(c >= 0);
      assert(c < part.size());      
      if(dists[k] > 0 && (dist == 0 || dists[k] < dist)) {
        dist = dists[k];                              
      }
    }
    count++;
  }
  return sqrt(dist);
}

int getBest(const Point3f &seed, VoronoiPartition &part, 
	    vector<bool> &mark,
	    vector<unsigned int> &counts) {

  vector<int> nears;
  vector<float> dist;
  int nnear = 7;
  if(part.size() < 7) nnear = part.size()/2;
  if(!nnear) return -1;
  
  part.Closest(seed, nnear, nears, dist);
  int best = -1;
  int bestcount = -1;
  int bestdist = -1;
  for(int k = 0; k < nnear; k++) {
    int c = nears[k];
    assert(c >= 0);
    assert(c < part.size());    if(mark[c]) continue;
    if(part[c] == seed) continue;
    if(bestcount < 0 || 
       (counts[c] < bestcount)) {
      best = c;
      bestcount = counts[c];
    }
    /*if(bestdist < 0 ||
       Distance(seed, part[c]) < bestdist) {
      best = c;
      bestdist = Distance(seed, part[c]);
      }*/
  }
  return best;
}

//return false if still not ok
bool VoronoiChain::Optimize(int mean, VoronoiPartition &part, 
			    vector<Point3f> &centroids,
			    vector<unsigned int> &counts, 
			    bool join) {



    //remove small or really big patches.
    unsigned int failed = 0;
    vector<Point3f> seeds;
    vector<bool> mark;
    mark.resize(part.size(), false);

    //first pass we check only big ones
    for(unsigned int i = 0; i < part.size(); i++) {
      if(counts[i] > max_size || counts[i] > 2 * mean) {
	      failed++;
	      cerr << "Failed> " << counts[i] << endl;
	      float radius= getClosest(part[i], part);
        cerr << "RADIUS: " << radius << endl;
        if(radius == 0) {
          cerr << "Radius zero???\n";
          exit(0);
        }
	      radius /= 3;
	      if(radius < 0) continue;
	      seeds.push_back(part[i] + Point3f(1, 0, 0) * radius);
	      seeds.push_back(part[i] + Point3f(0, 1, 0) * radius);
	      seeds.push_back(part[i] + Point3f(0, 0, 1) * radius);
      
	      seeds.push_back(part[i] - Point3f(1, 0, 0) * radius);
	      seeds.push_back(part[i] - Point3f(0, 1, 0) * radius);
	      seeds.push_back(part[i] - Point3f(0, 0, 1) * radius);
	      mark[i];
      }
    }

    cerr << "Join now!" << endl;
    for(unsigned int i = 0; i < part.size(); i++) {
      if(mark[i]) continue;
      if(join && counts[i] < min_size) {
        failed++;
	      int best = getBest(part[i], part, mark, counts);
        if(best < 0) {
          cerr << "Best not found! how strange!\n";
          continue;
        }
        if(best >= part.size()) {
          cerr << "Invalid best!!!\n";
          exit(0);
        }
	      assert(mark[best] == false);
	      mark[best] = true;
	      mark[i] = true;
	      seeds.push_back((part[i] + part[best])/2);
      }
    }

    for(unsigned int i = 0; i < part.size(); i++) {
      if(mark[i]) continue;
      if(join) {
        if(counts[i] < min_size) {
          cerr << "Qualche problema serio!\n";          
        } else {
          part[i] = centroids[i]/(float)counts[i];      
        }
      }
      seeds.push_back(part[i]);      
    }

    part.clear();
    for(unsigned int i = 0; i < seeds.size(); i++)
      part.push_back(seeds[i]);

    if(part.size() == 0) part.push_back(Point3f(0,0,0));
    cerr << "Initing!\n";
    part.Init();    
    cerr << "Inited!\n";
    return failed == 0;
}

void VoronoiChain::Init(VFile<Point3f> &baricenters, 
			float scaling, int steps) {
  
  unsigned int f_cells = baricenters.Size() / mean_size;
  unsigned int c_cells = (unsigned int)(scaling * f_cells);

  levels.push_back(VoronoiPartition());
  levels.push_back(VoronoiPartition());
  VoronoiPartition &fine = levels[0];
  VoronoiPartition &coarse = levels[1];

  srand(0);
  
  float coarse_vmean = mean_size/scaling;

  for(unsigned int i = 0; i < baricenters.Size(); i++) {
    int f = (int)(mean_size * (float)rand()/(RAND_MAX + 1.0));
    int c = (int)(coarse_vmean * (float)rand()/(RAND_MAX + 1.0));
    if(f == 2) {
      Point3f &point = baricenters[i];
      fine.push_back(point);
    }
    if(c == 2) {
      Point3f &point = baricenters[i];
      coarse.push_back(point);
    }
  }
  //TODO! Check for duplicates (use the closest :P)
  //  cerr << "fine_seeds.size: " << fine.size() << endl;
  //cerr << "coarse_seeds.size: " << coarse.size() << endl;
  fine.Init();
  coarse.Init();

  //here goes some optimization pass.
  //Fine optimization.
  Report report;
  vector<Point3f> centroids;
  vector<unsigned int> counts;  
  for(int step = 0; step < steps; step++) {
    cerr << "Optimization step: " << step+1 << "/" << steps << endl;
    centroids.clear();
    counts.clear();
    centroids.resize(fine.size(), Point3f(0, 0, 0));
    counts.resize(fine.size(), 0);
    
    report.Init(baricenters.Size());
    for(unsigned int v = 0; v < baricenters.Size(); v++) {
      report.Step(v);
      unsigned int target = fine.Locate(baricenters[v]);
      centroids[target] += baricenters[v];
      counts[target]++;
    }    
    if(step == steps-1) {
      if(!Optimize(mean_size, fine, centroids, counts, false))
	step--;
    } else 
      Optimize(mean_size, fine, centroids, counts, true);
  }
  cerr << "Optimized (fine)!\n";
//here goes some optimization pass.
//Coarse optimization.
//vector<float> radius;
  for(int step = 0; step < steps; step++) {
    cerr << "Optimization step: " << step+1 << "/" << steps << endl;
    centroids.clear();
    counts.clear();
    centroids.resize(coarse.size(), Point3f(0, 0, 0));
    counts.resize(coarse.size(), 0);
    //radius.resize(coarse.size(), 0);
    
    report.Init(baricenters.Size());
    for(unsigned int v = 0; v < baricenters.Size(); v++) {
      if(v & 0xffff) report.Step(v);
      unsigned int ctarget = 0xffffffff;
      ctarget = coarse.Locate(baricenters[v]);
      //      float dist;
      //      coarse.Closest(crude.vert[v], ctarget, dist);
      assert(ctarget != 0xffffffff);
      centroids[ctarget] += baricenters[v];
      counts[ctarget]++;
      //if(dist > radius[ctarget]) radius[ctarget] = dist;
    }
    if(step == steps-1) {
      if(!Optimize((int)coarse_vmean, coarse, centroids, counts, false))
	      step --;
    } else 
      Optimize((int)coarse_vmean, coarse, centroids, counts, true);
  }    
  cerr << "Optimized coarse!\n";
}

unsigned int VoronoiChain::Locate(unsigned int level, 
				  const vcg::Point3f &p) {
  return levels[level].Locate(p);
}

//TODO move this to nxsbuild
void VoronoiChain::RemapFaces(VFile<Point3f> &baricenters,
                              VFile<unsigned int> &face_remap,
			                        vector<unsigned int> &patch_faces,
			                        float scaling, int steps) {
  
  Init(baricenters, scaling, steps);

  //TODO: improve quality of patches and implement threshold.
  typedef  map<pair<unsigned int, unsigned int>, unsigned int> FragIndex;

  //  map<pair<unsigned int, unsigned int>, unsigned int> patches;
  FragIndex patches;

  unsigned int totpatches = 0;

  Point3f bari;
  for(unsigned int i = 0; i < baricenters.Size(); i++) {
    bari = baricenters[i];
    
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

    while(patch_faces.size() <= patch) 
      patch_faces.push_back(0);      
    patch_faces[patch]++;
  }

  //prune faces (now only 0 faces);
  //TODO prune really small faces
  unsigned int tot_patches = 0;
  vector<int> patch_remap;
  for(unsigned int i = 0; i < patch_faces.size(); i++) {
    //if below threshold (and can join faces)
    if(patch_faces[i] == 0)
      patch_remap.push_back(-1);
    else
      patch_remap.push_back(tot_patches++);

    if(patch_faces[i] > 32000) {
      //TODO do something to reduce patch size... :P
      cerr << "Found a patch too big... sorry\n";
      exit(0);
    }
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
  patch_faces.clear();
  patch_faces.resize(totpatches, 0);
  for(unsigned int i = 0; i < face_remap.Size(); i++) {
    unsigned int patch = face_remap[i];
#ifdef CONTROLS
    if(patch == 0xffffffff) {
      cerr << "RESIGH\n";
      exit(0);
    }
    if(patch_remap[patch] == -1) {
        cerr << "SIGH!\n";
        exit(0);
    }
#endif
    unsigned int newpatch = patch_remap[patch];
    face_remap[i] = newpatch;
    patch_faces[newpatch]++;
  }

  
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
  vector<Point3f> centroids;
  vector<unsigned int> counts;

  for(int step = 0; step < steps; step++) {
    cerr << "Optimization step: " << step+1 << "/" << steps << endl;
    centroids.clear();
    counts.clear();
    centroids.resize(coarse.size(), Point3f(0, 0, 0));
    counts.resize(coarse.size(), 0);

    report.Init(nexus.index.size());
    for(unsigned int idx = offset; idx < nexus.index.size(); idx++) {
      report.Step(idx);
      Patch patch = nexus.GetPatch(idx);
      for(unsigned int i = 0; i < patch.nv; i++) {
	
	unsigned int ctarget = coarse.Locate(patch.Vert(i));
	assert(ctarget < coarse.size());
	centroids[ctarget] += patch.Vert(i);
	counts[ctarget]++;
      }
    }
    if(step == steps-1) {
      if(!Optimize((int)coarse_vmean, coarse, centroids, counts, false))
	step--;
    } else 
      Optimize((int)coarse_vmean, coarse, centroids, counts, true);
  }    
  
  newfragments.clear();
}
