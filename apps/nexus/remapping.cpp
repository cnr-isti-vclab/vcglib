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
Revision 1.9  2005/02/08 12:43:03  ponchio
Added copyright


****************************************************************************/

#include <iostream>

#include "remapping.h"
#include "watch.h"


using namespace std;
using namespace vcg;
using namespace nxs;

bool BlockIndex::Save(const string &file) {
  FILE *fp = fopen(file.c_str(), "wb+");
  if(!fp) {
    cerr << "Could not save: " << file << endl;
    return false;
  }
  unsigned int nsize = size();
  fwrite(&nsize, sizeof(unsigned int), 1, fp);
  fwrite(&*begin(), sizeof(BlockEntry), nsize, fp);
  fclose(fp);
  return true;
}

bool BlockIndex::Load(const string &file) {
  FILE *fp = fopen(file.c_str(), "rb");
  if(!fp) {
    cerr << "Could not load: " << file << endl;
    return false;
  }
  unsigned int nsize;
  fread(&nsize, sizeof(unsigned int), 1, fp);
  resize(nsize);
  fread(&*begin(), sizeof(BlockEntry), nsize, fp);
  fclose(fp);
  return true;
}

void nxs::Remap(VChain &chain,
		VFile<vcg::Point3f> &points,
		VFile<unsigned int> &remap,
		BlockIndex &index,
		unsigned int target_size,
		unsigned int min_size,
		unsigned int max_size,
		float scaling,
		int steps) {

  VPartition *finepart = new VPartition;  
  chain.push_back(finepart);
  BuildPartition(*finepart, points, target_size, min_size, max_size, steps);

  VPartition *coarsepart = new VPartition;  
  chain.push_back(coarsepart);
  BuildPartition(*coarsepart, points, 
		 (int)(target_size/scaling), min_size, max_size, steps);


  cerr << "Fine size: " << finepart->size() << endl;
  cerr << "Coarse size: " << coarsepart->size() << endl;


//  typedef  map<pair<unsigned int, unsigned int>, unsigned int> FragIndex;
  typedef map<unsigned int, map<unsigned int, unsigned int> > FragIndex;
  FragIndex patches;

  unsigned int totpatches = 0;
  vector<unsigned int> count;

  Point3f bari;
  for(unsigned int i = 0; i < points.Size(); i++) {
    bari = points[i];
    
    unsigned int fine = finepart->Locate(bari);
    unsigned int coarse = coarsepart->Locate(bari);

    unsigned int patch;
    if(!patches.count(coarse) || !patches[coarse].count(fine)) {    
      patch = totpatches;
      patches[coarse][fine] = totpatches++;
    } else
      patch = patches[coarse][fine];

    remap[i] = patch;

    while(count.size() <= patch) 
      count.push_back(0);      
    count[patch]++;
  }

  for(unsigned int i = 0; i < totpatches; i++) {
    if(count[i] > 32000) {
      //TODO do something to reduce patch size... :P
      cerr << "Found a patch too big... sorry\n";
      exit(0);
    }
  }

  unsigned int mean = 0;
  for(unsigned int i = 0; i < count.size(); i++)
    mean += count[i];
  mean /= count.size();
  
  min_size /= 4;  
  cerr << "Pruning small patches... < " << min_size << " mean: " << mean << endl;


  //prune small patches
  
  vector<int> patch_map;
  patch_map.resize(totpatches);
  for(unsigned int i = 0; i < totpatches; i++)
    patch_map[i] = i;

  for(FragIndex::iterator s = patches.begin(); s != patches.end(); s++) {
    map<unsigned int, unsigned int> &fines = (*s).second;
        
    while(1) {
      if(fines.size() <= 1) break;
      unsigned int inf_fine = 0xffffffff;
      unsigned int inf_count, min_count;
      unsigned int min_fine = 0xffffffff;
      unsigned int min_patch, inf_patch;
      map<unsigned int, unsigned int>::iterator t;
      for(t = fines.begin(); t != fines.end(); t++) {
        unsigned int c = count[(*t).second];        
        if(inf_fine == 0xffffffff || c < inf_count) {
          if(inf_fine != 0xffffffff) {
            min_fine = inf_fine;
            min_count = inf_count;
            min_patch = inf_patch;
          }
            inf_fine = (*t).first;
            inf_count = c;
            inf_patch = (*t).second;
        } else if(min_fine == 0xffffffff || c < min_count) {
            min_fine = (*t).first;
            min_count = c;
            min_patch = (*t).second;
        }
      }
      if(inf_count >= min_size || 
         inf_count + min_count > max_size) break;      
      
      count[min_patch] += count[inf_patch];
      patch_map[inf_patch] = min_patch;
      fines.erase(inf_fine);      
    }   
  }
  for(unsigned int i = 0; i < totpatches; i++) 
    while(patch_map[patch_map[i]] != patch_map[i])
      patch_map[i] = patch_map[patch_map[i]];
  
  //now we remap remaining patches into 0 - n.
  unsigned int new_totpatches = 0;
  vector<int> patch_remap;
  patch_remap.resize(totpatches, -1);
  for(unsigned int i = 0; i < totpatches; i++) {    
    unsigned int p = patch_map[i];
    if(patch_remap[p] == -1) 
      patch_remap[p] = new_totpatches++;
    patch_remap[i] = patch_remap[p];
  }  
  
  cerr << "Building fragments\n";

  //building fragments
  for(FragIndex::iterator s = patches.begin(); s != patches.end(); s++) {
    unsigned int coarse = (*s).first;
    map<unsigned int, unsigned int> &fines = (*s).second;
    map<unsigned int, unsigned int>::iterator t;
    for(t = fines.begin(); t != fines.end(); t++) {
      unsigned int fine = (*t).first;
      unsigned int oldpatch = (*t).second;
      assert(oldpatch < patch_remap.size());    
      unsigned int patch = patch_remap[oldpatch];
      if(patch != -1) //not deleted...
        chain.oldfragments[coarse].insert(patch);
    }  
  }

  cerr << "remapping faces again\n";
  //remapping faces
  index.resize(new_totpatches);
  for(unsigned int i = 0; i < remap.Size(); i++) {
    unsigned int patch = remap[i];
#ifdef CONTROLS
    if(patch == 0xffffffff) {
      cerr << "RESIGH\n";
      exit(0);
    }
    if(patch_remap[patch] == -1) {//must relocate this thing....
      //TODO
      cerr << "Could not do this\n";
      exit(0);
    }
#endif


    unsigned int newpatch = patch_remap[patch];
    assert(newpatch < index.size());
    remap[i] = newpatch;
    BlockEntry &entry = index[newpatch];
    entry.size++;
  }

  cerr << "fixing offsets in index\n";
  //Fixing offset
  int64 offset = 0;
  for(unsigned int i = 0; i < index.size(); i++) {
    assert(index[i].size < 65000);
    index[i].offset = offset;
    offset += index[i].size;
  }
  
}
  

void nxs::BuildPartition(VPartition &part,
			 VFile<vcg::Point3f> &points,
			 unsigned int target_size,
			 unsigned int min_size,
			 unsigned int max_size,
			 int steps) {

  //TODO: improve quality of patches and implement threshold.
  unsigned int ncells = points.Size()/target_size;
  cerr << "Target partition size: " << ncells 
       << " mean: " << points.Size()/ncells << endl;
  srand(0);

  for(unsigned int i = 0; i < points.Size(); i++) {
    int f = (int)(target_size * (float)rand()/(RAND_MAX + 1.0));
    if(f == 2) {
      Point3f &point = points[i];
      part.push_back(point);
    }
  }

  //TODO! Check for duplicates (use the closest :P)
  part.Init();


  vector<Point3f> centroids;
  vector<unsigned int> counts;  
  
  for(int step = 0; step < steps; step++) {
    cerr << "Optimization step: " << step+1 << "/" << steps << endl;

    centroids.clear();
    counts.clear();
    centroids.resize(part.size(), Point3f(0, 0, 0));
    counts.resize(part.size(), 0);

    Report report(points.Size());

    for(unsigned int v = 0; v < points.Size(); v++) {
      report.Step(v);

      unsigned int target = part.Locate(points[v]);
      centroids[target] += points[v];
      counts[target]++;
    }    

    for(unsigned int v = 0; v < centroids.size(); v++)
      if(counts[v] != 0)
	centroids[v]/= counts[v];
    
    double quality = 0;
    for(int i = 0; i < part.size(); i++)
      quality += (counts[i] - target_size) * (counts[i] - target_size);
    
    cerr << "Quality: " << quality << endl;

    if(step == steps-1) {
      if(!Optimize(part, ncells, target_size, min_size, max_size, 
		   centroids, counts, false))
	step--;
    } else 
      Optimize(part, ncells, target_size, min_size, max_size, 
	       centroids, counts, true);
  }
  cerr << "Partition size: " << part.size() 
       << " mean: " << (float)(points.Size()/part.size()) << endl << endl;
}

void nxs::BuildLevel(VChain &chain,
		     Nexus &nexus,
		     unsigned int offset, 
		     float scaling,
		     unsigned int target_size,
		     unsigned int min_size,
		     unsigned int max_size,
		     int steps) { 

  unsigned int totface = 0;
  unsigned int totvert = 0;
  for(unsigned int idx = offset; idx < nexus.size(); idx++) {
    totface += nexus[idx].nface;
    totvert += nexus[idx].nvert;
  }

  VPartition *fine = chain[chain.size()-1];  
  fine->Init();

  VPartition *coarse = new VPartition;
  chain.push_back(coarse);
  
  //unsigned int ncells = (unsigned int)(fine.size() * scaling);
  unsigned int ncells = (unsigned int)(scaling * totface/target_size);
  
  //TODO this method for selecting the seeds is ugly!
  float ratio = ncells/(float)(nexus.size() - offset);
  float cratio = 0;
  for(unsigned int idx = offset; idx < nexus.size(); idx++) {
    cratio += ratio;
    if(cratio > 1) {
      Patch patch = nexus.GetPatch(idx);
      Point3f &v = patch.Vert3f(0);
      coarse->push_back(v);
      cratio -= 1;
    }
  }
  
  if(coarse->size() == 0) {
    Patch patch = nexus.GetPatch(0);
    coarse->push_back(patch.Vert3f(0));
  }
  
  float coarse_vmean = totface/(float)coarse->size();
  
  coarse->Init();
  cerr << "Ncells: " << ncells << endl;
  cerr << "Coarse size: " << coarse->size() << endl;
  cerr << "Coarse mean: " << coarse_vmean << " mean_size: " << target_size << endl;

  //here goes some optimization pass.
  //Coarse optimization.
  vector<Point3f> centroids;
  vector<unsigned int> counts;

  for(int step = 0; step < steps; step++) {
    cerr << "Optimization step: " << step+1 << "/" << steps << endl;
    centroids.clear();
    counts.clear();
    centroids.resize(coarse->size(), Point3f(0, 0, 0));
    counts.resize(coarse->size(), 0);

    Report report(nexus.size());
    for(unsigned int idx = offset; idx < nexus.size(); idx++) {
      report.Step(idx);
      Patch patch = nexus.GetPatch(idx);
      for(unsigned int i = 0; i < patch.nf; i++) {
        unsigned short *face = patch.Face(i);	                                          
        Point3f bari = (patch.Vert3f(face[0]) + 
			patch.Vert3f(face[1]) + 
			patch.Vert3f(face[2]))/3;

	unsigned int target = coarse->Locate(bari);
	assert(target < coarse->size());
	centroids[target] += bari;
	counts[target]++;
      }
    }
    
    for(unsigned int v = 0; v < centroids.size(); v++)
      if(counts[v] != 0)
	centroids[v]/= counts[v];
    
    if(step == steps-1) {
      if(!Optimize(*coarse, ncells, (int)coarse_vmean, min_size, max_size, 
		   centroids, counts, false))
	      step--;
    } else 
      Optimize(*coarse, ncells, (int)coarse_vmean, min_size, max_size, 
	       centroids, counts, true);
  }    
  chain.newfragments.clear();
}

int nxs::GetBest(VPartition &part, unsigned int seed,
			  vector<bool> &mark,
			  vector<unsigned int> &counts) {
  
  vector<int> nears;
  vector<float> dists;
  int nnear = 7;
  if(part.size() < 7) nnear = part.size()/2;
  if(!nnear) return -1;
  
  part.Closest(part[seed], nnear, nears, dists);
  int best = -1;
  int bestcount = -1;
  int bestdist = -1;

  for(int k = 0; k < nnear; k++) {
    int c = nears[k];
    if(c == seed) continue;
    assert(c >= 0);
    assert(c < part.size());    
    if(mark[c]) continue;
    
    if(bestcount < 0 || 
       (counts[c] < bestcount)) {
      best = c;
      bestcount = counts[c];
    }
  }
  return best;
}

bool nxs::Optimize(VPartition &part, 
		   unsigned int target_cells,
		   unsigned int target_size,
		   unsigned int min_size,
		   unsigned int max_size,
		   vector<Point3f> &centroids,
		   vector<unsigned int> &counts, 
		   bool join) {

  if(max_size > target_size *3)
    max_size = target_size * 3;
  min_size = (unsigned int)(target_size * 0.3f);

  unsigned int toobig = 0;
  unsigned int toosmall = 0;
  for(unsigned int i = 0; i < part.size(); i++) {
    if(counts[i] > max_size) toobig++;
    if(counts[i] < min_size) toosmall--;
  }
  
  unsigned int close = part.size()/2;
  if(close < 1) close = 1;
  if(close > 10) close = 10;

  unsigned int failed = 0;
  vector<Point3f> seeds;
  vector<bool> mark;
  mark.resize(part.size(), false);
  
  vector<int> nears;
  vector<float> dists;
  //removing small ones.
  for(unsigned int i = 0; i < part.size(); i++) {
    if(counts[i] > max_size) {
      float radius;
      if(part.size() == 1) 
	radius = 0.00001;
      else
	radius = part.Radius(i)/4;
      seeds.push_back(centroids[i] + Point3f(1, -1, 1) * radius);
      seeds.push_back(centroids[i] + Point3f(-1, 1, 1) * radius);
      seeds.push_back(centroids[i] + Point3f(-1, -1, -1) * radius);
      seeds.push_back(centroids[i] + Point3f(1, 1, -1) * radius);
      continue;
    }
    if(counts[i] < min_size) 
      continue;

    part.Closest(part[i], close, nears, dists);
    Point3f dir(0,0,0);
    
    for(unsigned int k = 0; k < close; k++) {
      unsigned int n = nears[k];
      float c = (target_size - (float)counts[n])/
	((float)target_size * close);

      dir += (centroids[i] - part[n]) * c;
    }
    seeds.push_back(centroids[i] + dir);
  }
  part.clear();
  for(unsigned int i = 0; i < seeds.size(); i++)
    part.push_back(seeds[i]);
  
  if(part.size() == 0) {
    cerr << "OOOPS i accidentally deleted all seeds... backup :P\n";
    part.push_back(Point3f(0,0,0));
  }
  part.Init();    
  return true;
}
