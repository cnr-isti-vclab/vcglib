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

  chain.push_back(VPartition());
  BuildPartition(chain.back(), points, target_size, min_size, max_size, steps);



  chain.push_back(VPartition());
  BuildPartition(chain.back(), points, 
		 (int)(target_size/scaling), min_size, max_size, steps);



  VPartition &finepart = chain[0];
  finepart.Init();
  cerr << "Fine size: " << finepart.size() << endl;

  VPartition &coarsepart = chain[1];
  coarsepart.Init();
  cerr << "Coarse size: " << coarsepart.size() << endl;

  typedef  map<pair<unsigned int, unsigned int>, unsigned int> FragIndex;

  FragIndex patches;

  unsigned int totpatches = 0;
  vector<unsigned int> count;

  Point3f bari;
  for(unsigned int i = 0; i < points.Size(); i++) {
    bari = points[i];
    
    unsigned int fine = finepart.Locate(bari);
    unsigned int coarse = coarsepart.Locate(bari);

    unsigned int patch;
    
    if(!patches.count(make_pair(coarse, fine))) {
      patch = totpatches;
      patches[make_pair(coarse, fine)] = totpatches++;
    } else
      patch = patches[make_pair(coarse, fine)];

    remap[i] = patch;

    while(count.size() <= patch) 
      count.push_back(0);      
    count[patch]++;
  }

  cerr << "Prune faces...\n";

  //prune faces (now only 0 faces);
  unsigned int new_totpatches = 0;
  vector<int> patch_remap;
  for(unsigned int i = 0; i < totpatches; i++) {
    //if below threshold (and can join faces)
    //    if(count[i] < min_size)
    if(count[i] == 0)
      patch_remap.push_back(-1);
    else
      patch_remap.push_back(new_totpatches++);
    
    if(count[i] > 32000) {
      //TODO do something to reduce patch size... :P
      cerr << "Found a patch too big... sorry\n";
      exit(0);
    }
  }
  

  cerr << "BUilding fragmenbts\n";

  //building fragments
  FragIndex::iterator f;
  for(f = patches.begin(); f != patches.end(); f++) {
    unsigned int coarse = (*f).first.first;
    unsigned int fine = (*f).first.second;
    unsigned int oldpatch = (*f).second;
    assert(oldpatch < patch_remap.size());
    unsigned int patch = patch_remap[oldpatch];
    if(patch != -1) //not deleted...
      chain.oldfragments[coarse].insert(patch);
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
    assert(newpatch >= 0);
    remap[i] = newpatch;
    BlockEntry &entry = index[newpatch];
    entry.size++;
  }

  cerr << "fixing offsets in index\n";
  //Fixing offset
  int64 offset = 0;
  for(unsigned int i = 0; i < index.size(); i++) {
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
    
    if(step == steps-1) {
      if(!Optimize(part, target_size, min_size, max_size, 
		   centroids, counts, false))
	step--;
    } else 
      Optimize(part, target_size, min_size, max_size, 
	       centroids, counts, true);
  }
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
  for(unsigned int idx = offset; idx < nexus.index.size(); idx++) {
    totface += nexus.index[idx].nface;
    totvert += nexus.index[idx].nvert;
  }
  
  chain.push_back(VPartition());
  VPartition &coarse = chain[chain.size()-1];
  VPartition &fine = chain[chain.size()-2];
  fine.Init();
  
  unsigned int ncells = (unsigned int)(fine.size() * scaling);
  
  //TODO this method for selecting the seeds is ugly!
  float ratio = ncells/(float)(nexus.index.size() - offset);
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
  
  float coarse_vmean = totface/(float)coarse.size();
  
  coarse.Init();
  cerr << "Ncells: " << ncells << endl;
  cerr << "Coarse size: " << coarse.size() << endl;
  cerr << "Coarse mean: " << coarse_vmean << " mean_size: " << target_size << endl;

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

    Report report(nexus.index.size());
    for(unsigned int idx = offset; idx < nexus.index.size(); idx++) {
      report.Step(idx);
      Patch patch = nexus.GetPatch(idx);
      for(unsigned int i = 0; i < patch.nf; i++) {
        unsigned short *face = patch.Face(i);	                                          
        Point3f bari = (patch.Vert(face[0]) + 
			patch.Vert(face[1]) + 
			patch.Vert(face[2]))/3;

	unsigned int target = coarse.Locate(bari);
	assert(target < coarse.size());
	centroids[target] += bari;
	counts[target]++;
      }
    }
    
    for(unsigned int v = 0; v < centroids.size(); v++)
      if(counts[v] != 0)
	centroids[v]/= counts[v];
    
    if(step == steps-1) {
      if(!Optimize(coarse, (int)coarse_vmean, min_size, max_size, 
		   centroids, counts, false))
	      step--;
    } else 
      Optimize(coarse, (int)coarse_vmean, min_size, max_size, 
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
		   unsigned int target_size,
		   unsigned int min_size,
		   unsigned int max_size,
		   vector<Point3f> &centroids,
		   vector<unsigned int> &counts, 
		   bool join) {
  
    unsigned int failed = 0;
    vector<Point3f> seeds;
    vector<bool> mark;
    mark.resize(part.size(), false);

    //first pass we check only big ones
    for(unsigned int i = 0; i < part.size(); i++) {
      if(counts[i] > max_size || counts[i] > 2 * target_size) {
	failed++;
	float radius;

	if(part.size() == 1) 
	  radius = 0.00001;
	else
	  radius = part.Radius(i);
	
        if(radius == 0) {
          cerr << "Radius zero???\n";
          exit(0);
        }
	radius /= 4;
	seeds.push_back(centroids[i] + Point3f(1, -1, 1) * radius);
	seeds.push_back(centroids[i] + Point3f(-1, 1, 1) * radius);
	seeds.push_back(centroids[i] + Point3f(-1, -1, -1) * radius);
        seeds.push_back(centroids[i] + Point3f(1, 1, -1) * radius);
	mark[i];
      }
    }
    
    if(join) {
      for(unsigned int i = 0; i < part.size(); i++) {
	if(mark[i] || counts[i] >= min_size) continue;

	failed++;
	int best = GetBest(part, i, mark, counts);
	if(best < 0) {
	  cerr << "Best not found! while looking for: " << i << "\n";
	  continue;
	}
	assert(mark[best] == false);
	mark[best] = true;
	mark[i] = true;
	seeds.push_back((part[i] + part[best])/2);
      }
    }

    for(unsigned int i = 0; i < part.size(); i++) {
      if(mark[i]) continue;
      if(counts[i] == 0) continue;
      if(join) {
	//        if(counts[i] < min_size) {
	//          cerr << "Could not fix: " << i << endl;
	//        } else {
	part[i] = centroids[i];
	//        }
      }
      seeds.push_back(part[i]);      
    }
    
    part.clear();
    for(unsigned int i = 0; i < seeds.size(); i++)
      part.push_back(seeds[i]);
    
    if(part.size() == 0) {
      cerr << "OOOPS i accidentally deleted all seeds... backup :P\n";
      part.push_back(Point3f(0,0,0));
    }
    part.Init();    
    return failed == 0;
}
