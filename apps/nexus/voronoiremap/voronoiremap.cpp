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
Revision 1.1  2004/07/01 21:32:18  ponchio
Created


****************************************************************************/

#include <iostream>
#include "../pchain.h"
#include "../pvoronoi.h"
#include "../vert_remap.h"
#include "../crude.h"
#include "../stopwatch.h"
#include "../pintersect.h"

using namespace std;
using namespace nxs;
using namespace vcg;

StopWatch watch;

int main(int argc, char *argv[]) {
  if(argc < 3) {
    cerr << "Usage: " << argv[0] 
	 << " <crude input> <output (remap and pchain)> "
            "[patch_size] [threshold]\n\n";
    cerr << "Patch size and treshold expressed int faces "
            "and default to 1000 and 250\n";
    return -1;
  }



  Crude crude;
  if(!crude.Load(argv[1])) {
    cerr << "Could not open crude input\n";
    return -1;
  }

  string output = argv[2];

  unsigned int target_size = 1000;
  if(argc > 3)
    target_size = atoi(argv[3]);
  if(target_size <= 0) {
    cerr << "Invalid patch_size: " << argv[3] << endl;
    return -1;
  }

  unsigned int treshold = 250;
  if(argc > 4)
    treshold = atoi(argv[4]);
  if(treshold <= 0) {
    cerr << "Invalid patch_size: " << argv[4] << endl;
    return -1;
  }

  cerr << "Getting optimal radius...\n";
  watch.Start();
  


  vector<unsigned int> target;
  unsigned int patch_size = target_size;
  for(patch_size = target_size; patch_size < crude.vert.Size(); patch_size *=2)
    target.push_back(patch_size);
  vector<float> radius;
  radius = VoronoiPartition::OptimalRadii(crude.vert.Size(),
					  crude.vert.Begin(), 
					  crude.vert.End(),
					  crude.GetBox(),
					  target);
  for(unsigned int i = 0; i < radius.size(); i++) 
  cerr << "Radius: " << radius[i] << endl;

  watch.Stop();
  cerr << " ...done in " << watch.Elapsed() << " secs\n";
  watch.Reset();
  watch.Start();
  
  assert(radius.size() > 1);
  //TODO cosa succede se c'e' una sola patch?

  PChain<VoronoiPartition> chain;

  for(unsigned int i = 0; i < radius.size(); i++) {
    VoronoiPartition part;
    part.Init(crude.GetBox());
    chain.levels.push_back(part);
  }
  
  //TODO move this part to pvoronoi (and create Crude::vert_iterator
  VFile<Point3f>::iterator iter;
  for(iter = crude.vert.Begin(); iter != crude.vert.End(); ++iter) { 
    Point3f &v = *iter;
    VoronoiPartition::Key target;
    float dist;    
    for(unsigned int i = 0; i < radius.size(); i++) {
      VoronoiPartition &part = chain.levels[i];
      dist = part.Closest(v, target);
      if(dist >= radius[i] || dist == -1)
	part.Add(v, radius[i]);
    }
  }
  

  VFile<unsigned int> face_remap;
  if(!face_remap.Create(output + ".rmf")) {
    cerr << "Could not create remap files: " << output << ".frm\n";
    return -1;
  }
  face_remap.Resize(crude.face.Size());


  PIntersect<VoronoiPartition> inter(chain.levels[0], chain.levels[1]);

  cerr << "Splitting faces... ";
  
  Point3f bari;
  for(unsigned int i = 0; i < crude.face.Size(); i++) {
    bari = crude.GetBari(i);
    unsigned int patch = inter.Locate(bari);
    if(patch > 1000)
      cerr << "Patch: " << patch << endl;
    face_remap[i] = patch;
  }

  watch.Stop();
  cerr << "done in " << watch.Elapsed() << " secs\n";
  watch.Reset();
  watch.Start();

  //TODO Prune inter to threshold and relocate faces

  VertRemap vert_remap;
  if(!vert_remap.Create(output)) {
    cerr << "Could not create remap files: " << output << ".rmv and .rmb\n";
    return -1;
  }
  vert_remap.Resize(crude.vert.Size());

  cerr << "Splitting vertices... ";

  for(unsigned int i = 0; i < crude.face.Size(); i++) {
    Crude::Face &face = crude.GetFace(i);
    unsigned int patch = face_remap[i];
    if((i % 10000) == 0)
      cerr << "inserting: " << i << endl;
    for(int k = 0; k < 3; k++) {
      vert_remap.Insert(face[k], patch);
    }
  }
  watch.Stop();
  cerr << "done in " << watch.Elapsed() << " secs\n";
  chain.Save(output + ".chn");

  return 0;						 
}
