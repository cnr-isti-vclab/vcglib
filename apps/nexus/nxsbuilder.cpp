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
  if(!face_remap.Create(output + ".rmf")) {
    cerr << "Could not create remap files: " << output << ".rmf\n";
    exit(0);
  }
  face_remap.Resize(crude.Faces());

  VFile<Point3f> baricenters;
  if(!baricenters.Create(output + string(".bvr"))) {
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
  if(!face_index.Save(output + ".rmi")) {
    cerr << "Could not save file: " << output << ".rmi\n";
    exit(0);
  }

  baricenters.Delete();
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

  return 0;
}

