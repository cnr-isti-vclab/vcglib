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
#include "getopt.h"
#else
#include <unistd.h>
#endif

#include <iostream>
using namespace std;

#include "crude.h"
#include "nexus.h"

#include "voronoichain.h"
#include "pintersect.h"
#include "vert_remap.h"

#include "decimate.h"
#include "nxsbuild.h"
using namespace vcg;
using namespace nxs;

/*void RemapVertices(Crude &crude,
		   VertRemap &vert_remap,
		   VFile<unsigned int> &face_remap,	 
		   vector<unsigned int> &patch_verts);

void NexusAllocate(Crude &crude,
		   Nexus &nexus,
		   VFile<unsigned int> &face_remap,
		   vector<unsigned int> &patch_faces,
		   vector<unsigned int> &patch_verts);

void NexusFill(Crude &crude,
	       Nexus &nexus,
	       VertRemap &vert_remap,
	       VFile<RemapLink> &border_remap);

void NexusFixBorder(Nexus &nexus, 
VFile<RemapLink> &border_remap);*/

void NexusSplit(Nexus &nexus, VoronoiChain &vchain,
		unsigned int level,
		vector<Point3f> &newvert, 
		vector<unsigned int> &newface,
		vector<Link> &newbord,
		Nexus::Update &update, 
		float error);

/*float Decimate(unsigned int target_faces, 
	       vector<Point3f> &newvert, 
	       vector<unsigned int> &newface,
	       vector<Link> &newbord,
	       vector<int> &vert_remap);*/

void ReverseHistory(vector<Nexus::Update> &history);



int main(int argc, char *argv[]) {

  Decimation decimation = CLUSTER;
  unsigned int patch_size = 1000;
  unsigned int patch_threshold = 0xffffffff;
  unsigned int optimization_steps = 5;
  bool stop_after_remap = false;
  unsigned int max_level = 0xffffffff;
  float scaling = 0.5;

  int option;
  while((option = getopt(argc, argv, "f:t:l:s:d:ro:")) != EOF) {
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
    case 's': scaling = atof(optarg);
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
  if(!nexus.Create(output, NXS_FACES)) {
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
    cerr << "Could not create remap files: " << output << ".frm\n";
    return -1;
  }
  face_remap.Resize(crude.Faces());


  VertRemap vert_remap;
  if(!vert_remap.Create(output)) {
    cerr << "Could not create remap files: " << output << ".rmv and .rmb\n";
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

  vector<unsigned int> patch_faces;
  vchain.RemapFaces(crude, face_remap, patch_faces, scaling, optimization_steps);
  
  vector<unsigned int> patch_verts;
  patch_verts.resize(patch_faces.size(), 0);
  RemapVertices(crude, vert_remap, face_remap, patch_verts);
  
  if(stop_after_remap) return 0;

  //allocate chunks for patches and copy faces (absoklute indexing) into them.
  NexusAllocate(crude, nexus, face_remap, patch_faces, patch_verts);

  //insert vertices and remap faces, prepare borders
  NexusFill(crude, nexus, vert_remap, border_remap);

  NexusFixBorder(nexus, border_remap);

  //filling history
  Nexus::Update update;
  for(unsigned int i = 0; i < nexus.index.size(); i++) 
    update.created.push_back(i);
  nexus.history.push_back(update);

  //unify vertices otherwise you may get cracks.
  nexus.Unify();

  /* BUILDING OTHER LEVELS */
  unsigned int oldoffset = 0;

  for(unsigned int level = 1; level < max_level; level++) {
    cerr << "Level: " << level << endl;

    unsigned int newoffset = nexus.index.size();
    vchain.BuildLevel(nexus, oldoffset, scaling, optimization_steps);
    
    cerr << "Level built\n";
    vector<Nexus::Update> level_history;
    map<unsigned int, set<unsigned int> >::iterator fragment;
    for(fragment = vchain.oldfragments.begin(); 
	fragment != vchain.oldfragments.end(); fragment++) {


      update.created.clear();
      update.erased.clear();

      cerr << "Join ";
      set<unsigned int> &fcells = (*fragment).second;
      set<unsigned int>::iterator s;
      for(s = fcells.begin(); s != fcells.end(); s++) {
	update.erased.push_back(*s);
	cerr << *s << " ";
      }
      cerr << endl;
      
      vector<Point3f> newvert;
      vector<unsigned int> newface;
      vector<Link> newbord;
      
      nexus.Join((*fragment).second, newvert, newface, newbord);

      //simplyfy mesh
      vector<int> vert_remap;
      float error = Decimate(decimation,
			     (unsigned int)((newface.size()/3) * scaling), 
			     newvert, newface, newbord, vert_remap);


      NexusSplit(nexus, vchain, level, newvert, newface, newbord, 
		 update, error);


      level_history.push_back(update);
    }

    for(int i = 0; i < level_history.size(); i++)
      nexus.history.push_back(level_history[i]);
    //if(vchain.levels.back().size() == 1) break;
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

  //Clean up:
  nexus.Close();

  //TODO remove vert_remap, face_remap, border_remap
  return 0;
}

void NexusSplit(Nexus &nexus, VoronoiChain &vchain,
		unsigned int level,
		vector<Point3f> &newvert, 
		vector<unsigned int> &newface,
		vector<Link> &newbord, 
		Nexus::Update &update,
		float error) {

  //if != -1 remap global index to cell index (first arg)
  map<unsigned int, vector<int> > vert_remap;
  map<unsigned int, unsigned int> vert_count;
  
  //simply collects faces
  map<unsigned int, vector<int> > face_remap;
  map<unsigned int, unsigned int> face_count;
  
  for(unsigned int f = 0; f < newface.size(); f += 3) {
    Point3f bari = (newvert[newface[f]] + 
		    newvert[newface[f+1]] + 
		    newvert[newface[f+2]])/3;
    
    unsigned int cell = vchain.Locate(level+1, bari);
    vector<int> &f_remap = face_remap[cell];
    f_remap.push_back(newface[f]);
    f_remap.push_back(newface[f+1]);
    f_remap.push_back(newface[f+2]);
    face_count[cell]++;
    
    if(!vert_remap.count(cell)) {
      vert_remap[cell].resize(newvert.size(), -1);
      vert_count[cell] = 0;
    }
    
    vector<int> &v_remap = vert_remap[cell];
    
    for(int i = 0; i < 3; i++)
      if(v_remap[newface[f+i]] == -1)
	v_remap[newface[f+i]] = vert_count[cell]++;
  }

  //TODO prune small count cells 
  
  //lets count borders
  map<unsigned int, unsigned int> bord_count;

  map<unsigned int, unsigned int >::iterator c;
  for(c = vert_count.begin(); c != vert_count.end(); c++) {
    unsigned int cell = (*c).first;
    unsigned int &count = bord_count[cell];
    count = 0;

    vector<int> &v_remap = vert_remap[cell];

    //external borders
    for(unsigned int i = 0; i < newbord.size(); i++) {
      Link link = newbord[i];
      if(v_remap[link.start_vert] == -1) continue;
      count++;
    }

    //process internal borders;
    //TODO higly inefficient!!!
    map<unsigned int, unsigned int >::iterator t;  
    for(t = vert_count.begin(); t != vert_count.end(); t++) {
      if(cell == (*t).first) continue;
      vector<int> &vremapclose = vert_remap[(*t).first];
      for(unsigned int i = 0; i < newvert.size(); i++) {
	if(v_remap[i] != -1 && vremapclose[i] != -1) {
	  count++;
	}
      }
    }
  }

  map<unsigned int, unsigned int> cells2patches;

  //lets allocate space
  for(c = vert_count.begin(); c != vert_count.end(); c++) {
    unsigned int cell = (*c).first;
    unsigned int patch_idx = nexus.AddPatch(vert_count[cell],
					    face_count[cell],
					    3 * bord_count[cell]);
    
    //why double border space? because at next level
    //we will need to add those borders... 
    cells2patches[cell] = patch_idx;
    vchain.newfragments[cell].insert(patch_idx);
    update.created.push_back(patch_idx);
  }

  
  //fill it now.
  for(c = vert_count.begin(); c != vert_count.end(); c++) {
    unsigned int cell = (*c).first;
    unsigned int patch_idx = cells2patches[cell];
    
    //vertices first
    vector<int> &v_remap = vert_remap[cell];

    vector<Point3f> verts;
    verts.resize(vert_count[cell]);
    for(unsigned int i = 0; i < newvert.size(); i++) {
      if(v_remap[i] != -1)
	verts[v_remap[i]] = newvert[i];
    }

    //faces now
    vector<int> &f_remap = face_remap[cell];

    vector<unsigned short> faces;
    faces.resize(face_count[cell]*3);

    for(unsigned int i = 0; i < f_remap.size(); i++) {
      assert(v_remap[f_remap[i]] != -1);
      faces[i] = v_remap[f_remap[i]];
    }

    if(patch_idx == 68)
      cerr << "68 bord: " << bord_count[cell] << endl;
    //borders last
    vector<Link> bords;

    //process external borders 
    //for every esternal link we must update external patches!
    for(unsigned int i = 0; i < newbord.size(); i++) {
      Link link = newbord[i];
      if(v_remap[link.start_vert] == -1) continue;
      link.start_vert = v_remap[link.start_vert];
      bords.push_back(link);

      Nexus::Entry &rentry = nexus.index[link.end_patch];
      //TODO if !true reallocate borders.

      Border rborder = nexus.GetBorder(link.end_patch);
      if(rentry.border_used >= rentry.border_size) {
	cerr << "patch: " << link.end_patch << endl;
	cerr << "used: " << rentry.border_used << endl;
	cerr << "size: " << rentry.border_size << endl;
	unsigned int start = nexus.borders.Size();
	nexus.borders.Resize(nexus.borders.Size() + 2 * rentry.border_size);
	Link *tmp = new Link[rentry.border_size];
	memcpy(tmp, &rborder[0], sizeof(Link) * rentry.border_size);
	rentry.border_start = start;
	rentry.border_size *= 2;
	rborder = nexus.GetBorder(link.end_patch);
	memcpy(&rborder[0], tmp, sizeof(Link) * rentry.border_used);
	delete []tmp;
      }
      assert(rentry.border_used < rentry.border_size);


      Link &newlink = rborder[rentry.border_used++];
      newlink.start_vert = link.end_vert;
      newlink.end_vert = link.start_vert;
      newlink.end_patch = patch_idx;
    }

    //process internal borders;
    //TODO higly inefficient!!!
    map<unsigned int, unsigned int >::iterator t;  
    for(t = vert_count.begin(); t != vert_count.end(); t++) {
      if(cell == (*t).first) continue;
      vector<int> &vremapclose = vert_remap[(*t).first];
      for(unsigned int i = 0; i < newvert.size(); i++) {
	if(v_remap[i] != -1 && vremapclose[i] != -1) {
	  Link link;
	  link.end_patch = cells2patches[(*t).first];
	  link.start_vert = v_remap[i];
	  link.end_vert = vremapclose[i];
	  bords.push_back(link);
	}
      }
    }


    Nexus::Entry &entry = nexus.index[patch_idx];
    entry.error = error;
    
    Patch patch = nexus.GetPatch(patch_idx);
    memcpy(patch.FaceBegin(), &faces[0], 
	   faces.size() * sizeof(unsigned short));
    memcpy(patch.VertBegin(), &verts[0], verts.size() * sizeof(Point3f));
    
    
    for(int v = 0; v < verts.size(); v++) {
      entry.sphere.Add(verts[v]);
      nexus.sphere.Add(verts[v]);
    } 
    
    Border border = nexus.GetBorder(patch_idx);
    memcpy(&(border[0]), &(bords[0]), bords.size() * sizeof(Link));
    entry.border_used = bords.size();
  }  
}
 
void ReverseHistory(vector<Nexus::Update> &history) {
  reverse(history.begin(), history.end());
  vector<Nexus::Update>::iterator i;
  for(i = history.begin(); i != history.end(); i++)
    swap((*i).erased, (*i).created);
}
