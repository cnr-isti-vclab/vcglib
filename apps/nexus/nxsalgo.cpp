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
Revision 1.24  2005/02/22 14:20:44  ponchio
debug and mostly vertex unifying across borders
(still not perfect... :P)

Revision 1.23  2005/02/22 10:38:10  ponchio
Debug, cleaning and optimization.

Revision 1.22  2005/02/21 17:55:36  ponchio
debug debug debug

Revision 1.21  2005/02/20 18:07:01  ponchio
cleaning.

Revision 1.20  2005/02/19 10:45:04  ponchio
Patch generalized and small fixes.

Revision 1.19  2005/02/18 13:04:12  ponchio
Added patch reordering.

Revision 1.18  2005/02/17 16:40:35  ponchio
Optimized BuildLevels.

Revision 1.17  2005/02/08 12:43:03  ponchio
Added copyright


****************************************************************************/

#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <iostream>

//#include <wrap/strip/tristrip.h>

#include "nxsalgo.h"
#include "vpartition.h"
#include "vfile.h"
#include "nexus.h"
#include "zcurve.h"
#include "watch.h"

#include <vcg/space/line3.h>

using namespace std;
using namespace nxs;
using namespace vcg;


void nxs::Connect(Nexus &nexus, std::vector< set<unsigned int> > &close, 
		  float threshold) {

  VPartition grid;
  float max_radius = 0;
  
  for(unsigned int patch = 0; patch < nexus.size(); patch++) {
    Sphere3f &sphere = nexus[patch].sphere;
    grid.push_back(sphere.Center());
    float r = sphere.Radius();
    if(r > max_radius) max_radius = r;
  }
  grid.Init();
  close.clear();
  close.resize(nexus.size());

  vector<int> targets;
  vector<double> dists;
  for(unsigned int patch = 0; patch < nexus.size(); patch++) {
    float radius = nexus[patch].sphere.Radius();
    float max_distance = radius + max_radius + threshold;
    max_distance *= max_distance;
    grid.Closest(grid[patch], targets, dists, max_distance);

    for(unsigned int i = 0; i < targets.size(); i++) {
      unsigned int target = targets[i];
      if(target == patch) continue;
      float dist = radius + nexus[target].sphere.Radius() + threshold;
      dist *= dist;
      if(dist >= dists[i]) {
	close[patch].insert(target);     
      }
    }
  }
  
  //DOUBLECROSS CHECK
  for(unsigned int patch = 0; patch < nexus.size(); patch++) {
    set<unsigned int>::iterator i;
    for(i = close[patch].begin(); i != close[patch].end(); i++) {
      if(!close[*i].count(patch)) {
	cerr << "Some problem width sphere intersection. Have alook.\n";
	cerr << "Meanwhile i fix it.\n";
	close[*i].insert(patch);
      }
    }
  }
}

void nxs::ComputeNormals(Nexus &nexus) { 
  assert(nexus.signature.vnorm);
  assert(!nexus.borders.IsReadOnly()); 

  //first calculate level 0 normals
  //I load all patches in a fragment
  //calculate normals
  //fix external borders getting from level below

  //first try naive approach just load neighborough get normals and
  //fix border with lower level

  vector<int> levels;
  nexus.history.BuildLevels(levels); 
  Report report(nexus.size(), 15);

  //TODO check level 0 is the finer onr

  int current_level = 0;
  while(1) {
    int count = 0;
    for(unsigned int p = 0; p < nexus.size(); p++) {
      if(levels[p] != current_level) continue;
      count++;
      report.Step(p);

      Border &border = nexus.GetBorder(p);

      map<unsigned int, vector<Point3f> > normals;
      normals[p] = vector<Point3f>();

      for(unsigned int i = 0; i < border.Size(); i++) {
	Link &link = border[i];
	if(levels[link.end_patch] == current_level)
	  normals[link.end_patch] = vector<Point3f>();
      }
      map<unsigned int, vector<Point3f> >::iterator k;
      for(k = normals.begin(); k != normals.end(); k++) {
	Patch &patch = nexus.GetPatch((*k).first);

	vector<Point3f> &normal = (*k).second;
	normal.resize(patch.nv, Point3f(0, 0, 0));    
      
	if(nexus.signature.face == Signature::TRIANGLES) {
	  for(unsigned int i = 0; i < patch.nf; i++) {
	    unsigned short *f = patch.Face(i);
	    Point3f &v0 = patch.Vert3f(f[0]);
	    Point3f &v1 = patch.Vert3f(f[1]);
	    Point3f &v2 = patch.Vert3f(f[2]);
	    Point3f norm = (v1 - v0) ^ (v2 - v0);
	    norm.Normalize();
	    normal[f[0]] += norm;
	    normal[f[1]] += norm;
	    normal[f[2]] += norm;
	  }
	} else if(nexus.signature.face == Signature::STRIPS) {
	  for(int i = 0; i < patch.nf - 2; i++) {
	    unsigned short *f = patch.FaceBegin() + i;
	    Point3f &v0 = patch.Vert3f(f[0]);
	    Point3f &v1 = patch.Vert3f(f[1]);
	    Point3f &v2 = patch.Vert3f(f[2]);
	    Point3f norm = (v1 - v0) ^ (v2 - v0); 
	    norm.Normalize();
	    if(i%2) norm = -norm;
	    normal[f[0]] += norm;
	    normal[f[1]] += norm;
	    normal[f[2]] += norm;
	  }
	} else
	  assert(0);
      }

      //now fix borders
      map<unsigned int, vector<Link> > lowers;
      for(unsigned int i = 0; i < border.Size(); i++) {
	Link &link = border[i];
	if(levels[link.end_patch] == current_level) {
	  //TODO remove these asserts
	  assert(normals[p].size() > link.start_vert);
	  assert(normals.count(link.end_patch));
	  assert(normals[link.end_patch].size() > link.end_vert);
	  normals[p][link.start_vert] += normals[link.end_patch][link.end_vert];
	} else if (levels[link.end_patch] < current_level) {
	  lowers[link.end_patch].push_back(link);
	}
      }

      map<unsigned int, vector<Link> >::iterator s;
      for(s = lowers.begin(); s != lowers.end(); s++) {
	Patch &patch = nexus.GetPatch((*s).first);
	for(unsigned int i = 0; i < (*s).second.size(); i++) {
	  Link &link = (*s).second[i];
	  if(nexus.signature.vnorm == Encodings::FLOAT3) 
	    normals[p][link.start_vert] = 
	      ((Point3f *)patch.VNormBegin())[link.end_vert];
	  else if(nexus.signature.vnorm == Encodings::SHORT4) {
	    Point3f &nor = normals[p][link.start_vert];
	    short *n = ((short *)patch.VNormBegin()) + 4*link.end_vert;
	    nor[0] = n[0];
	    nor[1] = n[1];
	    nor[2] = n[2];
	  }
	}
      }
      //copy and normalize
      Patch &patch = nexus.GetPatch(p);
      Entry &entry = nexus[p];
      Point3f *norm = (Point3f *)patch.VNormBegin();
      vector<Point3f> &newnormals = normals[p];
      assert(newnormals.size() == patch.nv);
      for(unsigned int i = 0; i < patch.nv; i++) {
	newnormals[i].Normalize();
	if(nexus.signature.vnorm == Encodings::SHORT4) {
	  newnormals[i] *= 32766;
	  short *np = ((short *)norm) + 4 * i;
	  np[0] = (short)newnormals[i][0];
	  np[1] = (short)newnormals[i][1];
	  np[2] = (short)newnormals[i][2];
	  np[3] = 0;
	} else if(nexus.signature.vnorm == Encodings::FLOAT3) 
	  norm[i] = newnormals[i];
      }
    }
    if(count == 0) break;
    current_level++;
  }
}



/*void nxs::ComputeNormals(Nexus &nexus) {
  assert(nexus.signature.vnorm);

  //setting borders readonly:

  assert(!nexus.borders.IsReadOnly());
  nexus.borders.SetReadOnly(true);
  
  //TODO use a temporary file to store border normals
  unsigned int tmpb_offset = 0;
  vector<unsigned int> tmpb_start;
  VFile<Point3f> tmpb;
  if(!tmpb.Create("tmpb.tmp")) {
    cerr << "Could not create temporary border file\n";
    exit(0);
  }

  for(unsigned int p = 0; p < nexus.size(); p++) {
    Border &border = nexus.borders[p];
    tmpb_start.push_back(tmpb_offset);
    tmpb_offset += border.Size();
  }

  Point3f zero(0.0f, 0.0f, 0.0f);

  tmpb.Resize(tmpb_offset);
  for(unsigned int i = 0; i < tmpb.Size(); i++)
    tmpb[i] = zero;  
  tmpb.Flush();
  

  vector<int> levels;
  nexus.history.BuildLevels(levels);

  //first step normals in the same patch.
  cerr << "First Step\n";
  Report report(nexus.size(), 5);
  vector<Point3f> normals;

  for(unsigned int p = 0; p < nexus.size(); p++) {
    int current_level = levels[p];
    report.Step(p);
    Patch &patch = nexus.GetPatch(p);
    
    normals.clear();  
    normals.resize(patch.nv, Point3f(0, 0, 0));    

    if(nexus.signature.face == Signature::TRIANGLES)
      for(unsigned int i = 0; i < patch.nf; i++) {
	      unsigned short *f = patch.Face(i);
	      Point3f &v0 = patch.Vert3f(f[0]);
	      Point3f &v1 = patch.Vert3f(f[1]);
	      Point3f &v2 = patch.Vert3f(f[2]);
	
	      Point3f norm = (v1 - v0) ^ (v2 - v0); 
	      normals[f[0]] += norm;
	      normals[f[1]] += norm;
	      normals[f[2]] += norm;
      }
    if(nexus.signature.face == Signature::STRIPS)
      for(int i = 0; i < patch.nf - 2; i++) {
	      unsigned short *f = patch.FaceBegin() + i;
	      Point3f &v0 = patch.Vert3f(f[0]);
	      Point3f &v1 = patch.Vert3f(f[1]);
	      Point3f &v2 = patch.Vert3f(f[2]);
	
	      Point3f norm = (v1 - v0) ^ (v2 - v0); 
	      if(i%2) norm = -norm;
	      normals[f[0]] += norm;
	      normals[f[1]] += norm;
	      normals[f[2]] += norm;
      }
    
    if(nexus.signature.vnorm == Encodings::SHORT4) {
      short *n = (short *)patch.VNormBegin();
      for(unsigned int i = 0; i < patch.nv; i++, n += 4) {
	Point3f &norm = normals[i];
	norm.Normalize();
	for(int k = 0; k < 3; k++) 
	  n[k] = (short)(norm[k] * 32766);
	n[3] = 0;
      }
    } else if(nexus.signature.vnorm == Encodings::FLOAT3) {
      Point3f *n = (Point3f *)patch.VNormBegin();
      for(unsigned int i = 0; i < patch.nv; i++) {
	n[i] = normals[i];
	n[i].Normalize();
      }
    }


    Border &border = nexus.GetBorder(p);

    map<unsigned int, map<unsigned short, Point3f> > bnorm;
    map<unsigned int, Link> bcopy;

    unsigned int poff = tmpb_start[p];
    for(unsigned int i = 0; i < border.Size(); i++) {
      Link &link = border[i];
      if(link.IsNull()) continue;  //this should never happen now.
      Point3f pt = normals[link.start_vert];
      if(levels[link.end_patch] == current_level) {
	bnorm[link.end_patch][link.end_vert] = pt;
	tmpb[poff + i] += pt;                          
      } else if(levels[link.end_patch] > current_level) {
       	bcopy[i] = link;
      }
    }

    map<unsigned int, map<unsigned short, Point3f> >::iterator k;
    for(k = bnorm.begin(); k != bnorm.end(); k++) {
      unsigned int patch = (*k).first;
      Border &rborder = nexus.GetBorder(patch);
      unsigned int offset = tmpb_start[patch];
      for(unsigned int i = 0; i < rborder.Size(); i++) {
	      Link &link = rborder[i];
	      //assert(!link.IsNull());
	      //TODO not accurate
        if(link.end_patch != p) continue;
	      if((*k).second.count(link.start_vert))
	        tmpb[offset + i] += (*k).second[link.start_vert];
      }
    }

    //Uncomment this only when links are ok!
    map<unsigned int, Link>::iterator j;
    for(j = bcopy.begin(); j != bcopy.end(); j++) {
      unsigned int b = (*j).first;
      Link link = (*j).second;
      Border &rborder = nexus.GetBorder(link.end_patch, false);
      unsigned int offset = tmpb_start[link.end_patch];
      for(unsigned int i = 0; i < rborder.Size(); i++) {
	Link &rlink = rborder[i];
	if(rlink.end_patch == p && rlink.start_vert == link.end_vert) {
	  assert(rlink.end_vert == link.start_vert);
	  tmpb[poff + b] = tmpb[offset + i];
	}
      }
    }
    } 

  //Second step unify normals across borders
  cerr << "Second step\n";
  report.Init(nexus.size());
  for(unsigned int p = 0; p < nexus.size(); p++) {
    report.Step(p);
    Patch &patch = nexus.GetPatch(p);
    Border &border = nexus.GetBorder(p);

    Point3f *normf = (Point3f *)patch.VNormBegin();
    short *norms = (short *)patch.VNormBegin();

    for(unsigned int i = 0; i < border.Size(); i++) {
      Link &link = border[i];
      if(link.IsNull()) continue;
      unsigned int off = tmpb_start[p];
      //      Point3f &n = tmpb[off + i];
      Point3f n = tmpb[off + i];
      if(n == Point3f(0.0f,0.0f,0.0f)) continue;
      n.Normalize();
      if(nexus.signature.vnorm == Encodings::SHORT4) {
	n *= 32766;
	short *np = norms + 4 * link.start_vert;
	np[0] = (short)n[0];
	np[1] = (short)n[1];
	np[2] = (short)n[2];
	np[3] = 0;
      } else if(nexus.signature.vnorm == Encodings::FLOAT3) {
	normf[link.start_vert] = n;
      }
    }
  }
  tmpb.Close();
  tmpb.Delete();
  //TODO remove temporary file.
  nexus.borders.SetReadOnly(false);
  }*/



/*
  //TODO why i created this function? wonder...void nxs::Reorder(Signature &signature, Patch &patch) {
  vector<unsigned> remap;
  remap.resize(patch.nv, 0xffff);
  
  int nf = patch.nf;
  if(signature.face == Signature::TRIANGLES)
    nf *= 3;
  else if(signature.face != Signature::STRIPS) {
    assert(0); //mah...
  }
  
  //building remap
  unsigned short *f = patch.FaceBegin();
  unsigned int count = 0;
  for(int i = 0; i < nf; i++) {
    assert(f[i] < remap.size());
    if(remap[f[i]] == 0xffff) {
      remap[f[i]] = count++;
    }
  }
  //test no unreferenced vertices
  for(int i = 0; i < patch.nv; i++)
    if(remap[i] == 0xffff)
      remap[i] = i;
  
  //converting faces
  for(int i = 0; i < nf; i++)
    f[i] = remap[f[i]];
  
  vector<Point3f> vert;
  vert.resize(patch.nv);
  memcpy(&*vert.begin(), patch.Vert3fBegin(), patch.nv * sizeof(Point3f));
  for(int i = 0; i < patch.nv; i++)
    patch.Vert3f(remap[i]) = vert[i];
    }*/

//TODO actually use threshold

void nxs::Unify(vector<Point3f> &points, vector<unsigned short> &faces,
		vector<unsigned int> &remap, float threshold) {
  vector<unsigned short> newfaces = faces;

  VPartition grid;
  for(unsigned int i = 0; i < points.size(); i++)
    grid.push_back(points[i]);
  grid.Init();

  remap.resize(points.size());
  vector<int> targets;
  vector<double> dists;

  points.clear();
  unsigned int count = 0;
  for(unsigned int i = 0; i < grid.size(); i++) {

    grid.Closest(grid[i], targets, dists, threshold);
    
    if(targets.size() > 1) {
      unsigned int p;
      for(p = 0; p < targets.size(); p++) {
	if(targets[p] < i) {
	  remap[i] = remap[targets[p]];
	  break;
	}
      }
      if(p < targets.size()) continue;
    }
    remap[i] = count++;
    points.push_back(grid[i]);
  }
  
  //fixing faces now
  faces.clear();
  for(unsigned int i = 0; i < newfaces.size(); i += 3) {
    unsigned short f[3];
    f[0] = remap[newfaces[i]];
    f[1] = remap[newfaces[i+1]];
    f[2] = remap[newfaces[i+2]];
    if(f[0] == f[1] || f[0] == f[2] || f[1] == f[2])
      continue;

    for(int k = 0; k < 3; k++)
      faces.push_back(f[k]);
  }
}

/*void nxs::Unify(Nexus &nexus, float threshold) {
    threshold = 0.00001;
  unsigned int duplicated = 0;
  unsigned int degenerate = 0;

  for(unsigned int p = 0; p < nexus.size(); p++) {
    if(levels[p] != current_level) continue;
    count++;
    report.Step(p);
    
    Border &border = nexus.GetBorder(p);
    
      map<unsigned int, vector<Point3f> > normals;
      normals[p] = vector<Point3f>();
      
      for(unsigned int i = 0; i < border.Size(); i++) {
	Link &link = border[i];
	if(levels[link.end_patch] == current_level)
	  normals[link.end_patch] = vector<Point3f>();
      }
      map<unsigned int, vector<Point3f> >::iterator k;
      for(k = normals.begin(); k != normals.end(); k++) {
	Patch &patch = nexus.GetPatch((*k).first);
      }
      }
      }*/

/*
void nxs::Unify(Nexus &nexus, float threshold) {
  threshold = 0.001;
  //TODO what if colors or normals or strips?
  unsigned int duplicated = 0;
  unsigned int degenerate = 0;

  for(unsigned int p = 0; p < nexus.size(); p++) {
    Entry &entry = nexus[p];
    Patch &patch = nexus.GetPatch(p);


     
    VPartition part;
    for(unsigned int i = 0; i < patch.nv; i++) {
      Point3f &point = patch.Vert3f(i);
      part.push_back(point);
    }
    part.Init();

    unsigned int vcount = 0;
    vector<unsigned short> remap;
    remap.resize(patch.nv);

    int targets[8];
    double dists[8];

    //TODO CRITICAL FIX this unifying routine.
    for(unsigned int i = 0; i < patch.nv; i++) {
      Point3f &point = patch.Vert3f(i);
      part.Closest(point, 8, targets, dists);
      int k = 0;
      for(k = 0; k < 8; k++) {
	if(dists[k] > threshold) {
	  remap[i] = vcount++;
	  break;
	}
	if(targets[k] < i) {
	  remap[i] = remap[targets[k]];
	  duplicated++;
	  break;
	}
      } 
      if(k == 8)
	remap[i] = vcount++;

    }

    if(vcount == patch.nv) //no need to unify
      continue;

    vector<Point3f> newvert;
    newvert.resize(vcount);
    for(unsigned int i = 0; i < patch.nv; i++)
      newvert[remap[i]] = patch.Vert3f(i);

    vector<unsigned short> newface;
    //check no degenerate faces get created.
    for(unsigned int f = 0; f < entry.nface; f++) {
      unsigned short *face = patch.Face(f);
      if(face[0] != face[1] && face[1] != face[2] && face[0] != face[2] &&
	 newvert[remap[face[0]]] != newvert[remap[face[1]]] &&
	 newvert[remap[face[0]]] != newvert[remap[face[2]]] &&
	 newvert[remap[face[1]]] != newvert[remap[face[2]]]) {
	newface.push_back(remap[face[0]]);
	newface.push_back(remap[face[1]]);
	newface.push_back(remap[face[2]]);
      } else {
	degenerate++;
      }
    }

    //rewrite patch now.
    entry.nvert = newvert.size();
    entry.nface = newface.size()/3;
    patch.Init(nexus.signature, entry.nvert, entry.nface);

    memcpy(patch.Vert3fBegin(), &(newvert[0]), entry.nvert*sizeof(Point3f));
    memcpy(patch.FaceBegin(), &(newface[0]), entry.nface*3*sizeof(short));

    //testiamo il tutto...  TODO remove this of course
#ifdef NDEBUG
    for(unsigned int i =0; i < patch.nf; i++) {
      for(int k =0 ; k < 3; k++)
        if(patch.Face(i)[k] >= patch.nv) {
          cerr <<" Unify has problems\n";
          exit(0);
        }
    }
#endif
    //TODO CRITICAL FIX unify vertices across borders..... HOW??????
    
    //fix patch borders now
    set<unsigned int> close; //bordering pathes
    Border &border = nexus.GetBorder(p);
    for(unsigned int b = 0; b < border.Size(); b++) {
      if(border[b].IsNull()) continue;
      close.insert(border[b].end_patch);
      border[b].start_vert = remap[border[b].start_vert];
    }

    set<unsigned int>::iterator c;
    for(c = close.begin(); c != close.end(); c++) {
      Border &bord = nexus.GetBorder(*c);
      for(unsigned int b = 0; b < bord.Size(); b++) {
        if(bord[b].IsNull()) continue;
        if(bord[b].end_patch == p) {
          bord[b].end_vert = remap[bord[b].end_vert];
        }
      }
    }
  }
  //better to compact directly borders than setting them null.

  //finally: there may be duplicated borders
  for(unsigned int p = 0; p < nexus.size(); p++) {
    Border &border = nexus.GetBorder(p);
    set<Link> links;
    for(unsigned int b = 0; b < border.Size(); b++) {
      Link &link = border[b];
      assert(!link.IsNull());
      //if(border[b].IsNull()) continue;
      links.insert(link);
    }
    int count = 0;
    for(set<Link>::iterator k = links.begin(); k != links.end(); k++)
      border[count++] = *k;      
    
    nexus.borders[p].used = links.size();
  }
  
  nexus.totvert -= duplicated;
  if(duplicated)
    cerr << "Found " << duplicated << " duplicated vertices" << endl;
  if(degenerate)
    cerr << "Found " << degenerate << " degenerate face while unmifying\n";
}
*/
void nxs::ZSort(Nexus &nexus, vector<unsigned int> &forward,
		vector<unsigned int> &backward) {
  //lets get a bounding box from the sphere:
  ZCurve zcurve;
  float r = nexus.sphere.Radius();
  Point3f radius(r, r, r);
  zcurve.Set(nexus.sphere.Center());
  zcurve.Add(nexus.sphere.Center() - radius);
  zcurve.Add(nexus.sphere.Center() + radius);

  vector<int> levels;
  nexus.history.BuildLevels(levels);

  forward.clear();

  vector< vector<ZEntry> > entries;

  for(unsigned int i = 0; i < nexus.size(); i++) {
    int level = levels[i];
    while(level >= entries.size()) entries.push_back(vector<ZEntry>());

    ZEntry e;
    e.id = i;
    e.pos = zcurve.Pos(nexus[i].sphere.Center());
    entries[level].push_back(e);
  }

  for(unsigned int i = 0; i < entries.size(); i++) {
    vector<ZEntry> &lev = entries[i];
    std::sort(lev.begin(), lev.end());
    for(unsigned int k = 0; k < lev.size(); k++) 
      forward.push_back(lev[k].id);
  }

  backward.resize(forward.size());
  for(unsigned int i = 0; i < backward.size(); i++) 
    backward[forward[i]] = i;
}

bool nxs::LineIntersect(Nexus &nexus, Extraction &extraction,
			Line3f line, Point3f &hit) {
  //seguiamo la history;
  Point3f tmp;
  if(!Intersection(nexus.sphere, line, hit, tmp))
    return false;

  bool found = false;
  bool min_dist = -1;
  float bar1, bar2, dist;
  for(unsigned int i = 0; i < extraction.draw_size; i++) {
    unsigned int p = extraction.selected[i];
    if(!Intersection(nexus[p].sphere, line, hit, tmp))
      continue;
    Patch &patch = nexus.GetPatch(p);
    if(nexus.signature.face == Signature::TRIANGLES) {
      for(unsigned int i = 0; i < patch.nf; i++) {
	unsigned short *f = patch.Face(i);
	Point3f &v0 = patch.Vert3f(f[0]);
	Point3f &v1 = patch.Vert3f(f[1]);
	Point3f &v2 = patch.Vert3f(f[2]);
	if(Intersection(line, v0, v1, v2, bar1, bar2, dist) &&
	   dist > 0 && 
	   (min_dist == -1 || min_dist > dist)) {
	  hit = v0*(1-bar1-bar2)+v1*bar1+b2*bar2;
	  min_dist = dist;
	  found = true;
	}
      }
    } else if(nexus.signature.face == Signature::STRIPS) {
      for(int i = 0; i < patch.nf - 2; i++) {
	unsigned short *f = patch.FaceBegin() + i;
	Point3f &v0 = patch.Vert3f(f[0]);
	Point3f &v1 = patch.Vert3f(f[1]);
	Point3f &v2 = patch.Vert3f(f[2]);
	if(Intersection(line, v0, v1, v2, bar1, bar2, dist) &&
	   dist > 0 && 
	   (min_dist == -1 || min_dist > dist)) {
	  hit = v0*(1-bar1-bar2)+v1*bar1+b2*bar2;
	  min_dist = dist;
	  found = true;
	}
      }
    } else
      assert(0);
  }
  return found;
}
