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
Revision 1.1  2004/07/04 15:24:50  ponchio
Created


****************************************************************************/

#include <iostream>
#include <map>
#include <set>

#include "nexus.h"

using namespace std;
using namespace vcg;
using namespace nxs;

int main(int argc, char *argv[]) {
  if(argc != 2) {
    cerr << "Usage: " << argv[0] << " <nexus>\n";
    return -1;
  }
  
  Nexus nexus;
  if(!nexus.Load(argv[1])) {
    cerr << "Could not load " << argv[1] << endl;
    return -1;
  }

  for(unsigned int p = 0; p < nexus.index.size(); p++) {
    Nexus::Entry &entry = nexus.index[p];
    Patch patch = nexus.GetPatch(p);
    
    unsigned int vcount = 0;
    map<Point3f, unsigned short> vertices; 
    vector<unsigned short> remap;
    remap.resize(patch.VertSize());
    //    map<unsigned short, unsigned short> remap;
    for(unsigned int i = 0; i < patch.VertSize(); i++) {
      Point3f &point = patch.Vert(i);

      if(!vertices.count(point)) {
	vertices[point] = vcount++;
      } else {
	cerr << "Duplicated point\n";
      }

      remap[i] = vertices[point];
    }
    if(vertices.size() == patch.VertSize()) //no need to unify
      continue;

    vector<Point3f> newvert;
    newvert.resize(vertices.size());
    map<Point3f, unsigned short>::iterator k;
    for(k = vertices.begin(); k != vertices.end(); k++) {
      newvert[remap[(*k).second]] = (*k).first;
    }


    vector<unsigned short> newface;
    newface.resize(patch.FaceSize() * 3);
    for(unsigned int f = 0; f < (unsigned int)(patch.FaceSize() *3); f++) {
      newface[f] = remap[f];
    }
    //rewrite patch now.
    patch.Resize(newvert.size(), newface.size());
    entry.nvert = newvert.size();
    entry.nface = newface.size();
    memcpy(patch.VertBegin(), &(newvert[0]), 
	   patch.VertSize() * sizeof(Point3f));
    memcpy(patch.FaceBegin(), &(newface[0]), 
	   patch.FaceSize() * 3 * sizeof(unsigned short));
    
    //fix patch borders now
    set<unsigned int> close; //bordering pathes
    Border border = nexus.GetBorder(p);
    for(unsigned int b = 0; b < border.Size(); b++) {
      if(border[b].IsNull()) continue;
      close.insert(border[b].end_patch);
      border[b].start_vert = remap[border[b].start_vert];
    }
    
    set<unsigned int>::iterator c;
    for(c = close.begin(); c != close.end(); c++) {
      Border bord = nexus.GetBorder(*c);
      for(unsigned int b = 0; b < bord.Size(); b++) {
	if(border[b].IsNull()) continue;
	if(bord[b].end_patch == p) {
	  bord[b].end_vert = remap[border[b].end_vert];
	}
      }
    }
  }

  //finally: there may be duplicated borders
  for(unsigned int p = 0; p < nexus.index.size(); p++) {
    Border border = nexus.GetBorder(p);
    //Nexus::Entry &entry = nexus.index[p];
    
    set<Link> links;
    for(unsigned int b = 0; b < border.Size(); b++) {
      if(border[b].IsNull()) continue;
      if(links.count(border[b])) 
	border[b] = Link();
      else
	links.insert(border[b]);
    }
  }
  return 0;
}
