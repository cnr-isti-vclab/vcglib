#include <iostream>

#include "nexus.h"
#include "pvoronoi.h"
#include "pchain.h"
#include "pintersect.h"


using namespace std;
using namespace nxs;
using namespace vcg;

//prima tiriamo dentro il pchain. partizioni da 0 a n.
//prendiamo la partizione 1 e per ogni cell prendiamo tutte le cells
//che ci appartengono, facciamo il join, le semplifichiamo,
//e le risplittiamo secondo la partition 2.

int main(int argc, char *argv[]) {
  if(argc != 2) {
    cerr << "Usage: " << argv[0] << " <file>\n";
    return 0;
  }
  
  Nexus nexus;
  if(!nexus.Load(argv[1])) {
    cerr << "Could not load nexus.\n";
    return -1;
  }
  
  //write other level to other nexus (TEST)

  Nexus test;
  if(!test.Create("uffa")) {
    cerr << "Could not create testing nexus\n";
    return -1;
  }

  PChain<VoronoiPartition> pchain;
  if(!pchain.Load(argv[1])) {
    cerr << "Could not load partition chain: " << argv[1] << endl;
    return -1;
  }
  
  PIntersect<VoronoiPartition> pbase(&pchain.levels[0], 
				     &pchain.levels[1]);
  if(!pbase.Load(argv[1])) {
    cerr << "Could not load partition intersect\n";
    return -1;
  }
  
  //First thing we insert all of base cells into history
  map<IPair, ICell>::iterator i;
  for(i = pbase.pairs.begin(); i != pbase.pairs.end(); i++) {
    //ICell &cell = (*i).second;

    //    cerr << "Cell: " << cell.index << " -> " << cell.count << endl;
    //insert cell.index into history.
  }

  unsigned int cell_offset = pbase.cells.size();

  for(unsigned int level = 1; level < pchain.levels.size() -1; level++) {
    
    map<unsigned int, vector<unsigned int> > fragments;
    
    map<IPair, ICell>::iterator i;
    for(i = pbase.pairs.begin(); i != pbase.pairs.end(); i++) {
      IPair p = (*i).first;
      ICell c = (*i).second;
      fragments[p.coarse].push_back(c.index);
    }
    
    PIntersect<VoronoiPartition> height;
    height.Init(&pchain.levels[level], &pchain.levels[level+1], cell_offset);
    
    
    vector<Point3f> newvert;
    vector<unsigned int> newface;
    vector<Link> newbord;
    
    //now join fragment
    map<unsigned int, vector<unsigned int> >::iterator k;
    for(k = fragments.begin(); k != fragments.end(); k++) {
      cerr << "Joining: ";
      vector<unsigned int> &fcells = (*k).second;
      for(unsigned int i = 0; i < fcells.size(); i++)
	cerr << " " << fcells[i] << endl;
      cerr << endl;

      nexus.Join((*k).second, newvert, newface, newbord);
      
      //simplify(mesh);
      
      //if != -1 remap global index to cell index
      map<unsigned int, vector<int> > vert_remap;
      map<unsigned int, unsigned int> vert_count;

      //simply collects faces
      map<unsigned int, vector<int> > face_remap;
      map<unsigned int, unsigned int> face_count;

      for(unsigned int f = 0; f < newface.size(); f += 3) {
	Point3f bari = (newvert[newface[f]] + 
			newvert[newface[f+1]] + 
			newvert[newface[f+2]])/3;

	unsigned int cell = height.Locate(bari);

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
      
      //for every new cell
      map<unsigned int, unsigned int >::iterator c;
      for(c = vert_count.begin(); c != vert_count.end(); c++) {
	unsigned int cell = (*c).first;
	cerr << "Processing cell: " << cell << endl;

	vector<unsigned short> faces;
	vector<Point3f> verts;
	vector<Link> bords;
	
	vector<int> &v_remap = vert_remap[cell];
	vector<int> &f_remap = face_remap[cell];

	for(unsigned int i = 0; i < newvert.size(); i++) {
	  if(v_remap[i] != -1)
	    verts.push_back(newvert[v_remap[i]]);
	}

	assert(verts.size() == vert_count[cell]);
	assert(f_remap.size()/3 == face_count[cell]);

	faces.resize(face_count[cell]*3);

	for(unsigned int i = 0; i < f_remap.size(); i++) {
	  if(v_remap[f_remap[i]] == -1) {
	    cerr << "i: " << i << " f_remap[i]: " << f_remap[i] << endl;
	  }
	  assert(v_remap[f_remap[i]] != -1);
	  faces[i] = v_remap[f_remap[i]];
	}
	//process external borders
	for(unsigned int i = 0; i < newbord.size(); i++) {
	  Link link = newbord[i];
	  if(v_remap[link.start_vert] == -1) continue;
	  link.start_vert = v_remap[link.start_vert];
	  bords.push_back(link);
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
	      link.end_patch = (*t).first;
	      link.start_vert = v_remap[i];
	      link.end_vert = vremapclose[i];
	      bords.push_back(link);
	    }
	  }
	}
	unsigned int patch_idx = test.AddPatch(verts.size(),faces.size()/3,0);
	Patch patch = test.GetPatch(patch_idx);
	memcpy(patch.FaceBegin(), &faces[0], 
	       faces.size() * sizeof(unsigned short));
	memcpy(patch.VertBegin(), &verts[0], verts.size() * sizeof(Point3f));
	
	Nexus::Entry &entry = test.index[patch_idx];
	for(int v = 0; v < verts.size(); v++) {
	  entry.sphere.Add(verts[v]);
	  test.sphere.Add(verts[v]);
	}

	//create new nexus patch
	//collect external borders
      }
      


      //fix borders

      
      cell_offset += vert_remap.size();
      //and last move height -> pbase;
      pbase = height;
    }
    //break in the first cycle.
    break;
  }
  return 0;
}
