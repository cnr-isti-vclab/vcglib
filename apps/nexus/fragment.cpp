#include "fragment.h"
#include "border.h"
#include "pvoronoi.h"
#include <iostream>

using namespace std;
using namespace vcg;
using namespace nxs;
using namespace pt;

void NxsPatch::Write(outstm *out) {
  int vsize = vert.size();
  int fsize = face.size();
  int bsize = bord.size();

  out->write(&patch, sizeof(unsigned int));
  out->write(&vsize, sizeof(unsigned int));
  out->write(&fsize, sizeof(unsigned int));
  out->write(&bsize, sizeof(unsigned int));

  out->write(&*vert.begin(), vert.size() * sizeof(Point3f));
  out->write(&*face.begin(), face.size() * sizeof(unsigned short));
  out->write(&*bord.begin(), bord.size() * sizeof(Link));
}

void NxsPatch::Read(instm *in) {
  int vsize;
  int fsize;
  int bsize;

  in->read(&patch, sizeof(unsigned int));
  in->read(&vsize, sizeof(unsigned int));
  in->read(&fsize, sizeof(unsigned int));
  in->read(&bsize, sizeof(unsigned int));
  vert.resize(vsize);
  face.resize(fsize);
  bord.resize(bsize);
  in->read(&*vert.begin(), vert.size() * sizeof(Point3f));
  in->read(&*face.begin(), face.size() * sizeof(unsigned short));
  in->read(&*bord.begin(), bord.size() * sizeof(Link));
}

void Fragment::Write(outstm *out) {
  out->write(&id, sizeof(unsigned int));
  out->write(&error, sizeof(float));
  
  unsigned int ssize = seeds.size();
  out->write(&ssize, sizeof(unsigned int));
  
  out->write(&*seeds.begin(), ssize * sizeof(Seed));
  out->write(&*seeds_id.begin(), ssize * sizeof(unsigned int));

  unsigned int psize = pieces.size();
  out->write(&psize, sizeof(unsigned int));
  
  for(unsigned int i = 0; i < pieces.size(); i++)
    pieces[i].Write(out);
}

void Fragment::Read(instm *in) {

  in->read(&id, sizeof(unsigned int));
  in->read(&error, sizeof(float));

  unsigned int ssize;
  in->read(&ssize, sizeof(unsigned int));
  seeds.resize(ssize);
  seeds_id.resize(ssize);
  in->read(&*seeds.begin(), ssize * sizeof(Seed));
  in->read(&*seeds_id.begin(), ssize * sizeof(unsigned int));

  unsigned int psize;
  in->read(&psize, sizeof(unsigned int));
  pieces.resize(psize);

  for(unsigned int i = 0; i < psize; i++)
    pieces[i].Read(in);
}

void nxs::Join(Fragment &in, 
	       vector<Point3f> &newvert,
	       vector<unsigned int> &newface,
	       vector<BigLink> &newbord) {

  map<unsigned int, unsigned int> patch_remap;
  vector<unsigned int> offsets;

  unsigned int totvert = 0;
  for(unsigned int i = 0; i < in.pieces.size(); i++) {
    offsets.push_back(totvert);
    patch_remap[in.pieces[i].patch] = i;
    totvert += in.pieces[i].vert.size();
  }

  vector<unsigned int> remap;
  remap.resize(totvert, 0xffffffff);

  //TODO what if totvert > 1<<22?
  //todo we really need a set?
  //  set<Link> newborders;
  unsigned int vcount = 0;
  unsigned int fcount = 0;
  unsigned int bcount = 0;

  for(unsigned int i = 0; i < in.pieces.size(); i++) {
    unsigned int offset = offsets[i];
    vector<Point3f> &vert = in.pieces[i].vert;
    vector<unsigned short> &face = in.pieces[i].face;
    vector<Link> &bord = in.pieces[i].bord;

    fcount += face.size()/3;

    for(unsigned int k = 0; k < vert.size(); k++) {
      assert(offset + k < remap.size());
      if(remap[offset + k] == 0xffffffff)
	remap[offset + k] = vcount++;
    }

    for(unsigned int k = 0; k < bord.size(); k++) {
      Link link = bord[k];
      if(link.IsNull()) continue;
      
      if(patch_remap.count(link.end_patch)) {//internal
	unsigned int idx = patch_remap[link.end_patch];
	unsigned int extoffset = offsets[idx];

	assert(extoffset + link.end_vert < remap.size());
	if(remap[extoffset + link.end_vert] == 0xffffffff)  //first time
	  remap[extoffset + link.end_vert] = remap[offset + link.start_vert];
      }
    }
  }
  assert(vcount < (1<<16));

  //L(a, b): Exist link between a, b
  //An external link L(e, v) where v belongs to the patches (and e not)
  //is valid only if: for every x in patches L(v, x) => L(e, x)
  //this means the number of internal links for the same shared
  //vertex is E = (n * (n-1)) where n is the number of duplicated vertices
  //and n must be the number of externa links.
  vector<unsigned int> internal_links;
  internal_links.resize(vcount, 0);
  
  map<BigLink, unsigned int> newborders;
  for(unsigned int i = 0; i < in.pieces.size(); i++) {
    unsigned int offset = offsets[i];
    vector<Link> &bord = in.pieces[i].bord;
    for(unsigned int k = 0; k < bord.size(); k++) {
      Link llink = bord[k];
      if(llink.IsNull()) continue;
      if(!patch_remap.count(llink.end_patch)) {//external...may be erased though
	BigLink link;
	link.start_vert = remap[offset + llink.start_vert];
	link.end_patch = llink.end_patch;
	link.end_vert = llink.end_vert;
	if(!newborders.count(link))
	  newborders[link] = 1;
	else
	  newborders[link]++;
      } else { //internal
	internal_links[remap[offset + llink.start_vert]]++;
      }
    }
  }

  newvert.resize(vcount);
  newface.resize(fcount*3);
  newbord.resize(0);
  
  fcount = 0;
  for(unsigned int i = 0; i < in.pieces.size(); i++) {
    unsigned int offset = offsets[i];
    vector<Point3f> &vert = in.pieces[i].vert;
    vector<unsigned short> &face = in.pieces[i].face;
    vector<Link> &bord = in.pieces[i].bord;
    
    for(unsigned int i = 0; i < vert.size(); i++) {            
      assert(offset + i < remap.size());
      assert(remap[offset + i] < vcount);
      newvert[remap[offset + i]] = vert[i];
    }
    
    for(unsigned int i = 0; i < face.size(); i++) {
      assert(offset + face[i] < remap.size());
      assert(remap[offset + face[i]] < newvert.size());
      assert(fcount < newface.size());
      newface[fcount++] = remap[offset + face[i]];
    }
  }  
  
  map<BigLink, unsigned int>::iterator b;
  for(b = newborders.begin(); b != newborders.end(); b++) {
    //test that number of links on this vertex is equal to
    //number of internal links of the internal vertex 
    const BigLink &link = (*b).first;
    unsigned int n = (*b).second;
    if(n * (n-1) == internal_links[link.start_vert])
      newbord.push_back(link);
  }
}

void nxs::Split(Fragment &out, 
		vector<Point3f> &newvert,
		vector<unsigned int> &newface,
		vector<BigLink> &newbord, 
		VoronoiPartition &part) {

  unsigned int nseeds = out.seeds.size();
  vector<Seed> &seeds = out.seeds;
  vector<unsigned int> &seeds_id = out.seeds_id;
  //preliminary count
  vector<unsigned int> count;
  count.resize(nseeds, 0);
  for(unsigned int f = 0; f < newface.size(); f += 3) {
    Point3f bari = (newvert[newface[f]] + 
		    newvert[newface[f+1]] + 
		    newvert[newface[f+2]])/3;
    unsigned int seed = out.Locate(bari);
    assert(seed < nseeds);
    count[seed]++;
  }

  //pruning small patches
  float min_size = (newface.size()/3) / 20.0f;
  vector<Seed> newseeds;
  vector<unsigned int> newseeds_id;

  for(unsigned int seed = 0; seed < nseeds; seed++) {
    if(count[seed] > min_size) {
      newseeds.push_back(seeds[seed]);
      newseeds_id.push_back(seeds_id[seed]);
    }
    if(count[seed] > (1<<16)) {
      cerr << "Ooops a cell came too big... quitting\n";
      exit(0);
    }
  }
  seeds = newseeds;
  seeds_id = newseeds_id;

  nseeds = seeds.size();

  //if != -1 remap global index to cell index (first arg)
  vector< vector<int> > vert_remap;
  vector< vector<int> > face_remap;

  vector<int> vert_count;
  vector<int> face_count;

  vert_remap.resize(nseeds);
  face_remap.resize(nseeds);
  vert_count.resize(nseeds, 0);
  face_count.resize(nseeds, 0);
  
  for(unsigned int seed = 0; seed < nseeds; seed++) 
    vert_remap[seed].resize(newvert.size(), -1);
  
  for(unsigned int f = 0; f < newface.size(); f += 3) {
    Point3f bari = (newvert[newface[f]] + 
		    newvert[newface[f+1]] + 
		    newvert[newface[f+2]])/3;
    
    unsigned int seed = out.Locate(bari);

    vector<int> &f_remap = face_remap[seed];

    f_remap.push_back(newface[f]);
    f_remap.push_back(newface[f+1]);
    f_remap.push_back(newface[f+2]);
    face_count[seed]++;
    

    vector<int> &v_remap = vert_remap[seed];
    
    for(int i = 0; i < 3; i++)
      if(v_remap[newface[f+i]] == -1)
	v_remap[newface[f+i]] = vert_count[seed]++;
  }

  //TODO assure no big ones.
  
  out.pieces.resize(nseeds);
  
  for(unsigned int seed = 0; seed != nseeds; seed++) { 
    NxsPatch &patch = out.pieces[seed];
    patch.patch = seeds_id[seed];
    
    //vertices first
    vector<int> &v_remap = vert_remap[seed];
   
    assert(vert_count[seed] > 0);
    vector<Point3f> &verts = patch.vert;
    verts.resize(vert_count[seed]);
    for(unsigned int i = 0; i < newvert.size(); i++) {
      if(v_remap[i] != -1)
	verts[v_remap[i]] = newvert[i];
    }
    
    //faces now
    vector<int> &f_remap = face_remap[seed];
    
    vector<unsigned short> &faces = patch.face;
    faces.resize(face_count[seed]*3);
    for(unsigned int i = 0; i < f_remap.size(); i++) {
      assert(v_remap[f_remap[i]] != -1);
      faces[i] = v_remap[f_remap[i]];
    }

    //borders last
    vector<Link> &bords = patch.bord;

    //process external borders 
    //for every esternal link we must update external patches!
    for(unsigned int i = 0; i < newbord.size(); i++) {
      BigLink link = newbord[i];
      if(v_remap[link.start_vert] == -1) continue;
      link.start_vert = v_remap[link.start_vert];
      assert(link.start_vert < (1<<16));
      Link llink;
      llink.start_vert = link.start_vert;
      llink.end_patch = link.end_patch;
      llink.end_vert = link.end_vert;
      bords.push_back(llink);
    }
			   
			   //process internal borders;
    //TODO higly inefficient!!!
    for(unsigned int rseed = 0; rseed < nseeds; rseed++) {
      if(seed == rseed) continue;

      vector<int> &vremapclose = vert_remap[rseed];
      for(unsigned int i = 0; i < newvert.size(); i++) {
	if(v_remap[i] != -1 && vremapclose[i] != -1) {
	  Link link;
	  link.end_patch = rseed + (1<<31);
	  link.start_vert = v_remap[i];
	  link.end_vert = vremapclose[i];
	  bords.push_back(link);
	}
      }
    }
  }
}

unsigned int Fragment::Locate(const Point3f &p) {
  float max_dist = 1e20;
  unsigned int id = 0xffffffff;
  for(unsigned int i = 0; i < seeds.size(); i++) {
    float dist = seeds[i].Dist(p);
    if(dist < max_dist) {
      max_dist = dist;
      id = i;
    }
  }
  assert(id != 0xffffffff);
  return id;
}
