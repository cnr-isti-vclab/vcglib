#include "fragment.h"
#include "border.h"
#include "pvoronoi.h"
#include <iostream>

using namespace std;
using namespace vcg;
using namespace nxs;
using namespace pt;

void NxsPatch::write(outstm *out) {
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

void NxsPatch::read(instm *in) {
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

void Fragment::write(outstm *out) {
  out->write(&id, sizeof(unsigned int));
  out->write(&error, sizeof(float));

  unsigned int esize = update.erased.size();
  unsigned int csize = update.created.size();
  out->write(&esize, sizeof(unsigned int));
  out->write(&csize, sizeof(unsigned int));
  out->write(&*update.erased.begin(), esize * sizeof(unsigned int));
  out->write(&*update.created.begin(), csize * sizeof(unsigned int));

  unsigned int psize = pieces.size();
  out->write(&psize, sizeof(unsigned int));
  
  for(unsigned int i = 0; i < pieces.size(); i++)
    pieces[i].write(out);
}

void Fragment::read(instm *in) {

  in->read(&id, sizeof(unsigned int));
  in->read(&error, sizeof(float));

  unsigned int esize;
  unsigned int csize;
  in->read(&esize, sizeof(unsigned int));
  in->read(&csize, sizeof(unsigned int));
  update.erased.resize(esize);
  update.created.resize(csize);
  in->read(&*update.erased.begin(), esize * sizeof(unsigned int));
  in->read(&*update.created.begin(), csize * sizeof(unsigned int));

  unsigned int psize;
  in->read(&psize, sizeof(unsigned int));
  pieces.resize(psize);

  for(unsigned int i = 0; i < psize; i++)
    pieces[i].read(in);
}

void nxs::join(Fragment &in, 
	       vector<Point3f> &newvert,
	       vector<unsigned int> &newface,
	       vector<Link> &newbord) {

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

  //TODO what if totvert > 32768?
  cerr << "Totvert " << totvert << endl;
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
    /*  assert(link.start_vert < remap.size());
      assert(remap[offset + link.start_vert] != 0xffffffff);
      if(!patch_remap.count(link.end_patch)) { //external
	//test if erased in history... in wich case we do not add border
	//	if(!erased.count(link.end_patch)) { ?/????
	link.start_vert = remap[offset + link.start_vert];
	newborders.insert(link);
      } else {
	//internal
	//TODO unsigned int &rmpv = remap[link.end_patch][link.end_vert];
	unsigned int idx = patch_remap[link.end_patch];
	unsigned int extoffset = offsets[idx];
	if(extoffset + link.end_vert >= remap.size()) {
	  cerr << "extoffset: " << extoffset << endl;
	  cerr << "end_V: " << link.end_vert << endl;
	  cerr << "remsz: " << remap.size() << endl;
	  for(unsigned int i = 0; i < in.pieces.size(); i++)
	    cerr << "size: " << i << " =" << in.pieces[i].vert.size() << endl;
	}
	assert(extoffset + link.end_vert < remap.size());
	if(remap[extoffset + link.end_vert] == 0xffffffff)  //first time
	  remap[extoffset + link.end_vert] = remap[offset + link.start_vert];
      }
      }*/
  }

  //L(a, b): Exist link between a, b
  //An external link L(e, v) where v belongs to the patches (and e not)
  //is valid only if: for every x in patches L(v, x) => L(e, x)
  //this means the number of internal links for the same shared
  //vertex is E = (n * (n-1)) where n is the number of duplicated vertices
  //and n must be the number of externa links.
  vector<unsigned int> internal_links;
  internal_links.resize(vcount, 0);
  
  map<Link, unsigned int> newborders;
  for(unsigned int i = 0; i < in.pieces.size(); i++) {
    unsigned int offset = offsets[i];
    vector<Link> &bord = in.pieces[i].bord;
    for(unsigned int k = 0; k < bord.size(); k++) {
      Link link = bord[k];
      if(link.IsNull()) continue;
      if(!patch_remap.count(link.end_patch)) {//external...may be erased though
	link.start_vert = remap[offset + link.start_vert];
	if(!newborders.count(link))
	  newborders[link] = 1;
	else
	  newborders[link]++;
      } else { //internal
	internal_links[remap[offset + link.start_vert]]++;
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
  
  map<Link, unsigned int>::iterator b;
  for(b = newborders.begin(); b != newborders.end(); b++) {
    //test that number of links on this vertex is equal to
    //number of internal links of the internal vertex 
    const Link &link = (*b).first;
    unsigned int n = (*b).second;
    if(n * (n-1) == internal_links[link.start_vert])
      newbord.push_back(link);
  }
}

void nxs::split(Fragment &out, 
		vector<Point3f> &newvert,
		vector<unsigned int> &newface,
		vector<Link> &newbord, 
		VoronoiPartition &part) {
  
}
 
