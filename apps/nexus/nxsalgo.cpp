#include <vector>
#include <map>
#include <set>
#include <iostream>

//#include <wrap/strip/tristrip.h>

#include "nxsalgo.h"
#include "vfile.h"
#include "nexus.h"
#include "watch.h"

using namespace std;
using namespace nxs;
using namespace vcg;

#include "tristripper/tri_stripper.h"
using namespace triangle_stripper;

void nxs::ComputeNormals(Nexus &nexus) {
  assert(nexus.signature & NXS_NORMALS_SHORT ||
	 nexus.signature & NXS_NORMALS_FLOAT);

  //setting borders readonly:

  assert(!nexus.borders.IsReadOnly());
  nexus.borders.SetReadOnly(true);
  
  bool use_short = (nexus.signature & NXS_NORMALS_SHORT) != 0;

  //TODO use a temporary file to store border normals
  unsigned int tmpb_offset = 0;
  vector<unsigned int> tmpb_start;
  VFile<Point3f> tmpb;
  if(!tmpb.Create("tmpb.tmp")) {
    cerr << "Could not create temporary border file\n";
    exit(0);
  }

  for(unsigned int p = 0; p < nexus.size(); p++) {
    Border border = nexus.GetBorder(p);
    tmpb_start.push_back(tmpb_offset);
    tmpb_offset += border.Size();
  }

  Point3f zero(0.0f, 0.0f, 0.0f);

  tmpb.Resize(tmpb_offset);
  for(unsigned int i = 0; i < tmpb.Size(); i++)
    tmpb[i] = zero;
  tmpb.Flush();
  
  //first step normals in the same patch.
  cerr << "First Step\n";
  Report report(nexus.size(), 5);
  vector<Point3f> normals;

  for(unsigned int p = 0; p < nexus.size(); p++) {
    report.Step(p);
    Patch &patch = nexus.GetPatch(p);
    
    normals.clear();  
    normals.resize(patch.nv, Point3f(0, 0, 0));
    

    if(nexus.signature & NXS_FACES) 
      for(unsigned int i = 0; i < patch.nf; i++) {
	      unsigned short *f = patch.Face(i);
	      Point3f &v0 = patch.Vert(f[0]);
	      Point3f &v1 = patch.Vert(f[1]);
	      Point3f &v2 = patch.Vert(f[2]);
	
	      Point3f norm = (v1 - v0) ^ (v2 - v0); 
	      normals[f[0]] += norm;
	      normals[f[1]] += norm;
	      normals[f[2]] += norm;
      }

    if(nexus.signature & NXS_STRIP) 
      for(int i = 0; i < patch.nf - 2; i++) {
	      unsigned short *f = patch.FaceBegin() + i;
	      Point3f &v0 = patch.Vert(f[0]);
	      Point3f &v1 = patch.Vert(f[1]);
	      Point3f &v2 = patch.Vert(f[2]);
	
	      Point3f norm = (v1 - v0) ^ (v2 - v0); 
	      if(i%2) norm = -norm;
	      normals[f[0]] += norm;
	      normals[f[1]] += norm;
	      normals[f[2]] += norm;
      }
    
    if(use_short) {
      for(unsigned int i = 0; i < patch.nv; i++) {
	      Point3f &norm = normals[i];
	      norm.Normalize();
	      short *n = patch.Norm16(i);
	      for(int k = 0; k < 3; k++) {
	        n[k] = (short)(norm[k] * 32766);
	      }
	      n[3] = 0;
      }
    } else {
      memcpy(patch.Norm16Begin(), &*normals.begin(), 
	     normals.size() * sizeof(Point3f));
    }


    Border border = nexus.GetBorder(p);

    
    map<unsigned int, map<unsigned short, Point3f> > bnorm;


    unsigned int poff = tmpb_start[p];
    for(unsigned int i = 0; i < border.Size(); i++) {
      Link &link = border[i];
      if(link.IsNull()) continue;
      Point3f pt = normals[link.start_vert];
      //bnorm[p][link.start_vert] = pt;
      bnorm[link.end_patch][link.end_vert] = pt;
      tmpb[poff + i] += pt;                          
    }

    map<unsigned int, map<unsigned short, Point3f> >::iterator k;
    for(k = bnorm.begin(); k != bnorm.end(); k++) {
      unsigned int patch = (*k).first;
      Border border = nexus.GetBorder(patch);
      unsigned int offset = tmpb_start[patch];
      for(unsigned int i = 0; i < border.Size(); i++) {
	      Link &link = border[i];
	      //assert(!link.IsNull());
	      //TODO not accurate
        if(link.end_patch != p) continue;
	      if((*k).second.count(link.start_vert))
	        tmpb[offset + i] += (*k).second[link.start_vert];
      }
    }

    /*    set<unsigned int> close;
    for(unsigned int i = 0; i < border.Size(); i++) {
      Link &link = border[i];
      if(link.IsNull()) continue;
      unsigned int off = tmpb_start[p];
      Point3f p = tmpb.read(off + i);
      p += normals[link.start_vert];
      tmpb.write(off + i, p);
      //      tmpb[off + i] += normals[link.start_vert];
      close.insert(link.end_patch);
    }

    set<unsigned int>::iterator k;
    for(k = close.begin(); k != close.end(); k++) {
      Border remote = nexus.GetBorder(*k);
      unsigned int off = tmpb_start[*k];

      for(unsigned int i = 0; i < remote.Size(); i++) {
	Link &link = remote[i];
	if(link.IsNull()) continue;
	if(link.end_patch != p) continue;
	Point3f p = tmpb.read(off + i);
	p += normals[link.end_vert];
	tmpb.write(off + i, p);
	//	tmpb[off + i] += normals[link.end_vert];
      }
      }*/
  }

  //Second step unify normals across borders
  cerr << "Second step\n";
  report.Init(nexus.size());
  for(unsigned int p = 0; p < nexus.size(); p++) {
    report.Step(p);
    Patch &patch = nexus.GetPatch(p);
    Border border = nexus.GetBorder(p);

    for(unsigned int i = 0; i < border.Size(); i++) {
      Link &link = border[i];
      if(link.IsNull()) continue;
      unsigned int off = tmpb_start[p];
      //      Point3f &n = tmpb[off + i];
      Point3f n = tmpb[off + i];
      n.Normalize();
      if(use_short) {
	n *= 32766;
	short *np = patch.Norm16(link.start_vert);
	np[0] = (short)n[0];
	np[1] = (short)n[1];
	np[2] = (short)n[2];
	np[3] = 0;
      } else {
	patch.Norm32(link.start_vert) = n;
      }
    }
  }
  tmpb.Close();
  tmpb.Delete();
  //TODO remove temporary file.
  nexus.borders.SetReadOnly(false);
}

void nxs::ComputeTriStrip(unsigned short nfaces, unsigned short *faces, 
		  vector<unsigned short> &strip) {

  
  vector<unsigned int> index;
  index.resize(nfaces*3);
  for(int i = 0; i < nfaces*3; i++) {
    index[i] = faces[i];
  }
  int cache_size = 0;
  tri_stripper stripper(index);
  stripper.SetCacheSize(cache_size);		
  // = 0 will disable the cache optimizer
  stripper.SetMinStripSize(0);
  tri_stripper::primitives_vector primitives;
  stripper.Strip(&primitives);

  if(primitives.back().m_Indices.size() < 3) {
    primitives.pop_back();
  }
  //TODO spostare questo dentro il ciclo che rimonta le strip.
  if(primitives.back().m_Type == tri_stripper::PT_Triangles) {
    tri_stripper::primitives p;
    p = primitives.back();
    primitives.pop_back();		
    for(unsigned int i = 0; i < p.m_Indices.size(); i += 3) {
      tri_stripper::primitives s;
      s.m_Type = tri_stripper::PT_Triangle_Strip;
      s.m_Indices.push_back(p.m_Indices[i]);
      s.m_Indices.push_back(p.m_Indices[i+1]);
      s.m_Indices.push_back(p.m_Indices[i+2]);
      primitives.push_back(s);
    }
  }
  
  for(unsigned int i = 0; i < primitives.size(); i++) {
    tri_stripper::primitives &primitive = primitives[i];
    assert(primitive.m_Indices.size() != 0);
    int len = primitive.m_Indices.size();
    for(int l = 0; l < len; l++)  		
      strip.push_back(primitive.m_Indices[l]);
    
    
    if(i < primitives.size()-1) { //not the last primitive.
      strip.push_back(primitive.m_Indices[len-1]);
      //TODO optimize this!
      if((len%2) == 1) 	//do not change orientation....
	strip.push_back(primitive.m_Indices[len-1]);
      strip.push_back(primitives[i+1].m_Indices[0]);
    }			
  }
}

void nxs::Reorder(unsigned int signature, Patch &patch) {
  vector<unsigned> remap;
  remap.resize(patch.nv, 0xffff);
  
  int nf = patch.nf;
  if(signature & NXS_FACES)
    nf *= 3;
  
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
  memcpy(&*vert.begin(), patch.VertBegin(), patch.nv * sizeof(Point3f));
  for(int i = 0; i < patch.nv; i++)
    patch.Vert(remap[i]) = vert[i];
}

//TODO actually use threshold
void nxs::Unify(Nexus &nexus, float threshold) {
  //TODO what if colors or normals or strips?
  unsigned int duplicated = 0;
  unsigned int degenerate = 0;

  for(unsigned int p = 0; p < nexus.size(); p++) {
    Entry &entry = nexus[p];
    Patch &patch = nexus.GetPatch(p);

    unsigned int vcount = 0;
    map<Point3f, unsigned short> vertices;

    vector<unsigned short> remap;
    remap.resize(patch.nv);

    for(unsigned int i = 0; i < patch.nv; i++) {
      Point3f &point = patch.Vert(i);

      if(!vertices.count(point)) 
        vertices[point] = vcount++;
      else 
        duplicated++;

      remap[i] = vertices[point];
    }
    assert(vertices.size() <= patch.nv);
    if(vertices.size() == patch.nv) //no need to unify
      continue;

    vector<Point3f> newvert;
    newvert.resize(vertices.size());
    map<Point3f, unsigned short>::iterator k;
    for(k = vertices.begin(); k != vertices.end(); k++) {
      newvert[(*k).second] = (*k).first;
    }

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

    memcpy(patch.VertBegin(), &(newvert[0]), entry.nvert*sizeof(Point3f));
    memcpy(patch.FaceBegin(), &(newface[0]), entry.nface*3*sizeof(short));

    //testiamo il tutto...  TODO remove this of course
    for(unsigned int i =0; i < patch.nf; i++) {
      for(int k =0 ; k < 3; k++)
        if(patch.Face(i)[k] >= patch.nv) {
          cerr <<" Unify has problems\n";
          exit(0);
        }
    }
    
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
    Border border = nexus.GetBorder(p);
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

