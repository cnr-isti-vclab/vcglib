#include <iostream>

#include "nexus.h"  
#include "watch.h"

using namespace vcg;
using namespace std;
using namespace nxs;

int main(int argc, char *argv[]) {
  if(argc != 2) {
    cerr << "Usage: " << argv[0] << " <nexusfile>\n";
    return -1;
  }
  Nexus nexus;
  if(!nexus.Load(argv[1], true)) {
    cerr << "Could not open file: " << argv[1] << endl;
    return -1;
  }
  Report report(nexus.size());
  for(unsigned int patchid = 0; patchid < nexus.size(); patchid++) {    
    report.Step(patchid);
    Entry &info = nexus[patchid];
    Patch &patch = nexus.GetPatch(patchid);
    for(int f = 0; f < patch.nf; f++) {
      unsigned short *face = patch.Face(f);
      for(int k = 0; k < 3; k++) {
        if(face[k] > patch.nv) {
          cerr << "Invalid face number: " << face[k] << " > " << patch.nv << endl;
          cerr << "At patch: " << patchid << endl;
          //exit(0);
        }
      }      
    }
    Sphere3f &sphere = info.sphere;
    for(int v = 0; v < patch.nv; v++) {
      Point3f &p = patch.Vert(v);
      float dist = Distance(sphere, p);
      if(dist > 0.001) {
      //if(!info.sphere.IsIn(p)) {
        cerr << "Vertex outside bound: (" << p[0] << " " << p[1] << " " << p[2] << ")\n";
        Point3f &c = sphere.Center();
        cerr << "Sphere: (" << c[0] << " " << c[1] << " " << c[2] << ") R: " << sphere.Radius() << endl;;                
        cerr << "Distance: " << dist << endl;
        cerr << "At patch: " << patchid << endl;

      }
    }
  }
  report.Finish();

  cerr << "Testing borders\n";

  for(unsigned int patchid = 0; patchid < nexus.size(); patchid++) {
    Entry &info = nexus[patchid];
    Border &border = nexus.GetBorder(patchid);
    for(unsigned int i = 0; i < border.Size(); i++) {
      Link &link = border[i];
      if(link.start_vert == 0 && link.end_vert == 0 && link.end_patch == 0) {
        cerr << "patch: " << patchid << " corrupted memory?" << endl;
      }
      if(link.IsNull()) {
        cerr << "Null link: " << i << " at patch: " << patchid << endl;
        exit(0);
      }
      if(link.end_patch < 0 || link.end_patch >= nexus.size()) {
        cerr << "Invalid link end patch: " << link.end_patch << " at patch: " << patchid << endl;
        exit(0);
      }
      if(link.start_vert > info.nvert) {
        cerr << "Invalid link start_vert: " << link.start_vert << " at patch: " << patchid << endl;
        exit(0);
      }
      if(link.end_vert > nexus[link.end_patch].nvert) {
        cerr << "Invalid link end vert: " << link.end_vert << " at patch: " << patchid << endl;
        exit(0);
      }        
    }
  }

  cerr << "Reciprocity borders test\n";
  for(unsigned int patchid = 0; patchid < nexus.size(); patchid++) {
    Entry &info = nexus[patchid];
    Border &border = nexus.GetBorder(patchid);
    vector<Link> links;
    links.resize(border.Size());
    memcpy(&*links.begin(),&(border[0]),links.size() * sizeof(Link));    
    for(unsigned int i = 0; i < links.size(); i++) {      
      Link &link = links[i];      
      Border &rborder = nexus.GetBorder(link.end_patch, false);          

      bool found = false;
      for(unsigned int k = 0; k < rborder.Size(); k++) {
        Link rlink = rborder[k];
        if(rlink.end_patch == patchid) {

          if(rlink.end_vert == link.start_vert) {
            if(rlink.start_vert != link.end_vert) {
              cerr << "Something wrong with links!\n";
              exit(0);
            }
            found = true;
            break;
          }
          if(rlink.start_vert == link.end_vert) {
            if(rlink.end_vert != link.start_vert) {
              cerr << "Something wrong with links!\n";
              exit(0);
            }
            found = true;
            break;
          }
        }
      }
      if(!found) {
        cerr << "A link is one way from patch: " << patchid << " vert: " << link.start_vert
          << " to patch: " << link.end_patch << " vert: " << link.end_vert << endl;
        for(unsigned int t = 0; t < rborder.Size(); t++) {
          Link &rlink = rborder[t];
          cerr << rlink.start_vert << " -> p: " << rlink.end_patch << " v: " << rlink.end_vert << endl;
        }
        exit(0);
      }
    }
  }
  
  return 0;
}
