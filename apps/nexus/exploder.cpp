#include <stdlib.h>
#include <iostream>

#include "nexus.h"

using namespace nxs;
using namespace vcg;
using namespace std;

Point3f explode(Point3f &v, Point3f &center, Point3f &dir, float mag) {
  v = ((v - center) * mag + center) + dir;
  return v;
}

int main(int argc, char *argv[]) {
  if(argc < 3) {
    cerr << "Usage: " << argv[0] << " <input> <output> [factor]\n";
    return -1;
  }
  string input = argv[1];
  string output = argv[2];
  float factor = 2;
  if(argc == 4) {
    factor = atoi(argv[3]);
    if(factor == 0) {
      cerr << "Invalid factor: " << argv[3] << endl;
      return -1;
    }
  }
  Nexus in;
  if(!in.Load(input, true)) {
    cerr << "Could not load nexus: " << input << endl;
    return -1;
  }
  Nexus out;
  if(!out.Create(output, in.signature, in.chunk_size)) {
    cerr << "Could not create nexus: " << output << endl;
    return -1;
  }

  out.sphere = in.sphere;
  out.history = in.history;
  for(unsigned int i = 0; i < in.index.size(); i++) {
    unsigned int patch = i;
    Nexus::PatchInfo &src_entry = in.index[patch];
    Patch &src_patch = in.GetPatch(patch);
    Border src_border = in.GetBorder(patch);
    
    out.AddPatch(src_entry.nvert, src_entry.nface, src_border.Available());
  
    Nexus::PatchInfo &dst_entry = out.index[patch];

    Patch dst_patch = out.GetPatch(patch);
    
    Point3f dir = src_entry.sphere.Center() - in.sphere.Center();
    dir.Normalize();
    dir *= 10 * src_entry.error;

    memcpy(dst_patch.VertBegin(), src_patch.VertBegin(), 
	   src_patch.nv * sizeof(Point3f));
    
    Point3f *ptr = dst_patch.VertBegin();
    for(int i = 0; i < dst_patch.nv; i++) {
      ptr[i] =  explode(ptr[i], src_entry.sphere.Center(), dir, 
			0.5 *src_entry.error);
    } 

    
    if(in.signature & NXS_STRIP) {
      memcpy(dst_patch.FaceBegin(), src_patch.FaceBegin(), 
	     src_patch.nf * sizeof(unsigned short));
    } else {
      memcpy(dst_patch.FaceBegin(), src_patch.FaceBegin(), 
	     src_patch.nf * sizeof(unsigned short) * 3);
    }
  
    if((in.signature & NXS_COLORS) && (out.signature & NXS_COLORS))
      memcpy(dst_patch.ColorBegin(), src_patch.ColorBegin(), 
	     src_patch.nv * sizeof(unsigned int));
    
    if((in.signature & NXS_NORMALS_SHORT) && 
       (out.signature & NXS_NORMALS_SHORT))
      memcpy(dst_patch.Norm16Begin(), src_patch.Norm16Begin(), 
	     src_patch.nv * sizeof(short)*4);
    
    //reordering
    //WATCH OUT BORDERS!
    //    Reorder(out.signature, dst_patch);
    //copying entry information;
    dst_entry.sphere = src_entry.sphere;
    dst_entry.error = src_entry.error;
    
    //adding borders.
    for(unsigned int i = 0; i < src_border.Size(); i++) {
      Link &link = src_border[i];
      if(link.IsNull()) continue;
      assert(link.end_patch < in.index.size());
    }
    Border dst_border = out.GetBorder(patch);
    out.borders.ResizeBorder(patch, src_border.Size());
    memcpy(dst_border.Start(), src_border.Start(), 
	   src_border.Size() * sizeof(Link));
  }
  in.Close(); 
  out.Close();
  return 0;
}
