#ifndef NXS_EXPORT_H
#define NXS_EXPORT_H

#include <map>
#include <vector>
#include "nexus.h"
#include "extraction.h"

namespace nxs {

template <class MESH> void ExportTriMesh(Nexus &nexus, 
					 vector<unsigned int> &patches, 
					 MESH &mesh) {


  typedef typename MESH::VertexType VertexType;
  typedef typename MESH::FaceType FaceType;
  typedef typename MESH::ScalarType ScalarType;
  //for every patch  record the global position of the vertices
  std::map<unsigned int, std::vector<unsigned int> > remap;

  //Resize all remapping vectors
  for(unsigned int p = 0; p < patches.size(); p++) {
    unsigned int npatch = patches[p];
    remap[npatch].resize(nexus[npatch].nvert, 0xffffffff);
  }

  //Filling remapping vectors
  unsigned int vcount = 0;
  for(unsigned int p = 0; p < patches.size(); p++) {
    unsigned int npatch = patches[p];
    Patch &patch = nexus.GetPatch(npatch);
    
    std::vector<unsigned int> &rmp = remap[npatch];
    
    for(unsigned int v = 0; v < patch.nv; v++) 
      if(rmp[v] == 0xffffffff) 
	rmp[v] = vcount++;
    
    Border &border = nexus.GetBorder(npatch);
    for(unsigned int k = 0; k < border.Size(); k++) {
      Link link = border[k];
      if(link.IsNull()) continue;
      
      if(remap.count(link.end_patch)) //internal
	if(remap[link.end_patch][link.end_vert] == 0xffffffff)
	  remap[link.end_patch][link.end_vert] = rmp[link.start_vert];
    }
  }

  mesh.vert.resize(vcount);
  mesh.VertexNumber() = vcount;
  //copying vectors and faces
  for(unsigned int p = 0; p < patches.size(); p++) {
    unsigned int npatch = patches[p];
    Patch &patch = nexus.GetPatch(npatch);

    std::vector<unsigned int> &rmp = remap[npatch];

    //coping vertices
    VertexType vertex;
    vertex.ClearFlags();
    for(unsigned int v = 0; v < patch.nv; v++) {
      vertex.P()[0] = (ScalarType)patch.Vert3f(v)[0];
      vertex.P()[1] = (ScalarType)patch.Vert3f(v)[1];
      vertex.P()[2] = (ScalarType)patch.Vert3f(v)[2];
      if(mesh.HasPerVertexNormal()) {
	if(nexus.signature.vnorm == Encodings::SHORT4) {
	  vertex.N()[0] = (ScalarType)((short *)patch.VNormBegin())[v*4];
	  vertex.N()[1] = (ScalarType)((short *)patch.VNormBegin())[v*4 +1];
	  vertex.N()[2] = (ScalarType)((short *)patch.VNormBegin())[v*4 +2];
	} else if(nexus.signature.vnorm == Encodings::FLOAT3) {
	  vertex.N()[0] = (ScalarType)((float *)patch.VNormBegin())[v*3];
	  vertex.N()[0] = (ScalarType)((float *)patch.VNormBegin())[v*3 +1];
	  vertex.N()[0] = (ScalarType)((float *)patch.VNormBegin())[v*3 +2];
	} else if(nexus.signature.vnorm) {
	  //i should write other exporters
	  assert(0);
	}
      }
      if(mesh.HasPerVertexColor() && nexus.signature.vcolor) {
	if(nexus.signature.vcolor == Encodings::BYTE4) {
	  vertex.C()[0] = ((unsigned char *)patch.VColorBegin())[v*4];
	  vertex.C()[1] = ((unsigned char *)patch.VColorBegin())[v*4+1];
	  vertex.C()[2] = ((unsigned char *)patch.VColorBegin())[v*4+2];
	  vertex.C()[3] = ((unsigned char *)patch.VColorBegin())[v*4+3];
	} else if(nexus.signature.vcolor) {
	  //i should write other exporters
	  assert(0);
	}
      }
      if(mesh.HasPerVertexTexture() && nexus.signature.vtext) {
	//i should write other exporters
	assert(0);
      }
      assert(rmp[v] < mesh.vert.size());
      mesh.vert[rmp[v]] = vertex;
    }
    //remap faces now
    FaceType face;
    if(nexus.signature.face == Signature::TRIANGLES) {
      for(unsigned int f = 0; f < patch.nf; f++) {
	face.V(0) = &mesh.vert[rmp[patch.Face(f)[0]]];
	face.V(1) = &mesh.vert[rmp[patch.Face(f)[1]]];
	face.V(2) = &mesh.vert[rmp[patch.Face(f)[2]]];
	mesh.face.push_back(face);
	mesh.SimplexNumber()++;
  /*static bool HasPerFaceColor()  {return FaceType::HasFaceColor() ;}
    static bool HasPerFaceNormal() {return FaceType::HasFaceNormal();}
    static bool HasPerFaceMark()   {return FaceType::HasFaceMark()  ;}
    static bool HasPerFaceQuality(){return FaceType::HasFaceQuality();}*/
      }
    }
  }
}

}//namespace

#endif
