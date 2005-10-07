#ifndef __VCGLIB_VERTEXPLUS_AFVCVNVQ_TYPE
#define __VCGLIB_VERTEXPLUS_AFVCVNVQ_TYPE
namespace vcg {
  template <class FaceTemplate> 
  class VertexAFVCVNVQd  : public VertexSimp2< VertexAFVCVNVQd, DumET, FaceTemplate, vert::Normal3d, vert::VFAdj, vert::Flag, vert::Coord3d> {};
  template <class FaceTemplate> 
  class VertexAFVCVNVQf  : public VertexSimp2< VertexAFVCVNVQf, DumET, FaceTemplate, vert::Normal3f, vert::VFAdj, vert::Flag, vert::Coord3f> {};
}
#endif
