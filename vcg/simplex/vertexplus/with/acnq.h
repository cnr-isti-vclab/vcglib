#ifndef __VCGLIB_VERTEX_ACNQ_TYPE
#define __VCGLIB_VERTEX_ACNQ_TYPE
namespace vcg {
  template <class FaceTemplate> 
  class VertexACNQd  : public VertexSimp2< VertexACNQd, DumET, FaceTemplate, vert::Normal3d, vert::VFAdj, vert::Flag, vert::Coord3d> {};
  template <class FaceTemplate> 
  class VertexACNQf  : public VertexSimp2< VertexACNQf, DumET, FaceTemplate, vert::Normal3f, vert::VFAdj, vert::Flag, vert::Coord3f> {};
}
#endif
