#ifndef __VCGLIB_VERTEX_AN__TYPE
#define __VCGLIB_VERTEX_AN__TYPE
namespace vcg {
  template <class FaceTemplate> 
  class VertexANd  : public VertexSimp2< VertexANd, DumET, FaceTemplate, vert::Normal3d, vert::VFAdj, vert::Flag, vert::Coord3d> {};
  template <class FaceTemplate> 
  class VertexANf  : public VertexSimp2< VertexANf, DumET, FaceTemplate, vert::Normal3f, vert::VFAdj, vert::Flag, vert::Coord3f> {};
}
#endif
