#ifndef __VCGLIB_VERTEXPLUS_VN_TYPE
#define __VCGLIB_VERTEXPLUS_VN_TYPE
namespace vcg {

  template <class EdgeTemplate, class FaceTemplate> 
  class VertexVNf  : public VertexSimp2< VertexVNf, EdgeTemplate, FaceTemplate, vert::Normal3f, vert::Color4b, vert::Flag, vert::Coord3f> {};

  template <class EdgeTemplate, class FaceTemplate> 
  class VertexVNd  : public VertexSimp2< VertexVNd, EdgeTemplate, FaceTemplate, vert::Normal3d, vert::Color4b, vert::Flag, vert::Coord3d> {};

}

#endif
