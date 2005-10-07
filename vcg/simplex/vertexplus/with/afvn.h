#ifndef __VCGLIB_VERTEXPLUS_AFVN_TYPE
#define __VCGLIB_VERTEXPLUS_AFVN_TYPE
namespace vcg {

  template <class EdgeTemplate, class FaceTemplate> 
  class VertexAFVNd  : public VertexSimp2< VertexAFVNd, EdgeTemplate, FaceTemplate, vert::Normal3d, vert::VFAdj, vert::Flag, vert::Coord3d> {};
  
  template <class EdgeTemplate, class FaceTemplate> 
  class VertexAFVNf  : public VertexSimp2< VertexAFVNf, EdgeTemplate, FaceTemplate, vert::Normal3f, vert::VFAdj, vert::Flag, vert::Coord3f> {};
}
#endif
