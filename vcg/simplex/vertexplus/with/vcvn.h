#ifndef __VCGLIB_VERTEXPLUS_VCVN_TYPE
#define __VCGLIB_VERTEXPLUS_VCVN_TYPE
namespace vcg {

  template <class EdgeTemplate, class FaceTemplate> 
  class VertexVCVNf  : public VertexSimp1< VertexVCVNf, DumET, vert::Normal3f, vert::Color4b, vert::Flag, vert::Coord3f> {};
  template <class EdgeTemplate, class FaceTemplate> 
  class VertexVCVNd  : public VertexSimp1< VertexVCVNd, DumET, vert::Normal3d, vert::Color4b, vert::Flag, vert::Coord3d> {};

}

#endif
