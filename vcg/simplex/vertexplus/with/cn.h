#ifndef __VCGLIB_VERTEX__CN__TYPE
#define __VCGLIB_VERTEX__CN__TYPE
namespace vcg {

  class VertexCNf  : public VertexSimp1< VertexCNf, DumET, vert::Normal3f, vert::Color4b, vert::Flag, vert::Coord3f> {};
  class VertexCNd  : public VertexSimp1< VertexCNd, DumET, vert::Normal3d, vert::Color4b, vert::Flag, vert::Coord3d> {};

}

#endif
