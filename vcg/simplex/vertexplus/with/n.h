#ifndef __VCGLIB_VERTEX__N__TYPE
#define __VCGLIB_VERTEX__N__TYPE
namespace vcg {

  class VertexNd  : public VertexSimp1< VertexNd, DumET, vert::Normal3d, vert::Flag, vert::Coord3d> {};
  class VertexNf  : public VertexSimp1< VertexNf, DumET, vert::Normal3f, vert::Flag, vert::Coord3f> {};

}
#endif
