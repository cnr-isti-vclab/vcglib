#ifndef __VCGLIB_VERTEX__N__TYPE
#define __VCGLIB_VERTEX__N__TYPE

#define VERTEX_TYPE VertexN 

#define __VCGLIB_VERTEX_N

#include <vcg/simplex/vertex/base.h> 

#undef VERTEX_TYPE 

#undef __VCGLIB_VERTEX_N

using namespace vcg;

typedef VertexN<float>  VertexNf;
typedef VertexN<double> VertexNd;

#endif