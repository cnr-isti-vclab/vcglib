#ifndef __VCGLIB_VERTEX__CN__TYPE
#define __VCGLIB_VERTEX__CN__TYPE


#define VERTEX_TYPE VertexCN 

#define __VCGLIB_VERTEX_N
#define __VCGLIB_VERTEX_C

#include <vcg/simplex/vertex/base.h> 


#undef VERTEX_TYPE 
#undef __VCGLIB_VERTEX_N
#undef __VCGLIB_VERTEX_C


namespace vcg {
typedef VertexCN<float>  VertexCNf;
typedef VertexCN<double> VertexCNd;
}

#endif
