#ifndef __VCGLIB_VERTEX__VN__TYPE
#define __VCGLIB_VERTEX__VN__TYPE


#define VERTEX_TYPE VertexVN 

#define __VCGLIB_VERTEX_VN

#include <vcg/simplex/vertex/base.h> 


#undef VERTEX_TYPE 
#undef __VCGLIB_VERTEX_VN


namespace vcg {
typedef VertexVN<float>  VertexVNf;
typedef VertexVN<double> VertexVNd;
}

#endif
