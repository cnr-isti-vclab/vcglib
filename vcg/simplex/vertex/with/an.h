#ifndef __VCGLIB_VERTEX__AN__TYPE
#define __VCGLIB_VERTEX__AN__TYPE


#define VERTEX_TYPE VertexAN 

#define __VCGLIB_VERTEX_A
#define __VCGLIB_VERTEX_N


#include <vcg/simplex/vertex/base.h> 


#undef VERTEX_TYPE 
#undef __VCGLIB_VERTEX_A
#undef __VCGLIB_VERTEX_N


namespace vcg {
typedef VertexAN<float>  VertexANf;
typedef VertexAN<double> VertexANd;
}

#endif