#ifndef __VCGLIB_VERTEX__VQ__TYPE
#define __VCGLIB_VERTEX__VQ__TYPE


#define VERTEX_TYPE VertexVQ 

#define __VCGLIB_VERTEX_VQ

#include <vcg/simplex/vertex/base.h> 


#undef VERTEX_TYPE 

#undef __VCGLIB_VERTEX_VQ

namespace vcg {
typedef VertexVQ<float>  VertexVQf;
typedef VertexVQ<double> VertexVQd;
}

#endif
