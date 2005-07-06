#ifndef __VCGLIB_VERTEX__VNVQ__TYPE
#define __VCGLIB_VERTEX__VNVQ__TYPE


#define VERTEX_TYPE VertexVNVQ 

#define __VCGLIB_VERTEX_VN
#define __VCGLIB_VERTEX_VQ

#include <vcg/simplex/vertex/base.h> 


#undef VERTEX_TYPE 
#undef __VCGLIB_VERTEX_VN
#undef __VCGLIB_VERTEX_VQ


namespace vcg {
typedef VertexVNVQ<float>  VertexVNVQf;
typedef VertexVNVQ<double> VertexVNVQd;
}

#endif
