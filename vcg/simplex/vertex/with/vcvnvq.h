#ifndef __VCGLIB_VERTEX__VCVNVQ__TYPE
#define __VCGLIB_VERTEX__VCVNVQ__TYPE


#define VERTEX_TYPE VertexVCVNVQ 

#define __VCGLIB_VERTEX_VC
#define __VCGLIB_VERTEX_VN
#define __VCGLIB_VERTEX_VQ

#include <vcg/simplex/vertex/base.h> 


#undef VERTEX_TYPE 
#undef __VCGLIB_VERTEX_VC
#undef __VCGLIB_VERTEX_VN
#undef __VCGLIB_VERTEX_VQ


namespace vcg {
typedef VertexVCVNVQ<float>  VertexVCVNVQf;
typedef VertexVCVNVQ<double> VertexVCVNVQd;
}

#endif
