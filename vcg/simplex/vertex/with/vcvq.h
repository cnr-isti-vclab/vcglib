#ifndef __VCGLIB_VERTEX__VCVQ__TYPE
#define __VCGLIB_VERTEX__VCVQ__TYPE


#define VERTEX_TYPE VertexVCVQ 

#define __VCGLIB_VERTEX_VC
#define __VCGLIB_VERTEX_VQ

#include <vcg/simplex/vertex/base.h> 


#undef VERTEX_TYPE 

#undef __VCGLIB_VERTEX_VC
#undef __VCGLIB_VERTEX_VQ

namespace vcg {
typedef VertexVCVQ<float>  VertexVCVQf;
typedef VertexVCVQ<double> VertexVCVQd;
}

#endif
