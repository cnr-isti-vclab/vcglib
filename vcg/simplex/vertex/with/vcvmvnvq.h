#ifndef __VCGLIB_VERTEX__VCVMVNVQ__TYPE
#define __VCGLIB_VERTEX__VCVMVNVQ__TYPE


#define VERTEX_TYPE VertexVCVMVNVQ 

#define __VCGLIB_VERTEX_VC
#define __VCGLIB_VERTEX_VM
#define __VCGLIB_VERTEX_VN
#define __VCGLIB_VERTEX_VQ

#include <vcg/simplex/vertex/base.h> 


#undef VERTEX_TYPE 
#undef __VCGLIB_VERTEX_VC
#undef __VCGLIB_VERTEX_VM
#undef __VCGLIB_VERTEX_VN
#undef __VCGLIB_VERTEX_VQ


namespace vcg {
typedef  VertexVCVMVNVQ<float> VertexVCVMVNVQf;
typedef  VertexVCVMVNVQ<double> VertexVCVMVNVQd;
}

#endif
