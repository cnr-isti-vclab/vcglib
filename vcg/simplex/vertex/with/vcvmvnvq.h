#ifndef __VCGLIB_VERTEX__VCVN__TYPE
#define __VCGLIB_VERTEX__VCVN__TYPE


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

template < class VETYPE, class VFTYPE, class VTTYPE>
class VertexVCVMVNVQf : public VertexVCVMVNVQ<float,VETYPE,VFTYPE,VTTYPE> {};

template < class VETYPE, class VFTYPE, class VTTYPE>
class VertexVCVMVNVQd : public VertexVCVMVNVQ<double,VETYPE,VFTYPE,VTTYPE> {};

}

#endif
