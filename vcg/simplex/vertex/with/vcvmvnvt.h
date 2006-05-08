#ifndef __VCGLIB_VERTEX__VCVMVNVT__TYPE
#define __VCGLIB_VERTEX__VCVMVNVT__TYPE


#define VERTEX_TYPE VertexVCVMVNVT

#define __VCGLIB_VERTEX_VC
#define __VCGLIB_VERTEX_VM
#define __VCGLIB_VERTEX_VN
#define __VCGLIB_VERTEX_VT

#include <vcg/simplex/vertex/base.h> 


#undef VERTEX_TYPE 
#undef __VCGLIB_VERTEX_VC
#undef __VCGLIB_VERTEX_VM
#undef __VCGLIB_VERTEX_VN
#undef __VCGLIB_VERTEX_VT


namespace vcg {

template < class VETYPE, class VFTYPE, class VTTYPE>
class VertexVCVMVNVTf : public VertexVCVMVNVT<float,VETYPE,VFTYPE,VTTYPE> {};

template < class VETYPE, class VFTYPE, class VTTYPE>
class VertexVCVMVNVTd : public VertexVCVMVNVT<double,VETYPE,VFTYPE,VTTYPE> {};

}

#endif /* __VCGLIB_VERTEX__VCVMVNVT__TYPE */
