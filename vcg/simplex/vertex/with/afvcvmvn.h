#ifndef __VCGLIB_VERTEX__VCVN__TYPE
#define __VCGLIB_VERTEX__VCVN__TYPE


#define VERTEX_TYPE VertexAFVCVMVN 

#define __VCGLIB_VERTEX_AF
#define __VCGLIB_VERTEX_VN
#define __VCGLIB_VERTEX_VC
#define __VCGLIB_VERTEX_VM

#include <vcg/simplex/vertex/base.h> 


#undef VERTEX_TYPE 
#undef __VCGLIB_VERTEX_VM
#undef __VCGLIB_VERTEX_AF
#undef __VCGLIB_VERTEX_VN
#undef __VCGLIB_VERTEX_VC


namespace vcg {

template < class VETYPE, class VFTYPE, class VTTYPE>
class VertexAFVCVMVNf : public VertexAFVCVMVN<float,VETYPE,VFTYPE,VTTYPE> {};

template < class VETYPE, class VFTYPE, class VTTYPE>
class VertexAFVCVMVNd : public VertexAFVCVMVN<double,VETYPE,VFTYPE,VTTYPE> {};

}

#endif
