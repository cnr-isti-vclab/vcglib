#ifndef __VCGLIB_VERTEX__VMVNVT__TYPE
#define __VCGLIB_VERTEX__VMVNVT__TYPE


#define VERTEX_TYPE VertexVMVNVT

#define __VCGLIB_VERTEX_VM
#define __VCGLIB_VERTEX_VN
#define __VCGLIB_VERTEX_VT

#include <vcg/simplex/vertex/base.h> 

#undef __VCGLIB_VERTEX_VM
#undef __VCGLIB_VERTEX_VN
#undef __VCGLIB_VERTEX_VT

#undef VERTEX_TYPE 

namespace vcg {

template < class VETYPE, class VFTYPE, class VTTYPE>
class VertexVMVNVTf : public VertexVMVNVT<float,VETYPE,VFTYPE,VTTYPE> {};

template < class VETYPE, class VFTYPE, class VTTYPE>
class VertexVMVNVTd : public VertexVMVNVT<double,VETYPE,VFTYPE,VTTYPE> {};

}


#endif /* __VCGLIB_VERTEX__VMVNVT__TYPE */
