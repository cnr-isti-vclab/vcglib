#ifndef __VCGLIB_VERTEX__AF__TYPE
#define __VCGLIB_VERTEX__VM__TYPE
#define __VCGLIB_VERTEX__VN__TYPE
#define __VCGLIB_VERTEX__VQ__TYPE

#define VERTEX_TYPE VertexAFVMVNVQ

#define __VCGLIB_VERTEX_AF
#define __VCGLIB_VERTEX_VM
#define __VCGLIB_VERTEX_VN
#define __VCGLIB_VERTEX_VQ

#include <vcg/simplex/vertex/base.h> 

#undef VERTEX_TYPE 
#undef __VCGLIB_VERTEX__AF__TYPE
#undef __VCGLIB_VERTEX__VM__TYPE
#undef __VCGLIB_VERTEX__VN__TYPE
#undef __VCGLIB_VERTEX__VQ__TYPE

namespace vcg {

template < class VETYPE, class VFTYPE, class VTTYPE>
class VertexAFVMVNVQf : public VertexAFVMVNVQ<float,VETYPE,VFTYPE,VTTYPE> {};

template < class VETYPE, class VFTYPE, class VTTYPE>
class VertexAFVMVNVQd : public VertexAFVMVNVQ<double,VETYPE,VFTYPE,VTTYPE> {};

}

#endif
