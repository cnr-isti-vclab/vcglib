#ifndef __VCGLIB_VERTEX__AF__TYPE
#define __VCGLIB_VERTEX__VM__TYPE
#define __VCGLIB_VERTEX__VN__TYPE

#define VERTEX_TYPE VertexAFVMVN

#define __VCGLIB_VERTEX_AF
#define __VCGLIB_VERTEX_VM
#define __VCGLIB_VERTEX_VN

#include <vcg/simplex/vertex/base.h> 

#undef VERTEX_TYPE 
#undef __VCGLIB_VERTEX__AF__TYPE
#undef __VCGLIB_VERTEX__VM__TYPE
#undef __VCGLIB_VERTEX__VN__TYPE

namespace vcg {

template < class VETYPE, class VFTYPE, class VTTYPE>
class VertexAFVMVNf : public VertexAFVMVN<float,VETYPE,VFTYPE,VTTYPE> {};

template < class VETYPE, class VFTYPE, class VTTYPE>
class VertexAFVMVNd : public VertexAFVMVN<double,VETYPE,VFTYPE,VTTYPE> {};

}

#endif
