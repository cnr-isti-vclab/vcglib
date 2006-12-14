#ifndef __VCGLIB_VERTEX__VMVN__TYPE
#define __VCGLIB_VERTEX__VMVN__TYPE


#define VERTEX_TYPE VertexVMVN

#define __VCGLIB_VERTEX_VN
#define __VCGLIB_VERTEX_VM

#include <vcg/simplex/vertex/base.h> 


#undef VERTEX_TYPE 
#undef __VCGLIB_VERTEX_VM
#undef __VCGLIB_VERTEX_VN


namespace vcg {

template < class VETYPE, class VFTYPE, class VTTYPE>
class VertexVMVNf : public VertexVMVN<float,VETYPE,VFTYPE,VTTYPE> {};

template < class VETYPE, class VFTYPE, class VTTYPE>
class VertexVMVNd : public VertexVMVN<double,VETYPE,VFTYPE,VTTYPE> {};

}

#endif
