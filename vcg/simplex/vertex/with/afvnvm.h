#ifndef __VCGLIB_VERTEX__AFVNVM__TYPE
#define __VCGLIB_VERTEX__AFVNVM__TYPE


#define VERTEX_TYPE VertexAFVNVM 

#define __VCGLIB_VERTEX_AF
#define __VCGLIB_VERTEX_VM
#define __VCGLIB_VERTEX_VN


#include <vcg/simplex/vertex/base.h> 


#undef VERTEX_TYPE 
#undef __VCGLIB_VERTEX_AF
#undef __VCGLIB_VERTEX_VM
#undef __VCGLIB_VERTEX_VN


namespace vcg {

template < class VETYPE, class VFTYPE, class VTTYPE>
class VertexAFVNVMf : public VertexAFVNVM<float,VETYPE,VFTYPE,VTTYPE> {};

template < class VETYPE, class VFTYPE, class VTTYPE>
class VertexAFVNVMd : public VertexAFVNVM<double,VETYPE,VFTYPE,VTTYPE> {};

}

#endif