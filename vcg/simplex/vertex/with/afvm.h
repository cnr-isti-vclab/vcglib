#ifndef __VCGLIB_VERTEX__AFVM__TYPE
#define __VCGLIB_VERTEX__AFVM__TYPE


#define VERTEX_TYPE VertexAFVM 

#define __VCGLIB_VERTEX_AF
#define __VCGLIB_VERTEX_VM


#include <vcg/simplex/vertex/base.h> 


#undef VERTEX_TYPE 
#undef __VCGLIB_VERTEX_AF
#undef __VCGLIB_VERTEX_VM


namespace vcg {

template < class VETYPE, class VFTYPE, class VTTYPE>
class VertexAFVMf : public VertexAFVM<float,VETYPE,VFTYPE,VTTYPE> {};

template < class VETYPE, class VFTYPE, class VTTYPE>
class VertexAFVMd : public VertexAFVM<double,VETYPE,VFTYPE,VTTYPE> {};

}

#endif