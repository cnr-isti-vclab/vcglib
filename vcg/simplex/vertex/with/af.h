#ifndef __VCGLIB_VERTEX__AF__TYPE
#define __VCGLIB_VERTEX__AF__TYPE

#define VERTEX_TYPE VertexAF

#define __VCGLIB_VERTEX_AF

#include <vcg/simplex/vertex/base.h> 

#undef VERTEX_TYPE 
#undef __VCGLIB_VERTEX_AF

namespace vcg {

template < class VETYPE, class VFTYPE, class VTTYPE>
class VertexAFf : public VertexAF<float,VETYPE,VFTYPE,VTTYPE> {};

template < class VETYPE, class VFTYPE, class VTTYPE>
class VertexAFd : public VertexAF<double,VETYPE,VFTYPE,VTTYPE> {};

}

#endif
