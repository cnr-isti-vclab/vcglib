#ifndef __VCGLIB_VERTEX__AE__TYPE
#define __VCGLIB_VERTEX__AE__TYPE

#define VERTEX_TYPE VertexAE 

#define __VCGLIB_VERTEX_AE

#include <vcg/simplex/vertex/base.h> 

#undef VERTEX_TYPE 
#undef __VCGLIB_VERTEX_AE

namespace vcg {

template < class VETYPE, class VFTYPE, class VTTYPE>
class VertexAEf : public VertexAE<float,VETYPE,VFTYPE,VTTYPE> {};

template < class VETYPE, class VFTYPE, class VTTYPE>
class VertexAEd : public VertexAE<double,VETYPE,VFTYPE,VTTYPE> {};

}

#endif
