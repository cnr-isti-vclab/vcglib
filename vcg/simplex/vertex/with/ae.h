#ifndef __VCGLIB_VERTEX__AE__TYPE
#define __VCGLIB_VERTEX__AE__TYPE

#define VERTEX_TYPE VertexAE 

#define __VCGLIB_VERTEX_AE

#include <vcg/simplex/vertex/base.h> 

#undef VERTEX_TYPE 
#undef __VCGLIB_VERTEX_AE

namespace vcg {

template <class VFTYPE> 
class VertexAEf : public VertexAE<float,VFTYPE> {};

template <class VFTYPE> 
class VertexAEd : public VertexAE<double,VFTYPE> {};

}

#endif
