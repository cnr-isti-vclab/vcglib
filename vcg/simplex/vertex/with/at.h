#ifndef __VCGLIB_VERTEX__AT__TYPE
#define __VCGLIB_VERTEX__AT__TYPE

#define VERTEX_TYPE VertexAT

#define __VCGLIB_VERTEX_AT

#include <vcg/simplex/vertex/base.h> 

#undef VERTEX_TYPE 
#undef __VCGLIB_VERTEX_AT

namespace vcg {

template <class VTTYPE> 
class VertexATf : public VertexAT<float,VTTYPE> {};

template <class VFTYPE> 
class VertexATd : public VertexAT<double,VTTYPE> {};

}

#endif
