#ifndef __VCGLIB_VERTEX__A__TYPE
#define __VCGLIB_VERTEX__A__TYPE

#define VERTEX_TYPE VertexA 

#define __VCGLIB_VERTEX_A 

#include <vcg/simplex/vertex/base.h> 

#undef VERTEX_TYPE 
#undef __VCGLIB_VERTEX_A 

namespace vcg {

template <class VFTYPE> 
class VertexAf : public VertexA<float,VFTYPE> {};

template <class VFTYPE> 
class VertexAd : public VertexA<double,VFTYPE> {};

}

#endif
