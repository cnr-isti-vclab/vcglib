#ifndef __VCGLIB_VERTEX__EA__TYPE
#define __VCGLIB_VERTEX__EA__TYPE

#define VERTEX_TYPE VertexEA 

#define __VCGLIB_VERTEX_EA 

#include <vcg/simplex/vertex/base.h> 

#undef VERTEX_TYPE 
#undef __VCGLIB_VERTEX_EA 

namespace vcg {

template <class VFTYPE> 
class VertexEAf : public VertexEA<float,VFTYPE> {};

template <class VFTYPE> 
class VertexEAd : public VertexEA<double,VFTYPE> {};

}

#endif
