#ifndef __VCGLIB_VERTEX__AT__TYPE
#define __VCGLIB_VERTEX__AT__TYPE

#define VERTEX_TYPE VertexATVN

#define __VCGLIB_VERTEX_AT
#define __VCGLIB_VERTEX_VN

#include <vcg/simplex/vertex/base.h> 

#undef VERTEX_TYPE 
#undef __VCGLIB_VERTEX_AT
#undef __VCGLIB_VERTEX_VN

namespace vcg {

template <class VTTYPE> 
class VertexATVNf : public VertexATVN<float,DUMMYFACETYPE,VTTYPE> {};

template <class VTTYPE> 
class VertexATVNd : public VertexATVN<double,DUMMYFACETYPE,VTTYPE> {};

}

#endif
