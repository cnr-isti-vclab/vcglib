#ifndef __VCGLIB_VERTEX__AT__TYPE
#define __VCGLIB_VERTEX__AT__TYPE

#define VERTEX_TYPE VertexATVN

#define __VCGLIB_VERTEX_AT
#define __VCGLIB_VERTEX_VN
#define __VCGLIB_VERTEX_VM

#include <vcg/simplex/vertex/base.h> 

#undef VERTEX_TYPE 
#undef __VCGLIB_VERTEX_VM
#undef __VCGLIB_VERTEX_AT
#undef __VCGLIB_VERTEX_VN

namespace vcg {

template < class VETYPE, class VFTYPE, class VTTYPE>
class VertexATVNf : public VertexATVN<float,VETYPE,VFTYPE,VTTYPE> {};

template < class VETYPE, class VFTYPE, class VTTYPE>
class VertexATVNd : public VertexATVN<double,VETYPE,VFTYPE,VTTYPE> {};

}

#endif
