#ifndef __VCGLIB_VERTEX__VN__TYPE
#define __VCGLIB_VERTEX__VN__TYPE


#define VERTEX_TYPE VertexVN 

#define __VCGLIB_VERTEX_VN

#include <vcg/simplex/vertex/base.h> 


#undef VERTEX_TYPE 
#undef __VCGLIB_VERTEX_VN


namespace vcg {
template < class VETYPE, class VFTYPE, class VTTYPE>
class VertexVNf : public VertexVN<float,VETYPE,VFTYPE,VTTYPE> {};

template < class VETYPE, class VFTYPE, class VTTYPE>
class VertexVNd : public VertexVN<double,VETYPE,VFTYPE,VTTYPE> {};

}

#endif
