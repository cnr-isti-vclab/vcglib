#ifndef __VCGLIB_VERTEX__AFVN__TYPE
#define __VCGLIB_VERTEX__AFVN__TYPE


#define VERTEX_TYPE VertexAFVN 

#define __VCGLIB_VERTEX_AF
#define __VCGLIB_VERTEX_VN


#include <vcg/simplex/vertex/base.h> 


#undef VERTEX_TYPE 
#undef __VCGLIB_VERTEX_AF
#undef __VCGLIB_VERTEX_VN


namespace vcg {

template < class VETYPE, class VFTYPE, class VTTYPE>
class VertexAFVNf : public VertexAFVN<float,VETYPE,VFTYPE,VTTYPE> {};

template < class VETYPE, class VFTYPE, class VTTYPE>
class VertexAFVNd : public VertexAFVN<double,VETYPE,VFTYPE,VTTYPE> {};

}

#endif