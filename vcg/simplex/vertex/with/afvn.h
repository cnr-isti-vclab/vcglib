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

template <class VFTYPE> 
class VertexAFVNf : public VertexAFVN<float,VFTYPE> {};

template <class VFTYPE> 
class VertexAFVNd : public VertexAFVN<double,VFTYPE> {};

}

#endif