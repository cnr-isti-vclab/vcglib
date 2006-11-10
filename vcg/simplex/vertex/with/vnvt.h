#ifndef __VCGLIB_VERTEX__VN__TYPE
#define __VCGLIB_VERTEX__VN__TYPE


#define VERTEX_TYPE VertexVNVT 

#define __VCGLIB_VERTEX_VNVT
#define __VCGLIB_VERTEX_VTVT

#include <vcg/simplex/vertex/base.h> 


#undef VERTEX_TYPE 
#undef __VCGLIB_VERTEX_VN
#undef __VCGLIB_VERTEX_VT


namespace vcg {
template < class VETYPE, class VFTYPE, class VTTYPE>
class VertexVNVTf : public VertexVNVT<float,VETYPE,VFTYPE,VTTYPE> {};

template < class VETYPE, class VFTYPE, class VTTYPE>
class VertexVNVTd : public VertexVNVT<double,VETYPE,VFTYPE,VTTYPE> {};

}

#endif
