#ifndef __VCGLIB_VERTEX__VCVN__TYPE
#define __VCGLIB_VERTEX__VCVN__TYPE


#define VERTEX_TYPE VertexVC 

#define __VCGLIB_VERTEX_VC

#include <vcg/simplex/vertex/base.h> 


#undef VERTEX_TYPE 
#undef __VCGLIB_VERTEX_VC


namespace vcg {

template < class VETYPE, class VFTYPE, class VTTYPE>
class VertexVCf : public VertexVC<float,VETYPE,VFTYPE,VTTYPE> {};

template < class VETYPE, class VFTYPE, class VTTYPE>
class VertexVCd : public VertexVC<double,VETYPE,VFTYPE,VTTYPE> {};

}

#endif
