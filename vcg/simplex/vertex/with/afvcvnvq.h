#ifndef __VCGLIB_VERTEX__AFVCVNVQ__TYPE
#define __VCGLIB_VERTEX__AFVCVNVQ__TYPE


#define VERTEX_TYPE VertexAFVCVNVQ 

#define __VCGLIB_VERTEX_AF
#define __VCGLIB_VERTEX_VC
#define __VCGLIB_VERTEX_VN
#define __VCGLIB_VERTEX_VQ


#include <vcg/simplex/vertex/base.h> 


#undef VERTEX_TYPE 
#undef __VCGLIB_VERTEX_AF
#undef __VCGLIB_VERTEX_VC
#undef __VCGLIB_VERTEX_VN
#undef __VCGLIB_VERTEX_VQ


namespace vcg {

template < class VETYPE, class VFTYPE, class VTTYPE>
class VertexAFVCVNVQf : public VertexAFVCVNVQ<float,VETYPE,VFTYPE,VTTYPE> {};

template < class VETYPE, class VFTYPE, class VTTYPE>
class VertexAFVCVNVQd : public VertexAFVCVNVQ<double,VETYPE,VFTYPE,VTTYPE> {};

}

#endif
