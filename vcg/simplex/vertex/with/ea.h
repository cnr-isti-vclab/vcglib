#ifndef __VCGLIB_VERTEX__EA__TYPE
#define __VCGLIB_VERTEX__EA__TYPE

#define VERTEX_TYPE VertexEA 

#define __VCGLIB_VERTEX_EA 

#include <vcg/simplex/vertex/base.h> 

#undef VERTEX_TYPE 
#undef __VCGLIB_VERTEX_EA 

namespace vcg {

template < class VETYPE, class VFTYPE, class VTTYPE>
class VertexEAf : public VertexEA<float,VETYPE,VFTYPE,VTTYPE> {};

template < class VETYPE, class VFTYPE, class VTTYPE>
class VertexEAd : public VertexEA<double,VETYPE,VFTYPE,VTTYPE> {};

}

#endif
