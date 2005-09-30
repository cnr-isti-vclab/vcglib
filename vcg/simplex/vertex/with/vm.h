#ifndef __VCGLIB_VERTEX__VM__TYPE
#define __VCGLIB_VERTEX__VM__TYPE


#define VERTEX_TYPE VertexVM

#define __VCGLIB_VERTEX_VM

#include <vcg/simplex/vertex/base.h> 


#undef VERTEX_TYPE 
#undef __VCGLIB_VERTEX_VM



namespace vcg {

template < class VETYPE, class VFTYPE, class VTTYPE>
class VertexVMf : public VertexVM<float,VETYPE,VFTYPE,VTTYPE> {};

template < class VETYPE, class VFTYPE, class VTTYPE>
class VertexVMd : public VertexVM<double,VETYPE,VFTYPE,VTTYPE> {};

}

#endif