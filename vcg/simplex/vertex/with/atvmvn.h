#ifndef __VCGLIB_VERTEX__ATVMVN__TYPE
#define __VCGLIB_VERTEX__ATVMVN__TYPE

#define VERTEX_TYPE VertexATVMVN

#define __VCGLIB_VERTEX_AT
#define __VCGLIB_VERTEX_VM
#define __VCGLIB_VERTEX_VN

#include <vcg/simplex/vertex/base.h> 

#undef VERTEX_TYPE 
#undef __VCGLIB_VERTEX_AT
#undef __VCGLIB_VERTEX_VM
#undef __VCGLIB_VERTEX_VN

namespace vcg {

template < class VETYPE, class VFTYPE, class VTTYPE,class TCTYPE = TexCoord2<float,1>, class CoordTYPE= Point3<float> >
class VertexATVMVNf : public VertexATVMVN<float,VETYPE,VFTYPE,VTTYPE,TCTYPE,CoordTYPE> {};

template < class VETYPE, class VFTYPE, class VTTYPE>
class VertexATVMVNd : public VertexATVMVN<double,VETYPE,VFTYPE,VTTYPE> {};

}

#endif
