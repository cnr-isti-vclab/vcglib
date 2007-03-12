#ifndef __VCGLIB_VERTEX__BASE__TYPE
#define __VCGLIB_VERTEX__BASE__TYPE


#define VERTEX_TYPE Vertex

#include <vcg/simplex/vertex/base.h> 


#undef VERTEX_TYPE 


namespace vcg {
template < class VETYPE, class VFTYPE, class VTTYPE,class TCTYPE = TexCoord2<float,1>, class CoordTYPE= Point3<float> >
class Vertexf : public Vertex<float,VETYPE,VFTYPE,VTTYPE, TCTYPE , CoordTYPE> {};

template < class VETYPE, class VFTYPE, class VTTYPE>
class Vertexd : public Vertex<double,VETYPE,VFTYPE,VTTYPE> {};

}

#endif
