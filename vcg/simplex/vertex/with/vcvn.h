#ifndef __VCGLIB_VERTEX__VCVN__TYPE
#define __VCGLIB_VERTEX__VCVN__TYPE


#define VERTEX_TYPE VertexVCVN 

#define __VCGLIB_VERTEX_VN
#define __VCGLIB_VERTEX_VC

#include <vcg/simplex/vertex/base.h> 


#undef VERTEX_TYPE 
#undef __VCGLIB_VERTEX_VN
#undef __VCGLIB_VERTEX_VC


namespace vcg {
typedef VertexVCVN<float>  VertexVCVNf;
typedef VertexVCVN<double> VertexVCVNd;
}

#endif
