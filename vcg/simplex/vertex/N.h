/*#**************************************************************************
  History

 2000	Jan 31 First Working release
						 
****************************************************************************/
#ifndef __VCGLIB_VERTEX__N__TYPE
#define __VCGLIB_VERTEX__N__TYPE

#define VERTEX_TYPE VertexN 

#define __VCGLIB_VERTEX_N

#include <vcg/simplex/vertex/base.h> 

#undef VERTEX_TYPE 

#undef __VCGLIB_VERTEX_N

using namespace vcg;

typedef VertexN<short>  VertexNs;
typedef VertexN<int>	  VertexNi;
typedef VertexN<float>  VertexNf;
typedef VertexN<double> VertexNd;

#endif