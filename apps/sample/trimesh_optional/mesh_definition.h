#ifndef _MESH_DEF_
#define _MESH_DEF_

#include<vcg/complex/complex.h>
#include<vcg/simplex/vertex/base.h>
#include<vcg/simplex/face/base.h>
#include<vcg/simplex/face/topology.h>

#include<vcg/simplex/vertex/component_ocf.h>
#include<vcg/simplex/face/component_ocf.h>

 

class CFace;
class CFaceOcf;
 
class CVertex;
class CVertexOcf;
 

struct MyUsedTypes:		 public	vcg::UsedTypes<vcg::Use<CVertex>::AsVertexType,vcg::Use<CFace>::AsFaceType>{};
struct MyUsedTypesOcf: public vcg::UsedTypes<vcg::Use<CVertexOcf>::AsVertexType,vcg::Use<CFaceOcf>::AsFaceType>{};

// Optional stuff has two suffixes:
// OCF Optional Component Fast 
// OCC Optional Component Compact

class CVertex     : public vcg::Vertex<	MyUsedTypes,  vcg::vertex::Coord3f, vcg::vertex::BitFlags,vcg::vertex::Normal3f >{};
class CVertexOcf  : public vcg::Vertex< MyUsedTypesOcf,vcg::vertex::InfoOcf,vcg::vertex::Coord3f,vcg::vertex::QualityfOcf, vcg::vertex::BitFlags,vcg::vertex::Normal3f,vcg::vertex::RadiusfOcf >{};
class CFace       : public vcg::Face< MyUsedTypes,    vcg::face::FFAdj,    vcg::face::VertexRef, vcg::face::BitFlags, vcg::face::Normal3f > {};
class CFaceOcf    : public vcg::Face< MyUsedTypesOcf, vcg::face::InfoOcf, vcg::face::FFAdjOcf, vcg::face::VertexRef, vcg::face::BitFlags, vcg::face::Normal3fOcf > {};

class CMesh       : public vcg::tri::TriMesh<     std::vector<CVertex   >,           std::vector<CFace   > > {};
class CMeshOcf    : public vcg::tri::TriMesh<     vcg::vertex::vector_ocf<CVertexOcf>, vcg::face::vector_ocf<CFaceOcf> > {};

#endif
