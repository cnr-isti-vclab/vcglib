
#ifndef MESH_TYPE_H
#define MESH_TYPE_H

///vcg imports
#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/update/bounding.h>
#include <vcg/complex/algorithms/update/normal.h>
#include <vcg/complex/algorithms/update/topology.h>
#include <vcg/simplex/face/topology.h>
#include <vcg/complex/algorithms/parametrization/tangent_field_operators.h>
///wrapper imports
#include <wrap/io_trimesh/import.h>
#include <wrap/io_trimesh/export.h>
#include <wrap/miq/core/param_stats.h>

using namespace vcg;
class CFace;
class CVertex;

struct MyUsedTypes : public UsedTypes<	Use<CVertex>::AsVertexType, Use<CFace>::AsFaceType >{};

/// compositing wanted proprieties
class CVertex : public vcg::Vertex< MyUsedTypes,
        vcg::vertex::Coord3d, vcg::vertex::Normal3d,
        vcg::vertex::BitFlags,vcg::vertex::VFAdj,
        vcg::vertex::TexCoord2d,vcg::vertex::Qualityd>{};

class CFace   : public vcg::Face<  MyUsedTypes, vcg::face::VertexRef,
        vcg::face::VFAdj, vcg::face::FFAdj,vcg::face::Normal3d,
        vcg::face::WedgeTexCoord2d,vcg::face::BitFlags ,
        vcg::face::CurvatureDird,vcg::face::Qualityd,vcg::face::Color4b,
        vcg::face::Mark>{};


class CMesh   : public vcg::tri::TriMesh< std::vector<CVertex>, std::vector<CFace> >{};


class MyPolyFace;
class MyPolyVertex;
struct PolyUsedTypes: public vcg::UsedTypes<vcg::Use<MyPolyVertex>	::AsVertexType,
                                            vcg::Use<MyPolyFace>	::AsFaceType
                                            >{};

//class DummyEdge: public vcg::Edge<PolyUsedTypes>{};
class MyPolyVertex:public vcg::Vertex<	PolyUsedTypes,
                                        vcg::vertex::Coord3f,
                                        vcg::vertex::Normal3f,
                                        vcg::vertex::BitFlags>{} ;

class MyPolyFace:public vcg::Face<
     PolyUsedTypes
    ,vcg::face::PolyInfo // this is necessary  if you use component in vcg/simplex/face/component_polygon.h
                         // It says "this class is a polygon and the memory for its components (e.g. pointer to its vertices
                         // will be allocated dynamically")
    ,vcg::face::PFVAdj	 // Pointer to the vertices (just like FVAdj )
    ,vcg::face::BitFlags // bit flags
    ,vcg::face::Normal3f // normal
> {};

class MyPolyMesh: public
    vcg::tri::TriMesh<
    std::vector<MyPolyVertex>,	// the vector of vertices
    std::vector<MyPolyFace > 						// the vector of faces
    >{};


#endif
