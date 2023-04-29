#ifndef MY_TRI_MESH_TYPE
#define MY_TRI_MESH_TYPE

#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/update/bounding.h>
#include <vcg/complex/algorithms/smooth.h>
#include <wrap/io_trimesh/import.h>
#include <wrap/gl/trimesh.h>

class MyTriFace;
class MyTriEdge;
class MyTriVertex;

struct TriUsedTypes : public vcg::UsedTypes<vcg::Use<MyTriVertex>::AsVertexType,
                                            vcg::Use<MyTriEdge>::AsEdgeType,
                                            vcg::Use<MyTriFace>::AsFaceType>
{
};

class MyTriVertex : public vcg::Vertex<TriUsedTypes,
                                       vcg::vertex::Coord3f,
                                       vcg::vertex::Normal3f,
                                       vcg::vertex::BitFlags>
{
};

class MyTriEdge : public vcg::Edge<
                      TriUsedTypes,
                      vcg::edge::VertexRef,
                      vcg::edge::BitFlags>
{
};

class MyTriFace : public vcg::Face<TriUsedTypes,
                                   vcg::face::VertexRef,
                                   vcg::face::BitFlags,
                                   vcg::face::Normal3f>
{
};

class MyTriMesh : public vcg::tri::TriMesh<std::vector<MyTriVertex>,
                                           std::vector<MyTriEdge>,
                                           std::vector<MyTriFace>>
{
};


#endif
