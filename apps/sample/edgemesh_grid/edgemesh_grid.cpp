#include <stdio.h>
#include <vcg/space/color4.h>
#include <vcg/space/index/grid_static_ptr.h>
#include <vcg/simplex/vertex/base.h>
#include <vcg/simplex/vertex/component.h>
#include <vcg/simplex/edge/base.h>
#include <vcg/simplex/edge/component.h>
#include <vcg/complex/edgemesh/base.h>
#include <vcg/complex/edgemesh/allocate.h>
#include <vcg/complex/edgemesh/update/bounding.h>
#include <vcg/complex/edgemesh/closest.h>
#include <vcg/complex/trimesh/closest.h>
#include <vcg/complex/used_types.h>

//
//using namespace std;
//
class MyFace;
class MyEdge;
class MyVertex;
struct MyUsedTypes : public vcg::UsedTypes<	vcg::Use<MyVertex>		::AsVertexType,
											vcg::Use<MyEdge>		::AsEdgeType,
											vcg::Use<MyFace>		::AsFaceType>{};

class MyVertex  : public vcg::Vertex<MyUsedTypes,vcg::vertex::BitFlags,vcg::vertex::Coord3f> {};
class MyEdge    : public vcg::Edge<MyUsedTypes,vcg::edge::Mark,vcg::edge::VertexRef,vcg::edge::BitFlags> {};
class MyEdgeMesh: public vcg::edg::EdgeMesh< std::vector<MyVertex>, std::vector<MyEdge> > {};

typedef vcg::GridStaticPtr<MyEdge, MyEdge::ScalarType> EdgeMeshGrid;

#define VERT_NUMB 100

int main(int , char **)
{
  //create a random mesh of edge
  MyEdgeMesh em;
  srand(1000);
  em.vert.reserve(VERT_NUMB);
  for (int i=0;i<VERT_NUMB;i=i+2)
  {
	float x0=((float)rand()/(float)RAND_MAX)*1000.f;
	float y0=((float)rand()/(float)RAND_MAX)*1000.f;
	float z0=((float)rand()/(float)RAND_MAX)*1000.f;
	float x1=((float)rand()/(float)RAND_MAX)*1000.f;
	float y1=((float)rand()/(float)RAND_MAX)*1000.f;
	float z1=((float)rand()/(float)RAND_MAX)*1000.f;
	em.vert.push_back(MyVertex());
	MyVertex *v0=&em.vert.back();
	em.vert.push_back(MyVertex());
	MyVertex *v1=&em.vert.back();
	v0->P().X()=x0;
	v0->P().Y()=y0;
	v0->P().Z()=z0;
	v1->P().X()=x1;
	v1->P().Y()=y1;
	v1->P().Z()=z1;
	em.edges.push_back(MyEdge());
	MyEdge *e=&em.edges.back();
	e->V(0)=v0;
	e->V(1)=v1;
  }
  vcg::edg::UpdateBounding<MyEdgeMesh>::Box(em);
  EdgeMeshGrid static_grid;
  static_grid.Set(em.edges.begin(), em.edges.end());
  float dist;
  vcg::Point3f p;
  MyEdge *e=vcg::edgemesh::GetClosestEdge<MyEdgeMesh,EdgeMeshGrid>(em,static_grid,vcg::Point3f(500,500,500),1000,dist,p);
  std::vector<MyEdge*> ret;
  int num=vcg::edgemesh::GetInBoxEdge<MyEdgeMesh,EdgeMeshGrid,std::vector<MyEdge*> >(em,static_grid,em.bbox,ret);
  return 0;
}
