#include <vector>

#include <vcg/simplex/vertexplus/base.h>   
#include <vcg/simplex/vertexplus/component.h>   
#include <vcg/simplex/faceplus/base.h>   
#include <vcg/simplex/faceplus/component.h>   

#include <vcg/complex/trimesh/base.h>   

class MyEdge;
class MyFace;

class MyVertex: public vcg::VertexSimp2<MyVertex,MyEdge,MyFace, vcg::vert::Coord3d, vcg::vert::Normal3f>{};
class MyFace: public vcg::FaceSimp2<MyVertex,MyEdge,MyFace, vcg::face::VertexRef>{};

class MyMesh: public vcg::tri::TriMesh< std::vector<MyVertex>, std::vector<MyFace> > {};

int main()
{
	MyMesh m;
	return 0;
}