#include <vector>

#include <vcg/simplex/vertexplus/base.h>   
#include <vcg/simplex/vertexplus/component.h>   
#include <vcg/simplex/faceplus/base.h>   
#include <vcg/simplex/faceplus/component.h>   

#include <vcg/complex/trimesh/base.h>   
#include<vcg/complex/trimesh/create/platonic.h>

#include <vcg/simplex/face/pos.h> 

class MyEdge;
class MyFace;

class MyVertex: public vcg::VertexSimp2<MyVertex,MyEdge,MyFace, vcg::vert::Coord3d, vcg::vert::Normal3f>{};
class MyFace: public vcg::FaceSimp2<MyVertex,MyEdge,MyFace, vcg::face::VertexRef,vcg::face::FFAdj>{};

class MyMesh: public vcg::tri::TriMesh< std::vector<MyVertex>, std::vector<MyFace> > {};

void OneRingNeighborhood( MyFace * f)
{
	MyVertex * v = f->V(0);  
	MyFace* start = f;
	vcg::face::Pos<MyFace> p(f,0,0);// constructor that takes face, edge and vertex
	do
	{
		p.FlipF();
		p.FlipE();
	}while(p.f!=start);
}

int main()
{
	MyMesh m;
	vcg::tri::Tetrahedron(m);
	OneRingNeighborhood(&(*m.face.begin()));
	return 0;
}