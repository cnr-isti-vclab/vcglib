#include <vector>

#include<vcg/simplex/vertex/vertex.h>

#include<vcg/simplex/face/with/af.h>

#include<vcg/complex/trimesh/base.h>

// loader 
#include<wrap/io_trimesh/import_ply.h>

// topology computation
#include<vcg/complex/trimesh/update/topology.h>

// half edge iterators
#include<vcg/simplex/face/pos.h>



using namespace vcg;

class MyFace;
class MyEdge;
class MyVertex:public Vertex<float,MyEdge,MyFace>{};
class MyFace :public FaceAF<MyVertex,MyEdge,MyFace>{};
class MyMesh: public tri::TriMesh< std::vector<MyVertex>, std::vector<MyFace > >{};


int main(int argc,char ** argv){

	MyMesh m;

	//load the mesh
	vcg::tri::io::ImporterPLY<MyMesh>::Open(m,argv[1]);

	//update the face-face topology 
	vcg::tri::UpdateTopology<MyMesh>::FaceFace(m);

	// declare an iterator on the mesh
	vcg::face::Pos<MyMesh::FaceType> he;

	// set as the half edge: (first face of the mesh, edge 0,vertex 0)
	he.Set(&*m.face.begin(),0,(*m.face.begin()).V(0));

	// same vertex and edge..adjacent face
	he.FlipF();
	//again (it went back..it it manifold)	
	he.FlipF(); 

	}

