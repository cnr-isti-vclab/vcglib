

#include <vector>

// stuff to define the mesh
#include <vcg/simplex/vertex/with/atvmvn.h>
#include <vcg/complex/tetramesh/base.h>
#include <vcg/simplex/tetrahedron/with/atavtq.h>

// the trackball
#include <wrap/gui/trackball.h>

// io
#include <wrap/io_tetramesh/import_ply.h>
#include <wrap/io_tetramesh/export_ply.h>
#include <wrap/io_tetramesh/import_ts.h>


class MyEdge;
class MyTetrahedron;
class MyFace;
class MyVertex:public vcg::VertexATVMVNf<DUMMYEDGETYPE , MyFace, MyTetrahedron>{} ;
class MyTetrahedron : public vcg::TetraATAVTQ<MyVertex,MyTetrahedron>{};

class MyTMesh: public vcg::tetra::Tetramesh< std::vector<MyVertex>, std::vector<MyTetrahedron > >{};


#include <vcg/complex/local_optimization/base.h>
#include <vcg/complex/local_optimization/tetra_edge_collapse.h>


vcg::LocalOptimization<MyTMesh> loc;
vcg::tetra::TetraEdgeCollapse<MyTMesh> c;
MyTMesh mesh;

int main(int,char**argv){

	loc.m = & mesh;

	vcg::tetra::io::ImporterTS<MyTMesh>::Open(mesh,argv[1]);
	printf("mesh loaded %d %d \n",mesh.vn,mesh.tn);
//	vcg::tetra::io::ExporterPLY<MyTMesh>::Save(mesh,(string(argv[1])+string(".ply")).c_str());
	


	vcg::tetra::TetraEdgeCollapse<MyTMesh> *_ ;
	bool res;
	do{
			vcg::tetra::UpdateTetraTopology<MyTMesh::VertexContainer,MyTMesh::TetraContainer>
		::VTTopology(mesh.vert,mesh.tetra);

	vcg::tetra::UpdateTetraTopology<MyTMesh::VertexContainer,MyTMesh::TetraContainer>
		::TTTopology(mesh.vert,mesh.tetra);

	vcg::tetra::UpdateTetraTopology<MyTMesh::VertexContainer,MyTMesh::TetraContainer>
		::setExternalVertices(mesh.vert,mesh.tetra);

	 _=new vcg::tetra::TetraEdgeCollapse<MyTMesh>();
		loc.h.push_back(vcg::LocalOptimization<MyTMesh>::HeapElem(_));

		loc.Init();
		loc.SetTargetSimplices(10);
		
		 res = loc.DoOptimization();

		printf("ood %d\n bor %d\n vol %d \n lkv %d \n lke %d \n lkf %d \n",
						FAIL::OFD(),
						FAIL::BOR(),
						FAIL::VOL(),
						FAIL::LKV(),
						FAIL::LKE(),
						FAIL::LKF()
			);
		printf("mesh  %d %d \n",mesh.vn,mesh.tn);
	}while(!res);


	vcg::tetra::io::ExporterPLY<MyTMesh>::Save(mesh,"out.ply");
	
	return 0;

}
