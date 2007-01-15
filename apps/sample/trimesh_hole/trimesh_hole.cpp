#include <vector>
#include <iostream>

#include<vcg/simplex/vertexplus/base.h>
#include<vcg/simplex/faceplus/base.h>
#include<vcg/simplex/face/topology.h>
#include<vcg/complex/trimesh/base.h>
#include <vcg/complex/local_optimization.h>
//#include <relaxation.h>
#include<vcg/complex/local_optimization/tri_edge_flip.h>


// topology computation
#include<vcg/complex/trimesh/update/topology.h>
#include <vcg/complex/trimesh/update/flag.h>
#include <vcg/complex/trimesh/update/normal.h>

// half edge iterators
#include<vcg/simplex/face/pos.h>

#include<vcg/complex/trimesh/hole.h>

// input output
#include <wrap/io_trimesh/import_ply.h>
#include <wrap/io_trimesh/export_ply.h>



using namespace vcg;
using namespace std;

class MyEdge;    // dummy prototype never used
class MyFace;
class MyVertex;

class MyVertex  : public VertexSimp2< MyVertex, MyEdge, MyFace, vert::Coord3f, vert::BitFlags, vert::Normal3f, vert::Mark, vert::Color4b >{};
class MyFace    : public FaceSimp2  < MyVertex, MyEdge, MyFace, face::VertexRef,face::FFAdj, face::Mark, face::BitFlags, face::Normal3f> {};

class MyMesh : public tri::TriMesh< vector<MyVertex>, vector<MyFace > >{};



/*
//for coplanar mesh
class MyTriEdgeFlip: public vcg::tri::PlanarEdgeFlip< MyMesh, MyTriEdgeFlip > {
						public:
						typedef  vcg::tri::PlanarEdgeFlip< MyMesh,  MyTriEdgeFlip > TEF;
            inline MyTriEdgeFlip(  const TEF::PosType &p, int i) :TEF(p,i){}
};
/*/
//Delaunay
class MyTriEdgeFlip: public vcg::tri::TriEdgeFlip< MyMesh, MyTriEdgeFlip > {
						public:
						typedef  vcg::tri::TriEdgeFlip< MyMesh,  MyTriEdgeFlip > TEF;
            inline MyTriEdgeFlip(  const TEF::PosType &p, int i) :TEF(p,i){}
};
//*/
bool callback(int percent, const char *str) {
  cout << "str: " << str << " " << percent << "%\r";
  return true;
}

int main(int argc,char ** argv){

	if(argc<5)
	{
		printf(
			"\n     HoleFilling ("__DATE__")\n"
			"Visual Computing Group I.S.T.I. C.N.R.\n"
      "Usage: trimesh_hole #algorithm #size filein.ply fileout.ply \n"
			"#algorithm: \n"
			" 1) Trivial Ear \n"
			" 2) Leipa Ear \n"
			" 3) Selfintersection Ear \n"
			" 4) Minimum weight \n"
			);
		exit(0);
	}

	int algorithm = atoi(argv[1]);
	int holeSize  = atoi(argv[2]);
	if(algorithm < 0 && algorithm > 4)
	{
		printf("Error in algorithm's selection\n",algorithm);
		exit(0);
	}

	MyMesh m;

	if(tri::io::ImporterPLY<MyMesh>::Open(m,argv[3])!=0)
	{
		printf("Error reading file  %s\n",argv[2]);
		exit(0);
	}


	//update the face-face topology 
	tri::UpdateTopology<MyMesh>::FaceFace(m);
	tri::UpdateNormals<MyMesh>::PerVertexPerFace(m);
	tri::UpdateFlags<MyMesh>::FaceBorderFromFF(m);
  assert(tri::Clean<MyMesh>::IsFFAdjacencyConsistent(m));
	

  tri::Hole<MyMesh> holeFiller;
	switch(algorithm)
	{
  case 1:			tri::Hole<MyMesh>::EarCuttingFill<tri::TrivialEar<MyMesh> >(m,holeSize,false);                	        break;
  case 2:   	tri::Hole<MyMesh>::EarCuttingFill<tri::MinimumWeightEar< MyMesh> >(m,holeSize,false,callback);          break;
  case 3: 		tri::Hole<MyMesh>::EarCuttingIntersectionFill<tri::SelfIntersectionEar< MyMesh> >(m,holeSize,false);		break;
  case 4: 		tri::Hole<MyMesh>::MinimumWeightFill(m, false);		                                                      break;
	}
  printf("\nCompleted. Saving...\n");
  assert(tri::Clean<MyMesh>::IsFFAdjacencyConsistent(m));
	tri::io::ExporterPLY<MyMesh>::Save(m,argv[4],false);
  printf("\nStart flipping...\n");

  vcg::LocalOptimization<MyMesh> FlippingSession(m);
	FlippingSession.SetTargetMetric(-0.000000000001f);//
  FlippingSession.Init<MyTriEdgeFlip >();
  FlippingSession.DoOptimization();
  tri::io::ExporterPLY<MyMesh>::Save(m,"out2.ply",vcg::tri::io::Mask::IOM_VERTCOLOR);
  
	return 0;
}

