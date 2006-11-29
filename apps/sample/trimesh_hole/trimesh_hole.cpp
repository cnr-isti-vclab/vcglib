#include <vector>

#include<vcg/simplex/vertexplus/base.h>
#include<vcg/simplex/faceplus/base.h>
#include<vcg/simplex/face/topology.h>
#include<vcg/complex/trimesh/base.h>

// topology computation
#include<vcg/complex/trimesh/update/topology.h>
#include <vcg/complex/trimesh/update/flag.h>
#include <vcg/complex/trimesh/update/normal.h>

// half edge iterators
#include<vcg/simplex/face/pos.h>

// input output
#include <wrap/io_trimesh/import_ply.h>
#include <wrap/io_trimesh/export_ply.h>

#include<vcg/complex/trimesh/hole.h>

using namespace vcg;

class MyEdge;    // dummy prototype never used
class MyFace;
class MyVertex;

class MyVertex  : public VertexSimp2< MyVertex, MyEdge, MyFace, vert::Coord3f, vert::BitFlags, vert::Normal3f  >{};
class MyFace    : public FaceSimp2  < MyVertex, MyEdge, MyFace, face::VertexRef,face::FFAdj, face::Mark, face::BitFlags, face::Normal3f > {};

class MyMesh : public tri::TriMesh< std::vector<MyVertex>, std::vector<MyFace > >{};

int main(int argc,char ** argv){

	if(argc<4)
	{
		printf(
			"\n     HoleFilling ("__DATE__")\n"
			"Visual Computing Group I.S.T.I. C.N.R.\n"
			"Usage: trimesh_hole #algorithm filein.ply fileout.ply \n"
			"#algorithm: \n"
			" 1) Trivial Ear \n"
			" 2) Leipa Ear \n"
			" 3) Selfintersection Ear \n"
			" 4) Minimum weight \n"
			);
		exit(0);
	}

	int algorithm = atoi(argv[1]);
	if(algorithm < 0 && algorithm > 4)
	{
		printf("Error in algorithm's selection\n",algorithm);
		exit(0);
	}

	MyMesh m;

	if(vcg::tri::io::ImporterPLY<MyMesh>::Open(m,argv[2])!=0)
	{
		printf("Error reading file  %s\n",argv[2]);
		exit(0);
	}


	//update the face-face topology 
	vcg::tri::UpdateTopology<MyMesh>::FaceFace(m);
	vcg::tri::UpdateNormals<MyMesh>::PerVertex(m);
	tri::UpdateFlags<MyMesh>::FaceBorderFromFF(m);

vcg::tri::Hole<MyMesh> holeFiller;

	switch(algorithm)
	{
	case 1:
		
		holeFiller.EarCuttingFill<vcg::tri::TrivialEar<MyMesh> >(m,50,false);
		break;
	case 2: 
		holeFiller.EarCuttingFill<vcg::tri::MinimumWeightEar< MyMesh> >(m,500,false);
		break;
	case 3:
		holeFiller.EarCuttingIntersectionFill<vcg::tri::SelfIntersectionEar< MyMesh> >(m,500,false);
		break;
	case 4:
		holeFiller.MinimumWeightFill(m, false);
		break;
	}

	vcg::tri::io::ExporterPLY<MyMesh>::Save(m,argv[3],false);
	return 0;
}

