	
// mesh definition 
//#include <vcg/simplex/vertex/with/vn.h>
//#include <vcg/simplex/face/with/af.h>
//#include <vcg/complex/trimesh/base.h>

#include<vcg/simplex/vertexplus/base.h>
#include<vcg/simplex/faceplus/base.h>
#include<vcg/simplex/face/topology.h>

#include<vcg/complex/trimesh/base.h>

#include <vcg/complex/trimesh/update/bounding.h>
#include <vcg/complex/trimesh/update/topology.h>
#include <vcg/complex/trimesh/update/normal.h>
#include <vcg/complex/trimesh/update/flag.h>
#include <vcg/complex/trimesh/create/ball_pivoting.h>

// input output
#include <wrap/io_trimesh/import_ply.h>
#include <wrap/io_trimesh/export_ply.h>

// std
#include <vector>
#include <time.h>

using namespace vcg;
using namespace std;

class MyEdge;    // dummy prototype never used
class MyFace;
class MyVertex;

class MyVertex  : public VertexSimp2< MyVertex, MyEdge, MyFace, vert::Coord3f, vert::Normal3f, vert::BitFlags  >{};
class MyFace    : public FaceSimp2  < MyVertex, MyEdge, MyFace, face::VertexRef, face::Normal3f, face::BitFlags > {};
class MyMesh    : public vcg::tri::TriMesh< vector<MyVertex>, vector<MyFace> > {};

int  main(int argc, char **argv)
{
 if(argc<3)
	{
		printf(
		"\n                  trimesh_ball_pivoting ("__DATE__")\n"
			"						Visual Computing Group I.S.T.I. C.N.R.\n"
      "Usage: PlyRefine filein.ply fileout.ply [opt] \n"
      "options: \n"
      "-r <val> radius of the rolling ball\n"
			);
		exit(0);
	}

  int i=3;  
  int CellNum=100000; 
  float CellSize=0; 
  bool DupFace=false;

	while(i<argc)
		{
			if(argv[i][0]!='-')
				{printf("Error unable to parse option '%s'\n",argv[i]); exit(0);}
			switch(argv[i][1])
			{				
				case 'k' :	CellNum=atoi(argv[i+1]); ++i; printf("Using %i clustering cells\n",CellNum); break;
				case 's' :	CellSize=atof(argv[i+1]); ++i; printf("Using %5f as clustering cell size\n",CellSize); break;
				case 'd' :	DupFace=true; printf("Enabling the duplication of faces for double surfaces\n"); break;

				default : {printf("Error unable to parse option '%s'\n",argv[i]); exit(0);}
			}
			++i;
		}

	MyMesh m;

	if(vcg::tri::io::ImporterPLY<MyMesh>::Open(m,argv[1])!=0)
		{
      printf("Error reading file  %s\n",argv[1]);
			exit(0);
		}
  vcg::tri::UpdateBounding<MyMesh>::Box(m);
  vcg::tri::UpdateNormals<MyMesh>::PerFace(m);
  printf("Input mesh  vn:%i fn:%i\n",m.vn,m.fn);

  int t0=clock();
  // Initialization

  int t1=clock();
  // the main processing

  int t2=clock();

  printf("Output mesh vn:%i fn:%i\n",m.vn,m.fn);
  printf("Simplified in :%i msec (%i+%i)\n",t2-t0,t1-t0,t2-t1);
	
  vcg::tri::io::PlyInfo pi;
	vcg::tri::io::ExporterPLY<MyMesh>::Save(m,argv[2],pi.mask);
	return 0;

}
