	
// mesh definition 
//#include <vcg/simplex/vertex/with/vn.h>
//#include <vcg/simplex/face/with/af.h>
//#include <vcg/complex/trimesh/base.h>

#include<vcg/simplex/vertex/base.h>
#include<vcg/simplex/face/base.h>
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
#include <iostream>
#include <vector>
#include <time.h>

using namespace vcg;
using namespace std;

class MyFace;
class MyVertex;

struct MyUsedTypes : public UsedTypes<	Use<MyVertex>		::AsVertexType,
																				Use<MyFace>			::AsFaceType>{};

class MyVertex  : public Vertex< MyUsedTypes, vertex::Coord3f, vertex::Normal3f, vertex::BitFlags, vertex::Mark>{};
class MyFace    : public Face  < MyUsedTypes, face::VertexRef, face::Normal3f, face::BitFlags > {};
class MyMesh    : public vcg::tri::TriMesh< vector<MyVertex>, vector<MyFace> > {};

bool callback(int percent, const char *str) {
  cout << "str: " << str << " " << percent << "%\n";
  return true;
}

int  main(int argc, char **argv)
{
 if(argc<3)
	{
		printf(
		"\n                  trimesh_ball_pivoting ("__DATE__")\n"
			"						Visual Computing Group I.S.T.I. C.N.R.\n"
      "Usage: trimesh_ball_pivoting filein.ply fileout.ply [opt]\n"
      "options: \n"
      "-r <val> radius of the rolling ball\n"
      "-c <val> clustering radius (as fraction of radius) default: 0.05\n"
			);
		exit(0);
	}

   float radius = 0.0f;
   float clustering = 0.05;
   int i = 3;
	while(i<argc)
		{
			if(argv[i][0]!='-')
				{printf("Error unable to parse option '%s'\n",argv[i]); exit(0);}
			switch(argv[i][1])
			{				
				case 'r' :	radius = atof(argv[++i]); printf("Using %f sphere radius\n",radius);  break;
				case 'c' :	clustering = atof(argv[++i]); printf("Using %f clustering radius\n",clustering); break;
      
				default : {printf("Error unable to parse option '%s'\n",argv[i]); exit(0);}
			}
			++i;
		}
    if(radius == 0) 
      printf("Autodetecting ball radius...\n");
                
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
  tri::BallPivoting<MyMesh> pivot(m, radius, clustering); 
  printf("Ball radius: %f\nClustering points withing %f radii\n", pivot.radius, clustering);

  int t1=clock();
  // the main processing
  pivot.BuildMesh(callback);

  int t2=clock();

  printf("Output mesh vn:%i fn:%i\n",m.vn,m.fn);
  printf("Created in :%i msec (%i+%i)\n",t2-t0,t1-t0,t2-t1);
	
  vcg::tri::io::PlyInfo pi;
  vcg::tri::io::ExporterPLY<MyMesh>::Save(m,argv[2],pi.mask);
  return 0;

}
