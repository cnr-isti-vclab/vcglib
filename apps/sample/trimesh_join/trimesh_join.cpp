#include <vector>

#include<vcg/simplex/vertexplus/base.h>
#include<vcg/simplex/faceplus/base.h>
#include<vcg/simplex/face/topology.h>

#include<vcg/complex/trimesh/base.h>
#include<vcg/complex/trimesh/append.h>
#include<vcg/complex/trimesh/clean.h>

// input output
#include <wrap/io_trimesh/import_ply.h>
#include <wrap/io_trimesh/export_ply.h>

// std
#include <vector>

using namespace vcg;
using namespace std;

class MyEdge;    // dummy prototype never used
class MyFace;
class MyVertex;

class MyVertex  : public VertexSimp2< MyVertex, MyEdge, MyFace, vert::Coord3f, vert::BitFlags  >{};
class MyFace    : public FaceSimp2  < MyVertex, MyEdge, MyFace, face::VertexRef, face::BitFlags > {};
class MyMesh    : public vcg::tri::TriMesh< vector<MyVertex>, vector<MyFace> > {};



int main(int argc,char **argv )
{
 if(argc<2)
	{
		printf(		"\n trimesh_join ("__DATE__")\n"
		        	  "Visual Computing Group I.S.T.I. C.N.R.\n"
                "Usage: trimesh_join filename.ply [filename.ply | *] \n"
			);
		exit(0);
	}

  MyMesh ml,mr;

	int i=1; 
	while(i<argc)
		{
		  if(vcg::tri::io::ImporterPLY<MyMesh>::Open(mr,argv[i])!=0)
		  {
        printf("Error reading file  %s\n",argv[1]);
			  exit(0);
		  }
      printf("Input mesh %3i vn:%9i fn:%9i\n",i, mr.vn, mr.fn);
      tri::Append<MyMesh,MyMesh>::Mesh(ml,mr); // append mesh mr to ml
  		++i;
		}
  
  printf("Output mesh vn:%i fn:%i\n",ml.vn,ml.fn);
	
  tri::io::ExporterPLY<MyMesh>::Save(ml,"joined.ply");
  int dv=tri::Clean<MyMesh>::RemoveDuplicateVertex(ml); 
  printf("Removed %i duplicated vertices\n",dv);
  tri::io::ExporterPLY<MyMesh>::Save(ml,"joined_unif.ply");
}

