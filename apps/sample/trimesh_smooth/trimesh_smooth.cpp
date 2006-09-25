#include <vector>

#include<vcg/simplex/vertexplus/base.h>
#include<vcg/simplex/faceplus/base.h>
#include<vcg/simplex/face/topology.h>

#include<vcg/complex/trimesh/base.h>

#include <vcg/complex/trimesh/update/topology.h>
#include <vcg/complex/trimesh/update/normal.h>
// to clean up a mesh
#include<vcg/complex/trimesh/clean.h>
#include<vcg/complex/trimesh/smooth.h>

// input output
#include <wrap/io_trimesh/import.h>
#include <wrap/io_trimesh/export_ply.h>


using namespace vcg;
using namespace std;


class MyEdge;    // dummy prototype never used
class MyFace;
class MyVertex;

class MyVertex  : public VertexSimp2< MyVertex, MyEdge, MyFace, vert::VFAdj, vert::Coord3f, vert::Normal3f, vert::BitFlags  >{};
class MyFace    : public FaceSimp2  < MyVertex, MyEdge, MyFace, face::VFAdj, face::Normal3f, face::VertexRef, face::BitFlags > {};
class MyMesh    : public vcg::tri::TriMesh<vector<MyVertex>, vector<MyFace> > {};

//class MyFace;
//class MyEdge;
//class MyVertex:public Vertex<float,MyEdge,MyFace>{};
//class MyFace :public FaceFN<MyVertex,MyEdge,MyFace>{};
//class MyMesh: public tri::TriMesh< std::vector<MyVertex>, std::vector<MyFace > >{};
//

int main(int argc,char ** argv)
{
if(argc<3) 
{
  printf("Usage: trimesh_smooth <filename> <steps> <sigma> <fitstep>\n");
  return 0;
}

	MyMesh m;

	//open a mesh
	int err = tri::io::Importer<MyMesh>::Open(m,argv[1]);
  if(err) {
      printf("Error in reading %s: '%s'\n",argv[1],tri::io::Importer<MyMesh>::ErrorMsg(err));
      exit(-1);
    }
  // some cleaning to get rid of bad file formats like stl that duplicate vertexes..
  int dup = tri::Clean<MyMesh>::RemoveDuplicateVertex(m);
  int unref =  tri::Clean<MyMesh>::RemoveUnreferencedVertex(m);
  printf("Removed %i duplicate and %i unreferenced vertices from mesh %s\n",dup,unref,argv[1]);
  int Step= atoi(argv[2]);
  tri::UpdateTopology<MyMesh>::VertexFace(m);
  tri::UpdateNormals<MyMesh>::PerFaceNormalized(m);
 
  for(int i=0;i<Step;++i)
  {
    tri::UpdateNormals<MyMesh>::PerFaceNormalized(m);
    PasoDobleSmoothFast(m,atoi(argv[3]),atof(argv[4]),atoi(argv[5]));
  }

  //LaplacianSmooth(m,atoi(argv[2]));
  tri::io::ExporterPLY<MyMesh>::Save(m,"out.ply");

  return 0;
}

