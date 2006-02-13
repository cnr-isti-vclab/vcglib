#include <vector>

using namespace std;

// VCG headers for triangular mesh processing
#include <vcg/complex/trimesh/base.h>
#include <vcg/complex/trimesh/update/topology.h>
#include <vcg/complex/trimesh/update/edges.h>
#include <vcg/complex/trimesh/update/bounding.h>
#include <vcg/complex/trimesh/update/flag.h>
#include <vcg/complex/trimesh/clean.h>

// VCG File Format Importer
#include <wrap/io_trimesh/import.h>
#include <wrap/io_trimesh/export.h>

// VCG Vertex
#include <vcg/simplex/vertexplus/base.h>
#include <vcg/simplex/vertexplus/component.h>

// VCG Faces
#include <vcg/simplex/faceplus/base.h>
#include <vcg/simplex/faceplus/component.h>

#include <vcg/space/index/grid_static_ptr.h>

using namespace vcg;

class MyFace;
class MyEdge;
class MyVertex  : public VertexSimp2< MyVertex, MyEdge, MyFace, vert::VFAdj, vert::Coord3f, 
	vert::BitFlags, vert::Normal3f > {};
class MyFace    : public FaceSimp2< MyVertex, MyEdge, MyFace, face::FFAdj, face::VFAdj, 
	face::VertexRef, face::Normal3f, face::BitFlags, face::Mark > {};
class MyMesh    : public vcg::tri::TriMesh< vector<CVertex>, vector<CFace> > {};

typedef MyMesh::VertexPointer VertexPointer;
typedef MyMesh::VertexIterator VertexIterator;

typedef MyMesh::VertexPointer VertexPointer;
typedef MyMesh::VertexIterator VertexIterator;
typedef MyMesh::FaceContainer FaceContainer;


int main(int argc,char ** argv)
{
	if (argc<6) 
	{
		printf("Usage: trimesh_intersection <filename> <a> <b> <c> <d>\n\n");
		printf("       <filename>        Mesh model to intersect (PLY format).");
		printf("       <a> <b> <c> <d>   The coefficients that specifying a plane in the form:\n");
		printf("                             a*x + b*y + c*z + d = 0\n");

		printf("       Example: trimesh_intersection bunny.ply 1.0 0.0 0.0 0.0");

		return 0;
	}

	MyMesh m;

	// open a mesh
	int err = tri::io::Importer<MyMesh>::Open(m,argv[1]);
  if(err) 
	{
      printf("Error in reading %s: '%s'\n",argv[1],tri::io::Importer<MyMesh>::ErrorMsg(err));
      exit(-1);  
	}

  // some cleaning to get rid of bad file formats like stl that duplicate vertexes..
  int dup = tri::Clean<MyMesh>::RemoveDuplicateVertex(m);
  int unref =  tri::Clean<MyMesh>::RemoveUnreferencedVertex(m);

	if (dup > 0 || unref > 0)
		printf("Removed %i duplicate and %i unreferenced vertices from mesh %s\n",dup,unref,argv[1]);

  printf("");

  // Compute cross-intersection with the given plane
	/////////////////////////////////////////////////////////

	double a = atof(arg[2]);
	double b = atof(arg[3]);
	double c = atof(arg[4]);
	double d = atof(arg[5]);

	// export cross-section
  tri::io::ExporterPLY<MyMesh>::Save(m, "out.ply");

  return 0;
}

