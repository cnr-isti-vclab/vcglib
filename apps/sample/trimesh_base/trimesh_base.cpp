#include <vector>

#include<vcg/simplex/vertex/base.h>//
#include<vcg/simplex/vertex/component.h>
#include<vcg/simplex/face/base.h>//
#include<vcg/simplex/face/component.h>
#include<vcg/simplex/face/topology.h>//
#include<vcg/complex/trimesh/base.h>//

// input output
#include<wrap/io_trimesh/import.h>
#include<wrap/io_trimesh/export.h>//just in case

// topology computation
#include<vcg/complex/trimesh/update/topology.h>//
#include<vcg/complex/trimesh/update/flag.h>//

// half edge iterators
//#include<vcg/simplex/face/pos.h>

// normals and curvature
#include<vcg/complex/trimesh/update/normal.h> //class UpdateNormals 
#include<vcg/complex/trimesh/update/curvature.h> //class curvature

using namespace vcg;
using namespace std;

class MyEdge; // dummy prototype
class MyFace;
class MyVertex;

class MyVertex  : public VertexSimp2< MyVertex, MyEdge, MyFace, vertex::Coord3f, vertex::Normal3f, vertex::BitFlags  >{};
class MyFace    : public FaceSimp2  < MyVertex, MyEdge, MyFace, face::FFAdj,  face::VertexRef, face::BitFlags > {};
class MyMesh    : public vcg::tri::TriMesh< vector<MyVertex>, vector<MyFace> > {};

int main( int argc, char **argv ) {
MyMesh m;
// this is the section with problems
if(vcg::tri::io::ImporterPLY<MyMesh>::Open(m,argv[1])!=0)
{
printf("Error reading file  %s\n",argv[1]);
exit(0);
} // from here no problems

vcg::tri::UpdateTopology<MyMesh>::FaceFace(m);
vcg::tri::UpdateFlags<MyMesh>::FaceBorderFromFF(m);
vcg::tri::UpdateNormals<MyMesh>::PerVertexNormalized(m);
printf("Input mesh  vn:%i fn:%i\n",m.vn,m.fn);
printf( "Mesh has %i vert and %i faces\n", m.vn, m.fn );

return 0;
}