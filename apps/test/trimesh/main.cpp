#pragma warning(disable : 4267)
#include <stdio.h>
#include <vector>
#include<vcg/math/base.h>
#include<vcg/space/point3.h>
#include<vcg/space/point4.h>
#include<vcg/space/color4.h>
#include<vcg/simplex/vertex/with/n.h>
#include<vcg/simplex/face/with/fn.h>
#include<vcg/simplex/face/with/fa.h>
#include<vcg/complex/trimesh/base.h>
#include<vcg/complex/trimesh/platonic.h>
#include<vcg/complex/trimesh/update/flags.h>
#include<vcg/complex/trimesh/update/normal.h>
#include<vcg/complex/trimesh/update/topology.h>
#include<vcg/complex/trimesh/update/color.h>
#include<wrap/io_trimesh/export_stl.h>
#include<wrap/io_trimesh/import_ply.h>

using namespace vcg;
using namespace std;
//using namespace tri;

int main(int argc, char *argv[])
{
  printf("Hello Library!\n");
 
  VertexNf vnf;
  VertexNd vnd;

  FaceFN<VertexNf> fnf;
  FaceFN<VertexNd> fnd;

  typedef tri::TriMesh< vector<VertexNf>, vector< FaceFN<VertexNf> > > MyMesh; 
  MyMesh tm;

  typedef tri::TriMesh< vector<VertexNf>, vector< FaceFA<VertexNf> > > MyMeshAdj; 
  MyMeshAdj ta;

  printf("Sizeof(VertexNf) = %i (%i + %i + %i)\n",sizeof(VertexNf),sizeof(Point3f),sizeof(int),sizeof(Point3f));
  printf("Sizeof(VertexNd) = %i (%i + %i + %i) \n",sizeof(VertexNd),sizeof(Point3d),sizeof(int),sizeof(Point3d));
  printf("Sizeof(FaceFN<VertexNf>) = %i \n",sizeof(FaceFN<VertexNf>));
  printf("Sizeof(FaceFA<VertexNd>) = %i \n",sizeof(FaceFA<VertexNd>));
  
  tri::Tetrahedron(tm);

  tri::io::ExporterSTL<MyMesh>::Save(tm,"Tetra.stl",false);
  //tri::io::ExporterSTL<MyMesh>::Save(tm,"armawarp.stl",false);
  tri::UpdateNormals<MyMesh>::PerVertexNormalized(tm);
  tri::UpdateNormals<MyMesh>::PerFaceNormalized(tm);
  tri::UpdateNormals<MyMesh>::PerVertex(tm);
  tri::UpdateNormals<MyMesh>::PerFace(tm);

  tri::io::ImporterPLY<MyMeshAdj>::Open(ta,"armawarp.ply");
  printf("Loaded Mesh Has %i vn %i fn\n",ta.vn,ta.fn);
  tri::UpdateTopology<MyMeshAdj>::FaceFace(ta);
  tri::UpdateFlags<MyMeshAdj>::FaceBorderFromFF(ta);

  return 0;
}