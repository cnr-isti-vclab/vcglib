#pragma warning(disable : 4267)
#include <stdio.h>
#include <vector>
#include<vcg/math/base.h>
#include<vcg/space/point3.h>
#include<vcg/space/point4.h>
#include<vcg/space/color4.h>
#include<vcg/simplex/vertex/with/n.h>
#include<vcg/simplex/vertex/with/cn.h>
#include<vcg/simplex/face/with/fn.h>
#include<vcg/simplex/face/with/fafc.h>
#include<vcg/simplex/face/pos.h>
#include<vcg/complex/trimesh/base.h>
#include<vcg/complex/trimesh/platonic.h>
#include<vcg/complex/trimesh/update/flag.h>
#include<vcg/complex/trimesh/update/normal.h>
#include<vcg/complex/trimesh/update/topology.h>
#include<vcg/complex/trimesh/update/color.h>
#include<wrap/io_trimesh/export_stl.h>
#include<wrap/io_trimesh/import_ply.h>
#include<wrap/io_trimesh/import_stl.h>
#include<wrap/io_trimesh/export_ply.h>

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


  typedef tri::TriMesh< vector<VertexCNf>, vector< FaceFAFC<VertexCNf> > > MyMeshAdj; 
  MyMeshAdj ta;

  printf("Sizeof(VertexNf) = %i (%i + %i + %i)\n",sizeof(VertexNf),sizeof(Point3f),sizeof(int),sizeof(Point3f));
  printf("Sizeof(VertexNd) = %i (%i + %i + %i) \n",sizeof(VertexNd),sizeof(Point3d),sizeof(int),sizeof(Point3d));
  printf("Sizeof(FaceFN<VertexNf>) = %i \n",sizeof(FaceFN<VertexNf>));
  //printf("Sizeof(FaceFA<VertexCNd>) = %i \n",sizeof(FaceFA<VertexCNd>));
  
  typedef tri::TriMesh< vector<VertexNf>, vector< FaceFN<VertexNf> > > MyMesh; 
  MyMesh tm;
  tri::Tetrahedron(tm);
  tri::io::ExporterPLY<MyMesh>::Save(tm,"Tetrabin.ply");
  tri::Octahedron(tm);
  tri::io::ExporterPLY<MyMesh>::Save(tm,"Octbin.ply");
  tri::Icosahedron(tm);
  tri::io::ExporterPLY<MyMesh>::Save(tm,"Icobin.ply");
  
  //tri::io::ExporterSTL<MyMesh>::Save(tm,"Tetra.stl",false);
  //tri::io::ExporterPLY<MyMesh>::Save(tm,"Tetraascii.ply",false);
  //tri::io::ExporterSTL<MyMesh>::Save(tm,"armawarp.stl",false);
  //tri::UpdateNormals<MyMesh>::PerVertexNormalized(tm);
  //tri::UpdateNormals<MyMesh>::PerFaceNormalized(tm);
  //tri::UpdateNormals<MyMesh>::PerVertex(tm);
  //tri::UpdateNormals<MyMesh>::PerFace(tm);

  /*tri::io::ImporterPLY<MyMeshAdj>::Open(ta,"bigtest.ply");
  printf("Loaded Mesh Has %i vn %i fn\n",ta.vn,ta.fn);
  tri::UpdateTopology<MyMeshAdj>::FaceFace(ta);
  tri::UpdateFlags<MyMeshAdj>::FaceBorderFromFF(ta);
  tri::UpdateColor<MyMeshAdj>::FaceBF(ta);
  tri::io::PlyInfo pi;
  pi.mask=tri::io::PLYMask::PM_FACECOLOR;
  tri::io::ExporterPLY<MyMeshAdj>::Save(ta,"color.ply",true, pi);*/

  face::Pos<MyMeshAdj::FaceType> fp;

  tri::io::ImporterSTL<MyMeshAdj>::Open(ta,"knotbin.stl");
  tri::io::ExporterPLY<MyMeshAdj>::Save(ta,"knotbin.ply");

  tri::io::ImporterSTL<MyMeshAdj>::Open(ta,"knotasc.stl");
  tri::io::ExporterPLY<MyMeshAdj>::Save(ta,"knotasc.ply");

  return 0;
}