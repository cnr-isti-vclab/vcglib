/*
 echo $LD_LIBRARY_PATH
  LD_LIBRARY_PATH=/mnt/c/Users/super/Dropbox/3DProcessing/project3D/embree-3.13.3.x86_64.linux/lib:$LD_LIBRARY_PATH

  LD_LIBRARY_PATH=/mnt/e/UniversitàMagistrale/secondoSemestre/3DgeometricModelingProcessing/vcglib/wrap/embree/embree-3.13.3.x86_64.linux/lib:$LD_LIBRARY_PATH

  export LD_LIBRARY_PATH
  g++ ./wrap/embree/vcgForEmbree.cpp -o prova.o -lembree3 -I ./vcg -I ./ -I ./eigenlib -I ./wrap/embree/embree-3.13.3.x86_64.linux/include -L ./wrap/embree/embree-3.13.3.x86_64.linux/lib -std=c++11
  
  g++ ./wrap/embree/sample/embree_sample.cpp -o prova.o  -lembree3 -I ./vcg -I ./ -I ./eigenlib -I ./wrap/embree/embree-3.13.3.x86_64.linux/include -L ./wrap/embree/embree-3.13.3.x86_64.linux/lib -std=c++11 -fopenmp -O3
  ./prova.o ./wrap/embree/sample/ExampleMeshes/bunny10k.off 32 false
*/
#include <iostream>

#include <vcg/complex/complex.h>

//import export
#include<wrap/io_trimesh/import_ply.h>
#include<wrap/io_trimesh/export_ply.h>
#include <wrap/io_trimesh/export_off.h>
#include <wrap/io_trimesh/import_off.h>

#include <time.h>
#include <vcg/math/gen_normal.h>
#include <vcg/complex/allocate.h>

//vcgLibForEmbree
#include<wrap/embree/EmbreeAdaptor.h>

class MyVertex; class MyEdge; class MyFace;
struct MyUsedTypes : public vcg::UsedTypes<vcg::Use<MyVertex>   ::AsVertexType,
    vcg::Use<MyEdge>     ::AsEdgeType,
    vcg::Use<MyFace>     ::AsFaceType> {};

class MyVertex : public vcg::Vertex< MyUsedTypes, vcg::vertex::Coord3f, vcg::vertex::Normal3f, vcg::vertex::BitFlags, vcg::vertex::VFAdj, vcg::vertex::Qualityf, vcg::vertex::Color4b> {};
class MyFace : public vcg::Face<   MyUsedTypes, vcg::face::FFAdj, vcg::face::VFAdj, vcg::face::Normal3f, vcg::face::VertexRef, vcg::face::BitFlags, vcg::face::Color4b, vcg::face::Qualityf> {};
class MyEdge : public vcg::Edge<   MyUsedTypes> {};

class MyMesh : public vcg::tri::TriMesh< std::vector<MyVertex>, std::vector<MyFace>, std::vector<MyEdge>  > {};

using namespace vcg;
using namespace std;


int main( int argc, char **argv )
{
    cout << "start" << endl;
  MyMesh m;
  tri::io::ImporterOFF<MyMesh>::Open(m, "../ExampleMeshes/abominevole.off");//argv[1]);//metti il rferimento ad una mesh
  int ret = tri::io::ImporterOFF<MyMesh>::Open(m,"../ExampleMeshes/abominevole.off");
  if(ret!=tri::io::ImporterOFF<MyMesh>::NoError)
  {
    cout<<"Error reading file \n"<<endl;
    exit(0);
  }

  char *endptr;
  int nOfRays = 150; //strtof(argv[2], &endptr);
  
  MyMesh m2,m3,m4,m5,m6;
  vcg::tri::Append<MyMesh,MyMesh>::MeshCopy(m2,m);
  vcg::tri::Append<MyMesh,MyMesh>::MeshCopy(m3,m);
  vcg::tri::Append<MyMesh,MyMesh>::MeshCopy(m4, m);
  vcg::tri::Append<MyMesh, MyMesh>::MeshCopy(m5, m);
  vcg::tri::Append<MyMesh,MyMesh>::MeshCopy(m6,m);

  EmbreeAdaptor<MyMesh> adaptor = EmbreeAdaptor<MyMesh>(m,8);
  adaptor.computeAmbientOcclusion(m,nOfRays);
  tri::UpdateQuality<MyMesh>::VertexFromFace(m);
  tri::UpdateColor<MyMesh>::PerVertexQualityGray(m);
  tri::UpdateNormal<MyMesh>::NormalizePerVertex(m);
  tri::io::ExporterOFF<MyMesh>::Save(m,"testAO.off",tri::io::Mask::IOM_VERTCOLOR);
 
  cout << "Done AO" << endl;

  std::vector<Point3f> unifDirVec;
  std::vector<Point3f> ndir;
	GenNormal<float>::Fibonacci(nOfRays,unifDirVec);
    Point3f dir(0, 1, 0);

    for (int g = 0; g < nOfRays; g++) {
        if (unifDirVec.at(g) >= dir) {
            ndir.push_back(unifDirVec.at(g));
        }
    }
  adaptor = EmbreeAdaptor<MyMesh>(m2,8);
  adaptor.computeAmbientOcclusion(m2,ndir);
  tri::UpdateQuality<MyMesh>::VertexFromFace(m2);
  tri::UpdateColor<MyMesh>::PerVertexQualityGray(m2);
  tri::io::ExporterOFF<MyMesh>::Save(m2,"testAODir.off",tri::io::Mask::IOM_VERTCOLOR);

  cout << "Done AO Directioned" << endl;
  
  EmbreeAdaptor<MyMesh> adaptor2 = EmbreeAdaptor<MyMesh>(m4,8);
  adaptor2.computeSDF(m4,nOfRays,90);
  tri::UpdateQuality<MyMesh>::VertexFromFace(m4);
  tri::UpdateColor<MyMesh>::PerVertexQualityRamp(m4);
  tri::io::ExporterOFF<MyMesh>::Save(m4,"testSDF.off",tri::io::Mask::IOM_VERTCOLOR);
  
  cout << "Done SDF" << endl;

  adaptor = EmbreeAdaptor<MyMesh>(m5,8);
  adaptor.computeNormalAnalysis(m5, nOfRays);
  tri::io::ExporterOFF<MyMesh>::Save(m5, "testNormal.off", tri::io::Mask::IOM_FACENORMAL);
  //vector<Point3f> BentNormal = adaptor.AOBentNormal(m5,nOfRays);

  cout << "Done NormalAnlysis" << endl;

  adaptor = EmbreeAdaptor<MyMesh>(m6, 4);
  Point3f p(1, 0, 0);
  adaptor.selectVisibleFaces(m6, p);
  tri::io::ExporterOFF<MyMesh>::Save(m6, "testSelectS.off", tri::io::Mask::IOM_FACECOLOR);
  
  cout << "done face selection" << endl;
  cout << "Done All" << endl;
  return 0;
}