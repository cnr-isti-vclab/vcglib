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

//vcgLibForEmbree
#include<wrap/embree/EmbreeAdaptor.h>

using namespace vcg;
using namespace std;


int main( int argc, char **argv )
{
  MyMesh m;
  tri::io::ImporterOFF<MyMesh>::Open(m, "../ExampleMeshes/DragonHead.off");//argv[1]);//metti il rferimento ad una mesh
  int ret = tri::io::ImporterOFF<MyMesh>::Open(m,"../ExampleMeshes/DragonHead.off");
  if(ret!=tri::io::ImporterOFF<MyMesh>::NoError)
  {
    cout<<"Error reading file \n"<<endl;
    exit(0);
  }

  char *endptr;
  int nOfRays = 128; //strtof(argv[2], &endptr);
  
  MyMesh m2,m3,m4,m5;
  vcg::tri::Append<MyMesh,MyMesh>::MeshCopy(m2,m);
  vcg::tri::Append<MyMesh,MyMesh>::MeshCopy(m3,m);
  vcg::tri::Append<MyMesh,MyMesh>::MeshCopy(m4,m);

  EmbreeAdaptor<MyMesh> adaptor = EmbreeAdaptor<MyMesh>(m,1);
  adaptor.computeAmbientOcclusion(m,nOfRays);
  tri::UpdateQuality<MyMesh>::VertexFromFace(m);
  tri::UpdateColor<MyMesh>::PerVertexQualityGray(m);
  tri::io::ExporterOFF<MyMesh>::Save(m,"testAO.off",tri::io::Mask::IOM_VERTCOLOR);

  std::vector<Point3f> unifDirVec;
  std::vector<Point3f> ndir;
	GenNormal<float>::Fibonacci(nOfRays,unifDirVec);

  adaptor = EmbreeAdaptor<MyMesh>(m2,1);
  adaptor.computeAmbientOcclusion(m2,nOfRays,unifDirVec);
  tri::UpdateQuality<MyMesh>::VertexFromFace(m2);
  tri::UpdateColor<MyMesh>::PerVertexQualityGray(m2);
  tri::io::ExporterOFF<MyMesh>::Save(m2,"testAOM.off",tri::io::Mask::IOM_VERTCOLOR);

  adaptor = EmbreeAdaptor<MyMesh>(m3,4);
  adaptor.computeObscurance(m3,nOfRays,0.01);
  tri::UpdateQuality<MyMesh>::VertexFromFace(m3);
  tri::UpdateColor<MyMesh>::PerVertexQualityGray(m3);
  tri::io::ExporterOFF<MyMesh>::Save(m3,"testAODir.off",tri::io::Mask::IOM_VERTCOLOR);

  adaptor = EmbreeAdaptor<MyMesh>(m4,4);
  adaptor.computeSDF(m4,nOfRays);
  tri::UpdateQuality<MyMesh>::VertexFromFace(m4);
  tri::UpdateColor<MyMesh>::PerVertexQualityRamp(m4);
  tri::io::ExporterOFF<MyMesh>::Save(m4,"testSDF.off",tri::io::Mask::IOM_VERTCOLOR);

  adaptor = EmbreeAdaptor<MyMesh>(m5,4);
  vector<Point3f> BentNormal = adaptor.AOBentNormal(m5,nOfRays);

  
  cout << "done" << endl;
  return 0;
}