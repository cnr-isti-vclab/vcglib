/*
 echo $LD_LIBRARY_PATH
  LD_LIBRARY_PATH=/mnt/c/Users/super/Dropbox/3DProcessing/project3D/embree-3.13.3.x86_64.linux/lib:$LD_LIBRARY_PATH

  LD_LIBRARY_PATH=/mnt/e/UniversitàMagistrale/secondoSemestre/3DgeometricModelingProcessing/vcglib/wrap/embree/embree-3.13.3.x86_64.linux/lib:$LD_LIBRARY_PATH

  export LD_LIBRARY_PATH
  g++ ./wrap/embree/vcgForEmbree.cpp -o prova.o -lembree3 -I ./vcg -I ./ -I ./eigenlib -I ./wrap/embree/embree-3.13.3.x86_64.linux/include -L ./wrap/embree/embree-3.13.3.x86_64.linux/lib -std=c++11
  
  g++ ./wrap/embree/sample/embree_sample.cpp -o prova.o  -lembree3 -I ./vcg -I ./ -I ./eigenlib -I ./wrap/embree/embree-3.13.3.x86_64.linux/include -L ./wrap/embree/embree-3.13.3.x86_64.linux/lib -std=c++11
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
#include <math/gen_normal.h>

//vcgLibForEmbree
#include<wrap/embree/vcgForEmbree.h>

using namespace vcg;
using namespace std;
using namespace vcgEmb;

int main( int argc, char **argv )
{
  MyMesh m;
  tri::io::ImporterOFF<MyMesh>::Open(m,argv[1]);//metti il rferimento ad una mesh
  int ret = tri::io::ImporterOFF<MyMesh>::Open(m,argv[1]);
  if(ret!=tri::io::ImporterOFF<MyMesh>::NoError)
  {
    cout<<"Error reading file \n"<<endl;
    exit(0);
  }

  char *endptr;
  int nOfRays = strtof(argv[2], &endptr);
  
  MyMesh m2,m3,m4;
  vcg::tri::Append<MyMesh,MyMesh>::MeshCopy(m2,m);
  vcg::tri::Append<MyMesh,MyMesh>::MeshCopy(m3,m);
  vcg::tri::Append<MyMesh,MyMesh>::MeshCopy(m4,m);

  clock_t t;
  t = clock();
  vcgForEmbree::embreeAmbientOcclusion(m,nOfRays,1);
  cout<<"AO 01t "<< (float) (clock()-t)/CLOCKS_PER_SEC<<endl;
  tri::UpdateQuality<MyMesh>::VertexFromFace(m);
  tri::UpdateColor<MyMesh>::PerVertexQualityGray(m);
  tri::io::ExporterOFF<MyMesh>::Save(m,"testAO.off",tri::io::Mask::IOM_VERTCOLOR);

  vector<Point3f> BentNormal = vcgEmb::vcgForEmbree::AOBentNormal(m2,nOfRays,1);
  
  t = clock();
  vcgEmb::vcgForEmbree::embreeAmbientOcclusion(m2,nOfRays,8);
  cout<<"AO 08t "<< (float) (clock()-t)/CLOCKS_PER_SEC<<endl;
  tri::UpdateQuality<MyMesh>::VertexFromFace(m2);
  tri::UpdateColor<MyMesh>::PerVertexQualityGray(m2);
  tri::io::ExporterOFF<MyMesh>::Save(m2,"testAOM.off",tri::io::Mask::IOM_VERTCOLOR);

  t = clock();
  vcgEmb::vcgForEmbree::embreeAmbientOcclusion(m2,nOfRays,16);
  cout<<"AO 16t "<< (float) (clock()-t)/CLOCKS_PER_SEC<<endl;

  std::vector<Point3f> unifDirVec;
  std::vector<Point3f> ndir;
	GenNormal<float>::Fibonacci(nOfRays,unifDirVec);
  Point3f dir(0,1,0);
  
  for(int g = 0; g<nOfRays; g++){   
     if(unifDirVec.at(g)>=dir) {
      ndir.push_back(unifDirVec.at(g));       
     }
  }
  vcgEmb::vcgForEmbree::embreeAmbientOcclusion(m3,ndir.size(),ndir,1);
  tri::UpdateQuality<MyMesh>::VertexFromFace(m3);
  tri::UpdateColor<MyMesh>::PerVertexQualityGray(m3);
  tri::io::ExporterOFF<MyMesh>::Save(m3,"testAODir.off",tri::io::Mask::IOM_VERTCOLOR);


  vcgEmb::vcgForEmbree::embreeSDF(m4,nOfRays);
  tri::UpdateQuality<MyMesh>::VertexFromFace(m4);
  tri::UpdateColor<MyMesh>::PerVertexQualityRamp(m4);
  tri::io::ExporterOFF<MyMesh>::Save(m4,"testSDF.off",tri::io::Mask::IOM_VERTCOLOR);

  return 0;
  
}