#include <iostream>

// stuff to define the mesh
#include <vcg/simplex/vertex/with/afvmvn.h>
#include <vcg/simplex/edge/edge.h>
#include <vcg/math/quadric.h>
#include <vcg/complex/trimesh/base.h>
#include <vcg/simplex/face/with/av.h>

#include <vcg/complex/trimesh/update/topology.h>

#include <vcg/complex/local_optimization.h>
#include <vcg/complex/local_optimization/tri_edge_collapse_quadric.h>

#include <vcg/space/point3.h>

#include "vpartition.h"
#include "fragment.h"

#include "decimate.h"
#include <wrap/io_trimesh/export_ply.h>

using namespace vcg;
using namespace tri;
using namespace nxs;
using namespace std;

class MyEdge;
class MyFace;
class MyVertex: 
  public vcg::VertexAFVMVNf<MyEdge, MyFace,DUMMYTETRATYPE> {
public: 
  ScalarType w;
  vcg::math::Quadric<double> q;
  ScalarType & W() { return w; }
};

struct MyEdge: public Edge<double,MyEdge,MyVertex> {
  inline MyEdge():Edge<double,MyEdge,MyVertex>(){UberFlags()=0;}
  inline MyEdge(MyVertex* a,MyVertex* b):Edge<double,MyEdge,MyVertex>(a,b){
    UberFlags()=0;}
};
class MyFace : public vcg::FaceAV<MyVertex, MyEdge, MyFace> {};

class MyMesh: 
  public vcg::tri::TriMesh< std::vector<MyVertex>, std::vector<MyFace > > {};

class MyTriEdgeCollapse: public vcg::tri::TriEdgeCollapseQuadric< MyMesh, MyTriEdgeCollapse > {
						public:
						typedef  vcg::tri::TriEdgeCollapseQuadric< MyMesh,  MyTriEdgeCollapse > TECQ;
						typedef  TECQ::EdgeType EdgeType;
						inline MyTriEdgeCollapse(  EdgeType p, int i) :TECQ(p,i){}
};

float Cluster(MyMesh &mesh, unsigned int target_faces);
float Quadric(MyMesh &mesh, unsigned int target_faces);

float nxs::Decimate(Decimation mode,
		    unsigned int target_faces, 
		    vector<Point3f> &newvert, 
		    vector<unsigned int> &newface,
		    vector<BigLink> &newbord) {          

  for(unsigned int i = 0; i < newface.size(); i+= 3) {
    assert(newface[i] != newface[i+1]);
    assert(newface[i] != newface[i+2]);
    assert(newface[i+1] != newface[i+2]);
  }
  
  MyMesh mesh;
  
  //build mesh
  for(unsigned int i = 0; i < newvert.size(); i++) {
    MyVertex vertex;
    vertex.ClearFlags();
    vertex.P() = newvert[i];
    mesh.vert.push_back(vertex);
  }
  mesh.vn = mesh.vert.size();

  for(unsigned int i = 0; i < newface.size(); i+=3) {
    MyFace face;
    face.ClearFlags();
    for(int k = 0; k < 3; k++) {
      assert(newface[i+k] < mesh.vert.size());
      face.V(k) = &mesh.vert[newface[i+k]];
    }
    mesh.face.push_back(face);
  }
  mesh.fn = mesh.face.size();

  //mark borders 
  for(unsigned int i = 0; i < newbord.size(); i++) 
    mesh.vert[newbord[i].start_vert].ClearW();

  //    vcg::tri::io::ExporterPLY<MyMesh>::Save(mesh, "ribum.ply");

  float error;
  switch(mode) {
    case CLUSTER: error = Cluster(mesh, target_faces); break;
    case QUADRIC: error = Quadric(mesh, target_faces); break;
    default: cerr << "Unknown simplification mode: " << mode << endl;
             exit(0);
  }

  newvert.clear();
  newface.clear();

  unsigned int totvert = 0;
  vector<int> vert_remap;
  vert_remap.resize(mesh.vert.size(), -1);
  for(unsigned int i = 0; i < mesh.vert.size(); i++) {
    if(mesh.vert[i].IsD()) continue;
    newvert.push_back(mesh.vert[i].cP());
    vert_remap[i] = totvert++;
  }

  MyMesh::VertexPointer vert_start = &mesh.vert[0];
  for(unsigned int i = 0; i < mesh.face.size(); i++) {
    MyFace &face = mesh.face[i];
    if(face.IsD()) continue;
    for(int k = 0; k < 3; k++) 
      newface.push_back(vert_remap[face.V(k) - vert_start]);
  }

  for(unsigned int i = 0; i < newbord.size(); i++) {
    unsigned int &v = newbord[i].start_vert;
    assert(vert_remap[v] != -1);
    v = vert_remap[v];
  }

  //Temporary test again:
  /*  for(unsigned int i = 0; i < newface.size(); i+= 3) {
    assert(newface[i] != newface[i+1]);
    assert(newface[i] != newface[i+2]);
    assert(newface[i+1] != newface[i+2]);
    }*/
  
  return error;
}


float Quadric(MyMesh &mesh, unsigned int target_faces) {
  vcg::tri::UpdateTopology<MyMesh>::VertexFace(mesh);
  vcg::tri::UpdateBounding<MyMesh>::Box(mesh);

  vcg::LocalOptimization<MyMesh> DeciSession(mesh);
  
  MyTriEdgeCollapse::SetDefaultParams();
      
  DeciSession.Init<MyTriEdgeCollapse>();
      
  DeciSession.SetTargetSimplices(target_faces);
  DeciSession.DoOptimization(); 
  
  float error = 0;
  int count = 0;
  for(unsigned int i = 0; i < mesh.face.size(); i++) {
    MyFace &face = mesh.face[i];
    if(face.IsD()) continue;
    for(int k = 0; k < 3; k++) {
      error += (face.cV(k)->cP() - face.cV((k+1)%3)->cP()).Norm();
      count++;
    }
  }
  error /= count;
  return error;     
  return 0;
}

float Cluster(MyMesh &mesh, unsigned int target_faces) {
  unsigned int starting = mesh.vn;
  
  unsigned int nseeds = target_faces/2;
#ifndef NDEBUG
  if(nseeds >= mesh.vert.size()) {
    cerr << "Strange! nseeds > vert.size(): " << nseeds  << " >= "<< mesh.vert.size() << endl;
  }
#endif
  
  vector<unsigned int> remap;
  
  VPartition part;
  for(unsigned int i = 0; i < mesh.vert.size(); i++) {
    const Point3f &p = mesh.vert[i].cP();
    if(!mesh.vert[i].IsW()) {
      part.push_back(p);
      remap.push_back(i);
      nseeds--;
    }
  }
  unsigned int nborder = part.size();
  //Dovrei supersamplare prima....
  while(nseeds > 0 && part.size() < mesh.vn) {
    unsigned int i = rand() % mesh.vert.size();
    if(mesh.vert[i].IsW() && !mesh.vert[i].IsV()) {
      const Point3f &p = mesh.vert[i].cP();
      part.push_back(p);
      mesh.vert[i].SetV();
      remap.push_back(i);
      nseeds--;
    }
  }
  part.Init();

  vector<Point3f> centroid;
  vector<unsigned int> count;
  for(unsigned int i = 0; i < 3; i++) {
    centroid.clear();
    centroid.resize(mesh.vert.size(), Point3f(0, 0, 0));
    count.clear();
    count.resize(mesh.vert.size(), 0);
    for(unsigned int i = 0; i < mesh.vert.size(); i++) {
      unsigned int target = part.Locate(mesh.vert[i].cP());
      centroid[target] += mesh.vert[i].cP();
      count[target]++;
    }
    for(unsigned int i = nborder; i < part.size(); i++) {
      if(count[i] > 0)
      	part[i] = centroid[i]/count[i];
    }
  }

  for(unsigned int i = nborder; i < part.size(); i++) {
    assert(mesh.vert[remap[i]].IsV());
    mesh.vert[remap[i]].P() = part[i];
  }

  float error = 0;
  //rimappiamo le facce.....
  for(unsigned int i = 0; i < mesh.face.size(); i++) {
    MyFace &face = mesh.face[i];
    for(int k = 0; k < 3; k++) {
      unsigned int target = part.Locate(face.V(k)->cP());
      assert(target < remap.size());
      assert(remap[target] < mesh.vert.size());
      MyVertex &vert = mesh.vert[remap[target]];

      float dist = Distance(vert.cP(), face.V(k)->cP());
      if(dist > error) error = dist;

      face.V(k) = &vert;
    }
  }
  
  for(unsigned int i = 0; i < mesh.face.size(); i++) {
    MyFace &face = mesh.face[i];
    assert(!face.IsD());
    for(int k = 0; k < 3; k++) {
      assert(face.cV(k)->IsV() || !face.cV(k)->IsW());
    }
    if(face.cV(0) == face.cV(1) ||
       face.cV(0) == face.cV(2) ||
       face.cV(1) == face.cV(2)) {
      face.SetD();
      mesh.fn--;
    }
  }
  
  for(unsigned int i = 0; i < mesh.vert.size(); i++)
    if(!mesh.vert[i].IsV() && mesh.vert[i].IsW()) {
      mesh.vert[i].SetD();
      mesh.vn--;
    }
  return error;

}

