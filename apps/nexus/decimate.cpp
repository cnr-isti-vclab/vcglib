#include <vector>
#include <iostream>

// stuff to define the mesh
#include <vcg/simplex/vertex/with/afvmvn.h>
#include <vcg/math/quadric.h>
#include <vcg/complex/trimesh/base.h>
#include <vcg/simplex/face/with/av.h>

// io
//#include <wrap/io_trimesh/import_ply.h>
//#include <wrap/io_trimesh/export_ply.h>
// update
#include <vcg/complex/trimesh/update/topology.h>

#include <vcg/complex/local_optimization.h>
#include <vcg/complex/local_optimization/tri_edge_collapse_quadric.h>

#include <vcg/space/point3.h>

#include "pvoronoi.h"
#include "border.h"

class MyEdge;
class MyFace;
class MyVertex:public vcg::VertexAFVMVNf<DUMMYEDGETYPE , MyFace,DUMMYTETRATYPE>{public: 
ScalarType w;
vcg::math::Quadric<vcg::Plane3<ScalarType,false> >q;
ScalarType & W(){return w;}
} ;
class MyFace : public vcg::FaceAV<MyVertex,DUMMYEDGETYPE , MyFace>{};

class MyMesh: 
  public vcg::tri::TriMesh< std::vector<MyVertex>, std::vector<MyFace > >{};

class MyTriEdgeCollapse: 
  public vcg::tri::TriEdgeCollapseQuadric< MyMesh, MyTriEdgeCollapse >{
public:
  typedef  vcg::tri::TriEdgeCollapseQuadric<MyMesh, MyTriEdgeCollapse > TECQ;
  typedef  TECQ::PosType PosType;
  MyTriEdgeCollapse(PosType p, int i):TECQ(p,i){}
  ~MyTriEdgeCollapse(){}
};

using namespace vcg;
using namespace tri;
using namespace nxs;
using namespace std;

float Clustering(unsigned int target_faces, 
	       vector<Point3f> &newvert, 
	       vector<unsigned int> &newface,
	       vector<Link> &newbord,
	       vector<int> &vert_remap) {
}

float Cluster(MyMesh &mesh, unsigned int target_faces);

float Decimate(unsigned int target_faces, 
	       vector<Point3f> &newvert, 
	       vector<unsigned int> &newface,
	       vector<Link> &newbord,
	       vector<int> &vert_remap) {
  
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
      assert(newface[i+k] < mesh.vn);
      face.V(k) = &mesh.vert[newface[i+k]];
    }
    mesh.face.push_back(face);
  }
  mesh.fn = mesh.face.size();

  //emark borders 
  for(unsigned int i = 0; i < newbord.size(); i++) 
    mesh.vert[newbord[i].start_vert].ClearW();

  //  int FinalSize = mesh.face.size()/2;
  //  if(FinalSize > target_faces) FinalSize = target_faces;
  int FinalSize = target_faces;
  


  printf("mesh loaded %d %d \n",mesh.vn,mesh.fn);
  printf("reducing it to %i\n",FinalSize);


  /*  
      int t0=clock();	
      vcg::tri::UpdateTopology<MyMesh>::VertexFace(mesh);
      int t1=clock();	
      vcg::LocalOptimization<MyMesh> DeciSession(mesh);
      MyTriEdgeCollapse::SetDefaultParams();
      
      DeciSession.Init<MyTriEdgeCollapse>();
      
      FinalSize = mesh.fn - FinalSize; //number of faces to remove
      FinalSize/=2; //Number of vertices to remove
      DeciSession.SetTargetOperations(FinalSize);
      DeciSession.DoOptimization(); 
      float error = DeciSession.currMetric/4;//1; //get error; 
      int t3=clock();	
  */
  float error = Cluster(mesh, target_faces);


  
  /*  printf(" vol %d \n lkv %d \n lke %d \n lkf %d \n ood %d\n bor %d\n ",
	 MyTriEdgeCollapse::FailStat::Volume()           ,
	 MyTriEdgeCollapse::FailStat::LinkConditionFace(),
	 MyTriEdgeCollapse::FailStat::LinkConditionEdge(),
	 MyTriEdgeCollapse::FailStat::LinkConditionVert(),
	 MyTriEdgeCollapse::FailStat::OutOfDate()        ,
	 MyTriEdgeCollapse::FailStat::Border()           
	 );*/
 
  //  printf("Completed in %i+%i+%i msec\n",t1-t0,t2-t1,t3-t2);
  //  printf("mesh  %d %d \n",mesh.vn,mesh.fn);

  //recort vert start.

  
  newvert.clear();
  newface.clear();

  unsigned int totvert = 0;
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
    for(int k = 0; k < 3; k++) {
      assert(vert_remap[face.V(k) - vert_start] != -1);
      newface.push_back(vert_remap[face.V(k) - vert_start]);
    }
  }

  for(unsigned int i = 0; i < newbord.size(); i++) {
    unsigned short &v = newbord[i].start_vert;
    assert(vert_remap[v] != -1);
    v = vert_remap[v];
  }
  return error;
}


float Cluster(MyMesh &mesh, unsigned int target_faces) {
  unsigned int starting = mesh.vn;

  unsigned int nseeds = target_faces/2;
  assert(nseeds < mesh.vert.size());

  vector<unsigned int> remap;

  VoronoiPartition part;
  Box3f box;
  for(unsigned int i = 0; i < mesh.vert.size(); i++) {
    const Point3f &p = mesh.vert[i].cP();
    box.Add(p);
    if(!mesh.vert[i].IsW()) {
      part.push_back(Seed(p, 1));
      remap.push_back(i);
      nseeds--;
    }
  }
  unsigned int nborder = part.size();
  //Dovrei supersamplare prima....
  while(nseeds > 0) {
    unsigned int i = rand() % mesh.vert.size();
    if(mesh.vert[i].IsW() && !mesh.vert[i].IsV()) {
      const Point3f &p = mesh.vert[i].cP();
      part.push_back(Seed(p, 1));
      mesh.vert[i].SetV();
      remap.push_back(i);
      nseeds--;
    }
  }
  part.SetBox(box);
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
	part[i].p = centroid[i]/count[i];
    }
  }

  for(unsigned int i = nborder; i < part.size(); i++) {
    assert(mesh.vert[remap[i]].IsV());
    mesh.vert[remap[i]].P() = part[i].p;
  }

  float error = 0;
  //rimappiamo le facce.....
  for(unsigned int i = 0; i < mesh.face.size(); i++) {
    MyFace &face = mesh.face[i];
    for(int k = 0; k < 3; k++) {
      unsigned int target = part.Locate(face.V(k)->cP());
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

  cerr << "Error: " << error << endl;
  cerr << "faces: " << mesh.fn << endl;
  cerr << "verts: " << mesh.vn << endl;
  return error;
}

