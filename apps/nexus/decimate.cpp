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

  int t0=clock();	

  printf("mesh loaded %d %d \n",mesh.vn,mesh.fn);
  printf("reducing it to %i\n",FinalSize);
  vcg::tri::UpdateTopology<MyMesh>::VertexFace(mesh);

  //  cerr << "topology ok" << endl;
  int t1=clock();	
  
  //micro random semplificatore

  /*  for(unsigned int i = 0; i < mesh.face.size(); i+= 2) {
    MyFace &face = mesh.face[i];
    if(face.V(0)->IsW() && face.V(1)->IsW() && face.V(1)->IsW())
      mesh.face[i].SetD();
  }

  for(unsigned int i = 0; i < mesh.vert.size(); i++) {
    if(mesh.vert[i].IsW())
      mesh.vert[i].SetD();
  }

  for(unsigned int i = 0; i < mesh.face.size(); i++) {
    MyFace &face = mesh.face[i];
    if(face.IsD()) continue;
    face.V(0)->ClearD();
    face.V(1)->ClearD();
    face.V(2)->ClearD();
    }*/

  //  Cluster(mesh, target_faces);
  //  cerr << "simplified" << endl;

  vcg::LocalOptimization<MyMesh> DeciSession(mesh);
  MyTriEdgeCollapse::SetDefaultParams();

  DeciSession.Init<MyTriEdgeCollapse>();

  int t2=clock();	
  //  printf("Initial Heap Size %i\n",DeciSession.h.size());
  
  FinalSize = mesh.fn - FinalSize; //number of faces to remove
  FinalSize/=2; //Number of vertices to remove
  DeciSession.SetTargetOperations(FinalSize);
  //  DeciSession.DoOptimization(); 
  float error = Cluster(mesh, target_faces);
  //  float error = DeciSession.currMetric/4;//1; //get error;
  int t3=clock();	
  
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
  cerr << "starting face: " << mesh.fn << endl;
  //veramente brutale
  vector<int> remap;
  remap.resize(mesh.vert.size());
  for(int i = 0; i < mesh.vert.size(); i++)
    remap[i] = -1;

  int toremove = mesh.fn - target_faces;

  cerr << "counting" << endl;
  map<float, pair<int, int> > dist;
  for(int i = 0; i < mesh.vert.size(); i++) {
    if(mesh.vert[i].IsD()) continue;
    if(!mesh.vert[i].IsW()) continue;
    for(int k = i+1; k < mesh.vert.size(); k++) {
      if(mesh.vert[k].IsD()) continue;
      if(!mesh.vert[k].IsW()) continue;
      float d = (mesh.vert[i].P() - mesh.vert[k].P()).SquaredNorm();
      dist[d] = make_pair(i, k);
    }
  }
  
  float error = 0;
  cerr << "done" << endl;
  map<float, pair<int, int> >::iterator s;
  for(s = dist.begin(); s != dist.end(); s++) {
    if(toremove < 0) break;
    int target = (*s).second.first;
    int source = (*s).second.second;

    if(remap[target] != -1) continue;
    if(remap[source] != -1) continue;

    assert(!mesh.vert[target].IsD());
    assert(!mesh.vert[source].IsD());

    mesh.vert[source].SetD();
    error = (*s).first;
    remap[source] = target;
    remap[target] = target;
    toremove -= 2;
    mesh.vn--;
    //    if(mesh.vn < starting/2) break;
  }

  //PULIAMO LE FACCE
  for(int i = 0; i < mesh.face.size(); i++) {
    MyFace &face = mesh.face[i];
    if(face.IsD()) continue;
    for(int k = 0; k < 3; k++) {
      if(face.V(k)->IsD()) {
	face.V(k) = &mesh.vert[remap[face.V(k) - &mesh.vert[0]]];
      }
      assert(!face.V(k)->IsD());
    }
    if(face.V(0) == face.V(1) || face.V(0) == face.V(2) || 
       face.V(1) == face.V(2)) {
      face.SetD();
      mesh.fn--;
    }
  }
  cerr << "Ending faces: " << mesh.fn << endl;
  return error;
}
