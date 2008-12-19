#include <stdio.h>
#include <time.h>

#include<vcg/simplex/vertex/base.h>
#include<vcg/simplex/face/base.h>
#include<vcg/simplex/face/topology.h>

#include<vcg/simplex/vertex/component_ocf.h>
#include<vcg/simplex/face/component_ocf.h>

#include<vcg/simplex/vertex/component_occ.h>
#include<vcg/simplex/face/component_occ.h>

#include<vcg/complex/trimesh/base.h>
#include<vcg/complex/trimesh/create/platonic.h>
#include<vcg/complex/trimesh/update/topology.h>
#include<vcg/complex/trimesh/update/flag.h>
#include<vcg/complex/trimesh/update/normal.h>
#include <vcg/complex/trimesh/refine.h>

using namespace vcg;
using namespace std;

class CEdge;    // dummy prototype never used
class CFace;
class CFaceOcf;
class CFaceOcc;
class CVertex;
class CVertexOcf;

// Optional stuff has two suffixes:
// OCF Optional Component Fast 
// OCC Optional Component Compact

class CVertex     : public VertexSimp2< CVertex,    CEdge, CFace,    vertex::Coord3f, vertex::BitFlags,vertex::Normal3f >{};
class CVertexOcf  : public VertexSimp2< CVertexOcf, CEdge, CFaceOcf, vertex::Coord3f, vertex::BitFlags,vertex::Normal3f >{};
class CVertexOcc  : public VertexSimp2< CVertexOcc, CEdge, CFaceOcc, vertex::Coord3f, vertex::BitFlags,vertex::Normal3f >{};

class CFace       : public FaceSimp2< CVertex,    CEdge, CFace,                   face::FFAdj,    face::VertexRef, face::BitFlags, face::Normal3f > {};
class CFaceOcf    : public FaceSimp2< CVertexOcf, CEdge, CFaceOcf, face::InfoOcf, face::FFAdjOcf, face::VertexRef, face::BitFlags, face::Normal3fOcf > {};
class CFaceOcc    : public FaceSimp2< CVertexOcc, CEdge, CFaceOcc,                face::FFAdjOcc, face::VertexRef, face::BitFlags, face::Normal3fOcc > {};

class CMesh       : public vcg::tri::TriMesh<     vector<CVertex   >,           vector<CFace   > > {};
class CMeshOcf    : public vcg::tri::TriMesh<     vector<CVertexOcf>, face::vector_ocf<CFaceOcf> > {};
class CMeshOcc    : public vcg::tri::TriMesh< vector_occ<CVertexOcc>,       vector_occ<CFaceOcc  > > {};



int main(int , char **)
{
  CMesh cm;
  CMeshOcf cmof;
 	CMeshOcc cmoc;


  tri::Tetrahedron(cm);
  tri::Tetrahedron(cmof);
 	tri::Tetrahedron(cmoc);
 
  printf("Generated mesh has %i vertices and %i triangular faces\n",cm.vn,cm.fn);
  
  /// Calculates both vertex and face normals.
  /// The normal of a vertex v is the weigthed average of the normals of the faces incident on v.
  /// normals are not normalized

	cmoc.face.EnableAttribute<CFaceOcc::FFAdjType>();  
	cmof.face.EnableFFAdjacency();


  printf("Size of CFace            %3i\n",sizeof(CFace));
  printf("Size of CFaceOcf         %3i\n",sizeof(CFaceOcf));
	printf("Size of CFaceOcc         %3i\n",sizeof(CFaceOcc));
 
  vcg::tri::UpdateTopology<CMesh   >::FaceFace(cm);
  vcg::tri::UpdateTopology<CMeshOcf>::FaceFace(cmof);
	vcg::tri::UpdateTopology<CMeshOcc>::FaceFace(cmoc);
	
  vcg::tri::UpdateFlags<CMesh   >::FaceBorderFromFF(cm);
	vcg::tri::UpdateFlags<CMeshOcf>::FaceBorderFromFF(cmof);
	vcg::tri::UpdateFlags<CMeshOcc>::FaceBorderFromFF(cmoc);

	vcg::tri::UpdateNormals<CMesh   >::PerVertexNormalized(cm);
	vcg::tri::UpdateNormals<CMeshOcf>::PerVertexNormalized(cmof);
	vcg::tri::UpdateNormals<CMeshOcc>::PerVertexNormalized(cmoc);


  printf("Normal of face 0 is %f %f %f\n\n",cm.face[0].N()[0],cm.face[0].N()[1],cm.face[0].N()[2]);
  int t0=0,t1=0,t2=0,t3=0;
  while(t3-t0<3000)
  {
    t0=clock();
    Refine(cm,MidPointButterfly<CMesh>(),0); 
    t1=clock();
    Refine(cmof,MidPointButterfly<CMeshOcf>(),0); 
    t2=clock();
		Refine(cmoc,MidPointButterfly<CMeshOcc>(),0);	
		t3=clock();
		printf("Mesh is %i %i in Std:%i Ocf:%i Occ:%i\n",cm.vn,cm.fn,t1-t0,t2-t1,t3-t2);
  }
  return 0;
}
