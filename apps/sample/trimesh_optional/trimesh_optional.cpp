#include <stdio.h>
#include <time.h>

#include<vcg/simplex/vertexplus/base.h>
#include<vcg/simplex/faceplus/base.h>
#include<vcg/simplex/face/topology.h>

#include<vcg/simplex/vertexplus/component_ocf.h>
#include<vcg/simplex/faceplus/component_ocf.h>

#include<vcg/simplex/vertexplus/component_occ.h>
#include<vcg/simplex/faceplus/component_occ.h>

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
class CVertex;
class CVertexOcf;

// Opt stuff

// OCF Optional Component Fast 
class CVertex     : public VertexSimp2< CVertex,    CEdge, CFace,    vert::Coord3f, vert::Normal3f >{};
class CVertexOcf  : public VertexSimp2< CVertexOcf, CEdge, CFaceOcf, vert::Coord3f, vert::Normal3f >{};

class CFace       : public FaceSimp2< CVertex,    CEdge, CFace,                   face::FFAdj,    face::VertexRef, face::Flag, face::Normal3f > {};
class CFaceOcf    : public FaceSimp2< CVertexOcf, CEdge, CFaceOcf, face::InfoOcf, face::FFAdjOcf, face::VertexRef, face::Flag, face::Normal3fOcf > {};

class CMeshOcf    : public vcg::tri::TriMesh< vector<CVertexOcf>, face::vector_ocf<CFaceOcf> > {};
class CMesh       : public vcg::tri::TriMesh< vector<CVertex   >,           vector<CFace   > > {};

// OCC Optional Component Compact
class CFaceOcc;
class CVertexOcc   : public VertexSimp2< CVertexOcc,CEdge, CFaceOcc,vert::Coord3f,vert::Normal3f >{};
class CFaceOcc     : public FaceSimp2<  CVertexOcc, CEdge, CFaceOcc,face::VertexRef,vcg::face::Normal3fOcc,face::FFAdjOcc,face::Flag> {};
class CMeshOcc     : public vcg::tri::TriMesh< vector_occ<CVertexOcc   >, vector_occ<CFaceOcc  > > {};


int main(int , char **)
{
  CMesh cm;
  CMeshOcf cmo;
	CMeshOcc cmoc;


  tri::Tetrahedron(cm);
  tri::Tetrahedron(cmo);
 	tri::Tetrahedron(cmoc);
 
  printf("Generated mesh has %i vertices and %i triangular faces\n",cm.vn,cm.fn);
  
  /// Calculates both vertex and face normals.
  /// The normal of a vertex v is the weigthed average of the normals of the faces incident on v.
  /// normals are not normalized

	cmoc.face.EnableAttribute<CFaceOcc::FFAdjType>();  
	cmo.face.EnableFFAdjacency();


  printf("Size of CFace            %3i\n",sizeof(CFace));
  printf("Size of CFaceOcf         %3i\n",sizeof(CFaceOcf));
  printf("Size of CFaceOcc         %3i\n",sizeof(CFaceOcc));
 
  vcg::tri::UpdateTopology<CMesh   >::FaceFace(cm);
  vcg::tri::UpdateTopology<CMeshOcf>::FaceFace(cmo);
  vcg::tri::UpdateTopology<CMeshOcc>::FaceFace(cmoc);
	
  vcg::tri::UpdateFlags<CMesh   >::FaceBorderFromFF(cm);
	vcg::tri::UpdateFlags<CMeshOcf>::FaceBorderFromFF(cmo);
 	vcg::tri::UpdateFlags<CMeshOcc>::FaceBorderFromFF(cmoc);

	vcg::tri::UpdateNormals<CMesh   >::PerVertexNormalized(cm);
	vcg::tri::UpdateNormals<CMeshOcf>::PerVertexNormalized(cmo);
	vcg::tri::UpdateNormals<CMeshOcc>::PerVertexNormalized(cmoc);


  printf("Normal of face 0 is %f %f %f\n\n",cm.face[0].N()[0],cm.face[0].N()[1],cm.face[0].N()[2]);
  int t0=0,t1=0,t2=0,t3=0;
  while(t3-t0<3000)
  {
    t0=clock();
    Refine(cm,MidPointButterfly<CMesh>(),0); 
    t1=clock();
    Refine(cmo,MidPointButterfly<CMeshOcf>(),0); 
    t2=clock();
    Refine(cmoc,MidPointButterfly<CMeshOcc>(),0);	
		t3=clock();
		printf("Mesh is %i %i in Std:%i Ocf:%i Occ:%i\n",cm.vn,cm.fn,t1-t0,t2-t1,t3-t2);
  }
  return 0;
}
