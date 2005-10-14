#include <stdio.h>
#include <time.h>

#include<vcg/simplex/vertexplus/base.h>
#include<vcg/simplex/faceplus/base.h>
#include<vcg/simplex/face/topology.h>

#include<vcg/simplex/vertexplus/component_ocf.h>
#include<vcg/simplex/faceplus/component_ocf.h>

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

class CVertex     : public VertexSimp2< CVertex,    CEdge, CFace,    vert::Coord3f, vert::Normal3f >{};
class CVertexOcf  : public VertexSimp2< CVertexOcf, CEdge, CFaceOcf, vert::Coord3f, vert::Normal3f >{};

class CFace       : public FaceSimp2< CVertex,    CEdge, CFace,                   face::FFAdj,    face::VertexRef, face::Flag, face::Normal3f > {};
class CFaceOcf    : public FaceSimp2< CVertexOcf, CEdge, CFaceOcf, face::InfoOcf, face::FFAdjOcf, face::VertexRef, face::Flag, face::Normal3fOcf > {};

class CMeshOcf    : public vcg::tri::TriMesh< vector<CVertexOcf>, face::vector_ocf<CFaceOcf> > {};
class CMesh       : public vcg::tri::TriMesh< vector<CVertex   >,           vector<CFace   > > {};


int main(int , char **)
{
  CMesh cm;
  CMeshOcf cmo;
  tri::Tetrahedron(cm);
  tri::Tetrahedron(cmo);
  
  printf("Generated mesh has %i vertices and %i triangular faces\n",cm.vn,cm.fn);
  
  /// Calculates both vertex and face normals.
  /// The normal of a vertex v is the weigthed average of the normals of the faces incident on v.
  /// normals are not normalized

 

  cmo.face.EnableFFAdjacency();
 
  vcg::tri::UpdateTopology<CMesh   >::FaceFace(cm);
  vcg::tri::UpdateTopology<CMeshOcf>::FaceFace(cmo);
	
  vcg::tri::UpdateFlags<CMesh   >::FaceBorderFromFF(cm);
	vcg::tri::UpdateFlags<CMeshOcf>::FaceBorderFromFF(cmo);

	vcg::tri::UpdateNormals<CMesh   >::PerVertexNormalized(cm);
	vcg::tri::UpdateNormals<CMeshOcf>::PerVertexNormalized(cmo);

  printf("Size of CFace            %3i\n",sizeof(CFace));
  printf("Size of CFaceOcf         %3i\n",sizeof(CFaceOcf));

  printf("Normal of face 0 is %f %f %f\n\n",cm.face[0].N()[0],cm.face[0].N()[1],cm.face[0].N()[2]);
  int t0=0,t1=0,t2=0;
  while(t2-t0<3000)
  {
    t0=clock();
    Refine(cm,MidPointButterfly<CMesh>(),0); 
    t1=clock();
    Refine(cmo,MidPointButterfly<CMeshOcf>(),0); 
    t2=clock();
    printf("Mesh is %i %i in Std:%i Ocf:%i\n",cm.vn,cm.fn,t1-t0,t2-t1);
  }
  return 0;
}
