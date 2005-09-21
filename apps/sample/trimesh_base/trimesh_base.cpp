#include <stdio.h>
#include<vcg/simplex/vertex/vertex.h>
#include<vcg/simplex/vertex/with/vn.h>
#include<vcg/simplex/face/with/fn.h>
#include<vcg/simplex/face/with/fcfn.h>
#include<vcg/complex/trimesh/base.h>
#include<vcg/complex/trimesh/create/platonic.h>
#include<vcg/complex/trimesh/update/normal.h>

using namespace vcg;
using namespace std;

class AEdge;    // dummy prototype never used
class AFace;
class AVertex   : public vcg::Vertex< double,AEdge,AFace > {};
class AFace     : public vcg::FaceFCFN< AVertex,AEdge,AFace > {};
class AMesh     : public vcg::tri::TriMesh< std::vector<AVertex>, std::vector<AFace> > {};

class CEdge;    // dummy prototype never used
class CFace;
class CVertex   : public vcg::VertexVN< double,CEdge,CFace > {};
class CFace     : public vcg::FaceFN< CVertex,CEdge,CFace > {};
class CMesh     : public vcg::tri::TriMesh< std::vector<CVertex>, std::vector<CFace> > {};


int main(int , char **)
{
  AMesh am;
  CMesh cm;
  tri::Tetrahedron(cm);
  tri::Tetrahedron(am);
  
  printf("Generated mesh has %i vertices and %i triangular faces\n",cm.vn,cm.fn);
  
  /// Calculates both vertex and face normals.
  /// The normal of a vertex v is the weigthed average of the normals of the faces incident on v.
  /// normals are not normalized
  tri::UpdateNormals<CMesh>::PerVertexPerFace(cm);
  printf("Normal of face 0 is %f %f %f",cm.face[0].N()[0],cm.face[0].N()[1],cm.face[0].N()[2]);
  
  tri::UpdateNormals<AMesh>::PerFace(am);  
  
  return 0;
}
