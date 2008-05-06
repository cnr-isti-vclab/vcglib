#include<iostream>

#include<vcg/complex/trimesh/base.h>
#include<vcg/simplex/vertexplus/base.h>
#include<vcg/simplex/vertexplus/component.h>
#include<vcg/simplex/faceplus/base.h>
#include<vcg/simplex/faceplus/component.h>
#include<vcg/complex/trimesh/update/normal.h>
#include<vcg/complex/trimesh/create/platonic.h>

#include<vcg/space/color4.h>

class AFace;
class AVertex;
class AEdge;   // dummy prototype never used
class AVertex  : public vcg::VertexSimp2< AVertex, AEdge, AFace, vcg::vert::Coord3f>{};
class AFace    : public vcg::FaceSimp2< AVertex, AEdge, AFace, vcg::face::VertexRef,vcg::face::Color4b,vcg::face::Normal3f> {};

class AMesh : public vcg::tri::TriMesh< std::vector<AVertex>, std::vector<AFace> > {};

class CFace;
class CVertex;
class CEdge;   // dummy prototype never used
class CVertex  : public vcg::VertexSimp2< CVertex, CEdge, CFace, vcg::vert::Coord3f,vcg::vert::Normal3f>{};
class CFace    : public vcg::FaceSimp2< CVertex, CEdge, CFace, vcg::face::VertexRef, vcg::face::Normal3f> {};

class CMesh : public vcg::tri::TriMesh< std::vector<CVertex>, std::vector<CFace> > {};

int main(int , char **)
{
  AMesh am;
  CMesh cm;
	vcg::tri::Tetrahedron(cm);
	vcg::tri::Tetrahedron(am);
  
	std::cout << "Generated mesh has " << cm.vn << " vertices and " << cm.fn << " triangular faces" << std::endl;
  
  /// Calculates both vertex and face normals.
  /// The normal of a vertex v is the weigthed average of the normals of the faces incident on v.
  /// normals are not normalized
	vcg::tri::UpdateNormals<CMesh>::PerVertexPerFace(cm);
	std::cout << "[cm mesh] Normal of face 0 is [" << cm.face[0].N()[0] << "," << cm.face[0].N()[1] << ","  << cm.face[0].N()[2] << "]" << std::endl;
	std::cout << "[cm mesh] Normal of vertex 0 is [" << cm.vert[0].N()[0] << "," << cm.vert[0].N()[1] << ","  << cm.vert[0].N()[2] << "]" << std::endl;
  
	/// Calculates face normals.
  /// normals are not normalized
	vcg::tri::UpdateNormals<AMesh>::PerFace(am); 
	std::cout << "[am mesh] Normal of face 0 is [" << cm.face[0].N()[0] << "," << cm.face[0].N()[1] << ","  << cm.face[0].N()[2] << "]" << std::endl;

  return 0;
}
