#include <stdio.h>
#include <time.h>

#include "mesh_definition.h"
#include<vcg/complex/trimesh/allocate.h>
#include<vcg/complex/trimesh/create/platonic.h>
#include<vcg/complex/trimesh/update/topology.h>
#include<vcg/complex/trimesh/update/flag.h>
#include<vcg/complex/trimesh/update/normal.h>
#include <vcg/complex/trimesh/refine.h>

using namespace vcg;
using namespace std;

int main(int , char **)
{

	vcg::tri::Allocator<CMesh>::NameTypeScope bounds;
	vcg::tri::Allocator<CMesh>::AddNameTypeBound<float>(bounds,"myfloat");


	CMesh cm;
	CMeshOcf cmof;
	CMeshOcc cmoc;

	CMesh::VertexPointer v = cm.face[0].V(0);


	cmoc.face.EnableAttribute<CFaceOcc::NormalType>();
	CMeshOcc::FaceIterator fi = vcg::tri::Allocator<CMeshOcc>::AddFaces(cmoc,1);
	(*fi).N() = vcg::Point3f(9,9,9);


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
