#include <stdio.h>
#include <time.h>

#include "mesh_definition.h"
#include<vcg/complex/allocate.h>
#include<vcg/complex/algorithms/create/platonic.h>
#include<vcg/complex/algorithms/update/topology.h>
#include<vcg/complex/algorithms/update/flag.h>
#include<vcg/complex/algorithms/update/normal.h>
#include<vcg/complex/algorithms/update/bounding.h>
#include <vcg/complex/algorithms/refine.h>

using namespace vcg;
using namespace std;

int main(int , char **)
{



	CMesh cm;
	CMeshOcf cmof;
 


  tri::Tetrahedron(cm);
  tri::Tetrahedron(cmof);
 
  printf("Generated mesh has %i vertices and %i triangular faces\n",cm.vn,cm.fn);
  
  /// Calculates both vertex and face normals.
  /// The normal of a vertex v is the weigthed average of the normals of the faces incident on v.
  /// normals are not normalized
  cmof.face.EnableFFAdjacency();

  tri::UpdateTopology<CMesh   >::FaceFace(cm);
  tri::UpdateTopology<CMeshOcf>::FaceFace(cmof);

  tri::UpdateFlags<CMesh   >::FaceBorderFromFF(cm);
  tri::UpdateFlags<CMeshOcf>::FaceBorderFromFF(cmof);

  tri::UpdateNormal<CMesh   >::PerVertexNormalized(cm);
  tri::UpdateNormal<CMeshOcf>::PerVertexNormalized(cmof);

  printf("Normal of face 0 is %f %f %f\n\n",cm.face[0].N()[0],cm.face[0].N()[1],cm.face[0].N()[2]);
  int t0=0,t1=0;
  while(t1-t0<200)
  {
    t0=clock();
    tri::Refine(cm,tri::MidPointButterfly<CMesh>(cm),0);
    t1=clock();
    tri::Refine(cmof,tri::MidPointButterfly<CMeshOcf>(cmof),0);
  }

	cmof.vert.EnableRadius();
	cmof.vert.EnableQuality();

	unsigned int hh = 0;
	for(CMeshOcf::VertexIterator vi = cmof.vert.begin(); vi!=cmof.vert.end();++vi,++hh){
		if(hh%3==0)
			vcg::tri::Allocator<CMeshOcf>::DeleteVertex(cmof,*vi);
	}

  for(CMeshOcf::VertexIterator vi = cmof.vert.begin(); vi!=cmof.vert.end();++vi)
      {
        if(!(*vi).IsD())
        {
          float q =vi->Q();
          float r =vi->R();
          assert(q==r);
        }
      }

      tri::Allocator<CMeshOcf>::CompactVertexVector(cmof);
      tri::UpdateBounding<CMeshOcf>::Box(cmof);
      for(CMeshOcf::VertexIterator vi = cmof.vert.begin(); vi!=cmof.vert.end();++vi)
      {
        if(!(*vi).IsD())
        {
          float q =vi->Q();
          float r =vi->R();
//          int ii = vcg::tri::Index(cmof, *vi);
          assert(q==r);
        }
      }


  return 0;
}
