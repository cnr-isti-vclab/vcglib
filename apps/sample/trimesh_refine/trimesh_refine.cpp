#include<vcg/simplex/vertex/base.h>
#include<vcg/simplex/face/base.h>
#include<vcg/simplex/face/component_ocf.h>
#include<vcg/simplex/face/topology.h>

#include<vcg/complex/trimesh/base.h>

#include <vcg/complex/trimesh/update/topology.h>
#include <vcg/complex/trimesh/update/normal.h>
#include <vcg/complex/trimesh/update/flag.h>
#include <vcg/complex/trimesh/refine.h>
#include <vcg/complex/trimesh/refine_loop.h>

#include <vcg/complex/trimesh/bitquad_creation.h>

// input output
#include <wrap/io_trimesh/import_ply.h>
#include <wrap/io_trimesh/export_ply.h>

// std
#include <vector>

using namespace vcg;
using namespace std;


class MyEdge;    // dummy prototype never used
class MyFace;
class MyVertex;

struct MyUsedTypes : public UsedTypes<	Use<MyVertex>::AsVertexType,
																				Use<MyFace>::AsFaceType>{};

class MyVertex  : public Vertex< MyUsedTypes, vertex::Coord3f, vertex::Normal3f, vertex::BitFlags  >{};
class MyFace    : public Face  < MyUsedTypes, face::InfoOcf, face::FFAdjOcf,  face::VertexRef, face::BitFlags > {};
class MyMesh    : public vcg::tri::TriMesh< vector<MyVertex>, face::vector_ocf<MyFace> > {};



#define FLAT	0
#define LOOP	1
#define CATMULL	2
#define BUTTERFLY 3

int  main(int argc, char **argv)
{
 if(argc<4)
	{
		printf(
		"\n                  PlyRefine ("__DATE__")\n"
			"						Visual Computing Group I.S.T.I. C.N.R.\n"
      "Usage: PlyRefine filein.ply fileout.ply ref_step [opt] \n"
			"Commands: \n"
			" Refinement rules:\n"
      "     -m  use simple midpoint subdivision (default) \n"
      "     -b  use butterfly subdivision scheme \n"
      "     -l  use butterfly subdivision scheme \n"
      "     -c  use butterfly subdivision scheme \n"
      "     -e# refine only if the the edge is longer than #(default 0.0)\n"
			);
		exit(0);
	}

	int RefMode = FLAT	;
  int i=4; int n_steps; float length=0;
	while(i<argc)
		{
			if(argv[i][0]!='-')
				{printf("Error unable to parse option '%s'\n",argv[i]); exit(0);}
			switch(argv[i][1])
			{				
				case 'm' :	RefMode=FLAT; break;
				case 'b' :	RefMode=BUTTERFLY; break;
      case 'l' :	RefMode=LOOP; break;
      case 'c' :	RefMode=CATMULL; break;
        case 'e' :	length=(float)atof(argv[i]+2); break;
				default : {printf("Error unable to parse option '%s'\n",argv[i]); exit(0);}
			}
			++i;
		}

	MyMesh m;
  if(tri::io::ImporterPLY<MyMesh>::Open(m,argv[1])!=0)
		{
      printf("Error reading file  %s\n",argv[1]);
			exit(0);
		}
  m.face.EnableFFAdjacency();
  tri::UpdateTopology<MyMesh>::FaceFace(m);
  tri::UpdateFlags<MyMesh>::FaceBorderFromFF(m);
  tri::UpdateNormals<MyMesh>::PerVertexNormalized(m);
  printf("Input mesh  vn:%i fn:%i\n",m.vn,m.fn);
	
  n_steps=atoi(argv[3]);
	
  for(i=0;i < n_steps;++i)			
  {
    switch(RefMode)
    {
    case FLAT:
      Refine<MyMesh, MidPoint<MyMesh> >(m,MidPoint<MyMesh>(&m),length);
      break;
    case LOOP:
      tri::RefineOddEven<MyMesh, tri::OddPointLoop<MyMesh>, tri::EvenPointLoop<MyMesh> >(m, tri::OddPointLoop<MyMesh>(), tri::EvenPointLoop<MyMesh>(), length);
      break;
    case CATMULL:
      tri::BitQuadCreation<MyMesh>::MakePureByRefine(m);
      assert(tri::BitQuadCreation<MyMesh>::IsBitTriQuadConventional(m));
      tri::UpdateTopology<MyMesh>::FaceFace(m);
      tri::BitQuadCreation<MyMesh>::MakePureByRefine(m);
      tri::UpdateNormals<MyMesh>::PerBitQuadFaceNormalized(m);
      break;
    case BUTTERFLY:
      Refine<MyMesh, MidPointButterfly<MyMesh> >(m,MidPointButterfly<MyMesh>(),length);
      break;
    }					
  }
  
  printf("Output mesh vn:%i fn:%i\n",m.vn,m.fn);
	
  vcg::tri::io::PlyInfo pi;
	vcg::tri::io::ExporterPLY<MyMesh>::Save(m,argv[2],pi.mask);
	return 0;
	}
