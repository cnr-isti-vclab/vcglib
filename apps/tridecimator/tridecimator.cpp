#include <vector>
#include <limits>

#include <stdio.h>
#include <stdlib.h>

using namespace std;

// stuff to define the mesh
#include <vcg/simplex/vertex/with/afvmvn.h>
#include <vcg/simplex/edge/edge.h>
#include <vcg/math/quadric.h>
#include <vcg/complex/trimesh/base.h>
#include <vcg/complex/trimesh/clean.h>
#include <vcg/simplex/face/with/av.h>

// io
#include <wrap/io_trimesh/import.h>
#include <wrap/io_trimesh/export_ply.h>

// update
#include <vcg/complex/trimesh/update/topology.h>
// local optimization
#include <vcg/complex/local_optimization.h>
#include <vcg/complex/local_optimization/tri_edge_collapse_quadric.h>

using namespace vcg;
using namespace tri;

class MyEdge;
class MyFace;
class MyVertex;

// for edge collpases we need verteses with:
// AF V->F adjacency
// VM per vertex incremental mark
// VN per vertex Normal
// Moreover for using this vertex also in a quadric based collapse it must have also a Quadric member Q();
class MyVertex:public vcg::VertexAFVMVNd<MyEdge , MyFace,DUMMYTETRATYPE>{
public: 
  ScalarType w;
  vcg::math::Quadric<double>q;
  ScalarType & W(){return w;}
};

class MyEdge : public Edge<MyEdge,MyVertex> {};
class MyFace : public vcg::FaceAV<MyVertex,Edge<MyEdge,MyVertex> , MyFace>{};
class MyMesh : public vcg::tri::TriMesh< std::vector<MyVertex>, std::vector<MyFace > >{};

class MyTriEdgeCollapse: public vcg::tri::TriEdgeCollapseQuadric< MyMesh, MyTriEdgeCollapse > {
						public:
						typedef  vcg::tri::TriEdgeCollapseQuadric< MyMesh,  MyTriEdgeCollapse > TECQ;
						typedef  TECQ::EdgeType EdgeType;
						inline MyTriEdgeCollapse(  const EdgeType &p, int i) :TECQ(p,i){}
};


// mesh to simplify
MyMesh mesh;

void Usage()
{
  	printf(
           "---------------------------------\n"
           "         TriSimp V.1.0 \n"
           "     http://vcg.isti.cnr.it\n"
           "    http://vcg.sourceforge.net\n"
           "   release date: "__DATE__"\n"
           "---------------------------------\n\n"
		  "TriDecimator 1.0 \n"__DATE__"\n"
			"Copyright 2003-2006 Visual Computing Lab I.S.T.I. C.N.R.\n"
      "\nUsage:  "\
      "tridecimator file1 file2 face_num [opt]\n"\
      "Where opt can be:\n"\
      "     -e# QuadricError threshold  (range [0,inf) default inf)\n"
			"     -b# Boundary Weight (default .5)\n"
			"     -q# Quality threshold (range [0.0, 0.866],  default .3 )\n"
			"     -n# Normal threshold  (degree range [0,180] default 90)\n"
			"     -E# Minimal admitted quadric value (default 1e-15, must be >0)\n"
			"     -Q[y|n]  Use or not Quality Threshold (default yes)\n"
			"     -N[y|n]  Use or not Normal Threshold (default no)\n"
			"     -A[y|n]  Use or not Area Weighted Quadrics (default yes)\n"
			"     -O[y|n]  Use or not vertex optimal placement (default yes)\n"
			"     -S[y|n]  Use or not Scale Independent quadric measure(default yes) \n"
			"     -B[y|n]  Preserve or not mesh boundary (default no)\n"
			"     -T[y|n]  Preserve or not Topology (default no)\n"
			"     -H[y|n]  Use or not Safe Heap Update (default no)\n"
		  "     -P       Before simplification, remove duplicate & unreferenced vertices\n"
                                       );
  exit(-1);
}


int main(int argc ,char**argv){
if(argc<4) Usage();

  int FinalSize=atoi(argv[3]);
	//int t0=clock();	
  int err=vcg::tri::io::Importer<MyMesh>::Open(mesh,argv[1]);
  if(err) 
  {
    printf("Unable to open mesh %s : '%s'\n",argv[1],vcg::tri::io::Importer<MyMesh>::ErrorMsg(err));
    exit(-1);
  }
	printf("mesh loaded %d %d \n",mesh.vn,mesh.fn);
  

MyTriEdgeCollapse::SetDefaultParams();
  MyTriEdgeCollapse::Params().QualityThr  =.3;
  double TargetError=numeric_limits<double>::max();
  bool CleaningFlag =false;
     // parse command line.
	  for(int i=4; i < argc;)
    {
      if(argv[i][0]=='-')
        switch(argv[i][1])
      { 
        case 'H' : MyTriEdgeCollapse::Params().SafeHeapUpdate=true; printf("Using Safe heap option\n"); break;
        case 'Q' : if(argv[i][2]=='y') { MyTriEdgeCollapse::Params().QualityCheck	= true;  printf("Using Quality Checking\n");	}
                                  else { MyTriEdgeCollapse::Params().QualityCheck	= false; printf("NOT Using Quality Checking\n");	}                break;		
				case 'N' : if(argv[i][2]=='y') { MyTriEdgeCollapse::Params().NormalCheck	= true;  printf("Using Normal Deviation Checking\n");	}
                                  else { MyTriEdgeCollapse::Params().NormalCheck	= false; printf("NOT Using Normal Deviation Checking\n");	}        break;		
				case 'O' : if(argv[i][2]=='y') { MyTriEdgeCollapse::Params().OptimalPlacement	= true;  printf("Using OptimalPlacement\n");	}
                                  else { MyTriEdgeCollapse::Params().OptimalPlacement	= false; printf("NOT Using OptimalPlacement\n");	}        break;		
				case 'S' : if(argv[i][2]=='y') { MyTriEdgeCollapse::Params().ScaleIndependent	= true;  printf("Using ScaleIndependent\n");	}
                                  else { MyTriEdgeCollapse::Params().ScaleIndependent	= false; printf("NOT Using ScaleIndependent\n");	}        break;		
				case 'B' : if(argv[i][2]=='y') { MyTriEdgeCollapse::Params().PreserveBoundary	= true;  printf("Preserving Boundary\n");	}
                                  else { MyTriEdgeCollapse::Params().PreserveBoundary	= false; printf("NOT Preserving Boundary\n");	}        break;		
				case 'T' : if(argv[i][2]=='y') { MyTriEdgeCollapse::Params().PreserveTopology	= true;  printf("Preserving Topology\n");	}
                                  else { MyTriEdgeCollapse::Params().PreserveTopology	= false; printf("NOT Preserving Topology\n");	}        break;		
				case 'q' :	MyTriEdgeCollapse::Params().QualityThr	= atof(argv[i]+2);	           printf("Setting Quality Thr to %f\n",atof(argv[i]+2)); 	 break;			
				case 'n' :	MyTriEdgeCollapse::Params().NormalThr		= atof(argv[i]+2)*M_PI/180.0;  printf("Setting Normal Thr to %f deg\n",atof(argv[i]+2)); break;	
				case 'b' :	MyTriEdgeCollapse::Params().BoundaryWeight  = atof(argv[i]+2);			printf("Setting Boundary Weight to %f\n",atof(argv[i]+2)); break;		
				case 'e' :	TargetError = atof(argv[i]+2);			printf("Setting TargetError to %g\n",atof(argv[i]+2)); break;		
				case 'P' :	CleaningFlag=true;  printf("Cleaning mesh before simplification\n",atof(argv[i]+2)); break;	

				default  :  printf("Unknown option '%s'\n", argv[i]);
          exit(0);
      }
      i++;
    }



  if(CleaningFlag){
      int dup = tri::Clean<MyMesh>::RemoveDuplicateVertex(mesh);
      int unref =  tri::Clean<MyMesh>::RemoveUnreferencedVertex(mesh);
      printf("Removed %i duplicate and %i unreferenced vertices from mesh \n",dup,unref);
  }


  printf("reducing it to %i\n",FinalSize);
	
	vcg::tri::UpdateBounding<MyMesh>::Box(mesh);
  
	// decimator initialization
  vcg::LocalOptimization<MyMesh> DeciSession(mesh);
	
	int t1=clock();		
	DeciSession.Init<MyTriEdgeCollapse >();
  int t2=clock();	
  printf("Initial Heap Size %i\n",DeciSession.h.size());

	DeciSession.SetTargetSimplices(FinalSize);
	DeciSession.SetTimeBudget(0.5f);
  if(TargetError< numeric_limits<double>::max() ) DeciSession.SetTargetMetric(TargetError);

  while(DeciSession.DoOptimization() && mesh.fn>FinalSize && DeciSession.currMetric < TargetError)
    printf("Current Mesh size %7i heap sz %9i err %9g \r",mesh.fn,DeciSession.h.size(),DeciSession.currMetric);
  int t3=clock();	
  printf("mesh  %d %d Error %g \n",mesh.vn,mesh.fn,DeciSession.currMetric);
  printf("\nCompleted in (%i+%i) msec\n",t2-t1,t3-t2);
	
  vcg::tri::io::ExporterPLY<MyMesh>::Save(mesh,argv[2]);
	return 0;

}
