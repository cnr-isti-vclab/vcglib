// stuff to define the mesh
#include <vcg/complex/complex.h>

// io
#include <wrap/io_trimesh/import.h>
#include <wrap/io_trimesh/export_ply.h>

// local optimization
#include <vcg/complex/algorithms/local_optimization.h>
#include <vcg/complex/algorithms/local_optimization/tri_edge_collapse_quadric.h>

using namespace vcg;
using namespace tri;

/**********************************************************
Mesh Classes for Quadric Edge collapse based simplification

For edge collpases we need verteses with:
- V->F adjacency
- per vertex incremental mark
- per vertex Normal


Moreover for using a quadric based collapse the vertex class
must have also a Quadric member Q();
Otherwise the user have to provide an helper function object
to recover the quadric.

******************************************************/
// The class prototypes.
class MyVertex;
class MyEdge;
class MyFace;

struct MyUsedTypes: public UsedTypes<Use<MyVertex>::AsVertexType,Use<MyEdge>::AsEdgeType,Use<MyFace>::AsFaceType>{};

class MyVertex  : public Vertex< MyUsedTypes,
  vertex::VFAdj,
  vertex::Coord3f,
  vertex::Normal3f,
  vertex::Mark,
  vertex::BitFlags  >{
public:
  vcg::math::Quadric<double> &Qd() {return q;}
private:
  math::Quadric<double> q;
  };

class MyEdge : public Edge< MyUsedTypes> {};

typedef BasicVertexPair<MyVertex> VertexPair;

class MyFace    : public Face< MyUsedTypes,
  face::VFAdj,
  face::VertexRef,
  face::BitFlags > {};

// the main mesh class
class MyMesh    : public vcg::tri::TriMesh<std::vector<MyVertex>, std::vector<MyFace> > {};


class MyTriEdgeCollapse: public vcg::tri::TriEdgeCollapseQuadric< MyMesh, VertexPair, MyTriEdgeCollapse, QInfoStandard<MyVertex>  > {
            public:
            typedef  vcg::tri::TriEdgeCollapseQuadric< MyMesh,  VertexPair, MyTriEdgeCollapse, QInfoStandard<MyVertex>  > TECQ;
            typedef  MyMesh::VertexType::EdgeType EdgeType;
            inline MyTriEdgeCollapse(  const VertexPair &p, int i, BaseParameterClass *pp) :TECQ(p,i,pp){}
};

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
      "Copyright 2003-2012 Visual Computing Lab I.S.T.I. C.N.R.\n"
      "\nUsage:  "\
      "tridecimator fileIn fileOut face_num [opt]\n"\
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

// mesh to simplify
MyMesh mesh;

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

  TriEdgeCollapseQuadricParameter qparams;
  qparams.QualityThr  =.3;
  float TargetError=std::numeric_limits<float>::max();
  bool CleaningFlag =false;
     // parse command line.
    for(int i=4; i < argc;)
    {
      if(argv[i][0]=='-')
        switch(argv[i][1])
      {
        case 'H' : qparams.SafeHeapUpdate=true; printf("Using Safe heap option\n"); break;
        case 'Q' : if(argv[i][2]=='y') { qparams.QualityCheck	= true;  printf("Using Quality Checking\n");	}
                                  else { qparams.QualityCheck	= false; printf("NOT Using Quality Checking\n");	}                break;
        case 'N' : if(argv[i][2]=='y') { qparams.NormalCheck	= true;  printf("Using Normal Deviation Checking\n");	}
                                  else { qparams.NormalCheck	= false; printf("NOT Using Normal Deviation Checking\n");	}        break;
        case 'O' : if(argv[i][2]=='y') { qparams.OptimalPlacement	= true;  printf("Using OptimalPlacement\n");	}
                                  else { qparams.OptimalPlacement	= false; printf("NOT Using OptimalPlacement\n");	}        break;
        case 'S' : if(argv[i][2]=='y') { qparams.ScaleIndependent	= true;  printf("Using ScaleIndependent\n");	}
                                  else { qparams.ScaleIndependent	= false; printf("NOT Using ScaleIndependent\n");	}        break;
        case 'B' : if(argv[i][2]=='y') { qparams.PreserveBoundary	= true;  printf("Preserving Boundary\n");	}
                                  else { qparams.PreserveBoundary	= false; printf("NOT Preserving Boundary\n");	}        break;
        case 'T' : if(argv[i][2]=='y') { qparams.PreserveTopology	= true;  printf("Preserving Topology\n");	}
                                  else { qparams.PreserveTopology	= false; printf("NOT Preserving Topology\n");	}        break;
        case 'q' :	qparams.QualityThr	= atof(argv[i]+2);	           printf("Setting Quality Thr to %f\n",atof(argv[i]+2)); 	 break;
        case 'n' :	qparams.NormalThrRad = math::ToRad(atof(argv[i]+2));  printf("Setting Normal Thr to %f deg\n",atof(argv[i]+2)); break;
        case 'b' :	qparams.BoundaryWeight  = atof(argv[i]+2);			printf("Setting Boundary Weight to %f\n",atof(argv[i]+2)); break;
        case 'e' :	TargetError = float(atof(argv[i]+2));			printf("Setting TargetError to %g\n",atof(argv[i]+2)); break;
        case 'P' :	CleaningFlag=true;  printf("Cleaning mesh before simplification\n"); break;

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
  vcg::LocalOptimization<MyMesh> DeciSession(mesh,&qparams);

  int t1=clock();
  DeciSession.Init<MyTriEdgeCollapse>();
  int t2=clock();
  printf("Initial Heap Size %i\n",int(DeciSession.h.size()));

  DeciSession.SetTargetSimplices(FinalSize);
  DeciSession.SetTimeBudget(0.5f);
  if(TargetError< std::numeric_limits<float>::max() ) DeciSession.SetTargetMetric(TargetError);

  while(DeciSession.DoOptimization() && mesh.fn>FinalSize && DeciSession.currMetric < TargetError)
    printf("Current Mesh size %7i heap sz %9i err %9g \r",mesh.fn, int(DeciSession.h.size()),DeciSession.currMetric);

  int t3=clock();
  printf("mesh  %d %d Error %g \n",mesh.vn,mesh.fn,DeciSession.currMetric);
  printf("\nCompleted in (%i+%i) msec\n",t2-t1,t3-t2);

  vcg::tri::io::ExporterPLY<MyMesh>::Save(mesh,argv[2]);
    return 0;

}
