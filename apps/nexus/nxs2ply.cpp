// mesh definition 
#include <vcg/simplex/vertex/with/afvn.h>
#include <vcg/simplex/face/with/af.h>
#include <vcg/complex/trimesh/base.h>

#include <vcg/complex/trimesh/update/normal.h>

// input output
#include <wrap/io_trimesh/import_ply.h>
#include <wrap/io_trimesh/export_ply.h>

#include "nxsexport.h"
#include "extraction.h"
#include "metric.h"

using namespace vcg;
using namespace nxs;

struct MyFace;
struct MyTetra;
struct MyEdge;
struct MyVertex: public VertexAFVNf<MyEdge,MyFace,MyTetra>{};
struct MyFace:   public FaceAF<MyVertex,MyEdge,MyFace>{};
struct MyMesh:   public tri::TriMesh< vector<MyVertex>, vector<MyFace> >{};

int main(int argc, char *argv[]) {

  string input;
  string output;
  float target_error = 0;

  int option;
  while((option = getopt(argc, argv, "e:")) != EOF) {
    switch(option) {
    case 'e': target_error = atof(optarg); break;
    case 'o': output = optarg; break;
    }
  }
  
  if(optind != argc -1) {
    cerr << "Usage: " << argv[0] << " <nexus> [options]\n"
	 << " -e <float>   : extraction target error\n"
	 << " -o <file.ply>: output name\n\n";
    return -1;
  }
  input = argv[optind];
  if(!output.size()) output = input + ".ply";

  Nexus nexus;
  if(!nexus.Load(input.c_str())) {
    cerr << "Could not open file: " << input.c_str() << ".nxs\n";
    return -1;
  }

  Extraction extr;
  extr.SetMetric(new FlatMetric);
  extr.target_error = target_error;
  extr.Extract(&nexus);

  vector<unsigned int> selected;
  for(unsigned int i = 0; i < extr.selected.size(); i++)
    selected.push_back(extr.selected[i].id);


  MyMesh m;
  ExportTriMesh<MyMesh>(nexus, selected, m);

  //write to ply:
  vcg::tri::io::PlyInfo pi;
  vcg::tri::io::ExporterPLY<MyMesh>::Save(m, output.c_str(), pi.mask);
  return 0;
}
