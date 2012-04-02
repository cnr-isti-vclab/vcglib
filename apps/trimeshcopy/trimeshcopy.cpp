#include <vcg/complex/append.h>

// stuff to define the mesh
#include <vcg/simplex/vertex/base.h>
#include <vcg/simplex/face/base.h>
#include <vcg/simplex/edge/base.h>
#include <vcg/complex/complex.h>
// io
#include <wrap/io_trimesh/import.h>
#include <wrap/io_trimesh/export_ply.h>

#include <cstdlib>

#include <sys/timeb.h>
#include <iostream>


class MyVertex;
class MyEdge;
class MyFace;

struct MyUsedTypes: public vcg::UsedTypes<vcg::Use<MyVertex>::AsVertexType,vcg::Use<MyEdge>::AsEdgeType,vcg::Use<MyFace>::AsFaceType>{};

class MyVertex  : public vcg::Vertex< MyUsedTypes,vcg::vertex::VFAdj,vcg::vertex::Coord3f,vcg::vertex::Normal3f,vcg::vertex::Mark,vcg::vertex::BitFlags  >
{
};

class MyEdge : public vcg::Edge< MyUsedTypes> {};

class MyFace    : public vcg::Face< MyUsedTypes,
	vcg::face::VFAdj,
	vcg::face::VertexRef,
	vcg::face::BitFlags > {};

// the main mesh class
class MyMesh    : public vcg::tri::TriMesh<std::vector<MyVertex>, std::vector<MyFace> > {};

void Usage()
{
	printf(
		"---------------------------------\n"
		"         TriMeshCopy V.1.0 \n"
		"     http://vcg.isti.cnr.it\n"
		"    http://vcg.sourceforge.net\n"
		"   release date: "__DATE__"\n"
		"---------------------------------\n\n"
		"TriMeshCopy 1.0 \n"__DATE__"\n"
		"Copyright 2003-2012 Visual Computing Lab I.S.T.I. C.N.R.\n"
		"\nUsage:  "\
		"trimeshcopy fileIn [fileOut]\n"\
		"trimeshcopy test vcg::MeshCopy efficiency.\nIt imports a fileIn file into a user defined mesh and test how long vcg::MeshCopy needs to copy the imported mesh in a second one.The copy time is expressed in milliseconds.\nA fileOut file can be passed to the tool in order to check if the mesh was successfully copied.\nThe file will be exported in PLY file format.\n"
		);
	exit(-1);
}

int main(int argc ,char**argv)
{
	MyMesh mesh;
	if(argc<2) 
		Usage();

	timeb start;
	timeb end;
	ftime(&start);
	int err=vcg::tri::io::Importer<MyMesh>::Open(mesh,argv[1]);
	if(err)
	{
		std::cerr << "Unable to open mesh " << argv[1] << " : " << vcg::tri::io::Importer<MyMesh>::ErrorMsg(err) << std::endl;
		exit(-1);
	}
	ftime(&end);
	int loadtime = (end.time * 1000 + end.millitm) - (start.time * 1000 + start.millitm);
	std::cout << "mesh loaded in " << loadtime << " msecs. Verts: " << mesh.vn << " Faces: " << mesh.fn << "\n";
	MyMesh mm;
	ftime(&start);
	vcg::tri::Append<MyMesh,MyMesh>::MeshCopy(mm,mesh);
	ftime(&end);
	int cptime = (end.time * 1000 + end.millitm) - (start.time * 1000 + start.millitm);
	std::cout << "mesh copied in " << cptime << " msecs." << std::endl;

	if (argc == 3)
		vcg::tri::io::ExporterPLY<MyMesh>::Save(mm,argv[2]);
	return 0;
}
