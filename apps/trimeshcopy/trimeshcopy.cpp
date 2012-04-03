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
#include <string>


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

class OcfVertex;
class OcfEdge;
class OcfFace;

// Declaration of the semantic of the used types
class OcfUsedTypes: public vcg::UsedTypes < vcg::Use<OcfVertex>::AsVertexType,
	vcg::Use<OcfEdge   >::AsEdgeType,
	vcg::Use<OcfFace  >::AsFaceType >{};


// The Main Vertex Class
// Most of the attributes are optional and must be enabled before use.
// Each vertex needs 40 byte, on 32bit arch. and 44 byte on 64bit arch.

class OcfVertex  : public vcg::Vertex< OcfUsedTypes,vcg::vertex::InfoOcf,vcg::vertex::Coord3f,vcg::vertex::BitFlags,vcg::vertex::Normal3fOcf,vcg::vertex::VFAdjOcf,vcg::vertex::MarkOcf>
{
};


// The Main Edge Class
// Currently it does not contains anything.
class OcfEdge : public vcg::Edge<OcfUsedTypes>
{
};

// Each face needs 32 byte, on 32bit arch. and 48 byte on 64bit arch.
class OcfFace    : public vcg::Face<  OcfUsedTypes,vcg::face::InfoOcf,vcg::face::VertexRef,vcg::face::BitFlags,vcg::face::VFAdjOcf> {};

class OcfMesh    : public vcg::tri::TriMesh< vcg::vertex::vector_ocf<OcfVertex>, vcg::face::vector_ocf<OcfFace> > 
{
};

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
		"trimeshcopy fileIn -(n|o) [fileOut]\n"\
		"trimeshcopy test vcg::MeshCopy efficiency.\nIt imports a fileIn file into a user defined mesh and test how long vcg::MeshCopy needs to copy the imported mesh in a second one.The copy time is expressed in milliseconds.\nIf the -n flag is used a non-optional attributes mesh will be tested, defining -o, instead, the target mesh will be an ocf one.\nA fileOut file can be passed to the tool in order to check if the mesh was successfully copied.\nThe file will be exported in PLY file format.\n"
		);
	exit(-1);
}

int main(int argc ,char**argv)
{
	MyMesh mesh;
	if(argc<3) 
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

	std::string tmp(argv[2]);
	if (tmp == "-n")
	{
		MyMesh mm;
		ftime(&start);
		vcg::tri::Append<MyMesh,MyMesh>::MeshCopy(mm,mesh);
		ftime(&end);
		int cptime = (end.time * 1000 + end.millitm) - (start.time * 1000 + start.millitm);
		std::cout << "mesh copied in " << cptime << " msecs." << std::endl;

		if (argc == 4)
			vcg::tri::io::ExporterPLY<MyMesh>::Save(mm,argv[3]);
		return 0;
	}

	//if (tmp == "-o")
	//{
	//	OcfMesh ocfm;
	//	ftime(&start);
	//	vcg::tri::Append<OcfMesh,MyMesh>::MeshCopy(ocfm,mesh);
	//	ftime(&end);
	//	cptime = (end.time * 1000 + end.millitm) - (start.time * 1000 + start.millitm);
	//	std::cout << "mesh copied in " << cptime << " msecs." << std::endl;

	//	if (argc == 4)
	//		vcg::tri::io::ExporterPLY<OcfMesh>::Save(ocfm,argv[3]);

	//	return 0;
	//}

	Usage();
	return 0;
}
