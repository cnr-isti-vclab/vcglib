/****************************************************************************
* VCGLib                                                            o o     *
* Visual and Computer Graphics Library                            o     o   *
*                                                                _   O  _   *
* Copyright(C) 2005                                                \/)\/    *
* Visual Computing Lab                                            /\/|      *
* ISTI - Italian National Research Council                           |      *
*                                                                    \      *
* All rights reserved.                                                      *
*                                                                           *
* This program is free software; you can redistribute it and/or modify      *   
* it under the terms of the GNU General Public License as published by      *
* the Free Software Foundation; either version 2 of the License, or         *
* (at your option) any later version.                                       *
*                                                                           *
* This program is distributed in the hope that it will be useful,           *
* but WITHOUT ANY WARRANTY; without even the implied warranty of            *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
* GNU General Public License (http://www.gnu.org/licenses/gpl.txt)          *
* for more details.                                                         *
*                                                                           *
****************************************************************************/
/****************************************************************************
History

$Log: not supported by cvs2svn $
Revision 1.13  2005/11/17 00:42:03  cignoni
Changed include order
removed clean::initialize

Revision 1.12  2005/11/16 16:45:51  rita_borgo
Minor changes

Revision 1.11  2005/11/16 15:59:46  cignoni
Changed name of the component from <Flag> to <BitFlags>

Revision 1.10  2005/11/14 09:21:07  cignoni
Heavily restructured the code of Trimeshinfo.
Now divided the collecting part from the reporting one (xml and ascii)

Revision 1.9  2005/11/04 15:37:57  rita_borgo
Removed Debug option

Revision 1.8  2005/10/11 16:03:40  rita_borgo
Moved all the main functions inside clean.h

Revision 1.7	2005/09/30 15:48:46	 rita_borgo
Fixed	manifold Test

Revision 1.6 2005/09/30	14:13:01 rita_borgo
Problem: Text	not	aligned

Revision 1.5 2005/09/30	13:29:40 rita_borgo
Fixed	Manifold Test

Revision 1.4 2005/09/29	14:48:15 rita_borgo
Fixed	code related to	creation of	the	XML	file

Revision 1.3 2005/09/28	13:57:09 rita_borgo
Fixed	some printout	not	alligned

Revision 1.2 2005/09/28	10:46:04 rita_borgo
Added	possibility	of saving	File in	OFF	format

Revision 1.1 2005/09/20	10:15:27 rita_borgo
Changed	file name	to uniform with	other	solution projects,
before was main.cpp

Revision 1.8 2005/02/15	12:26:06 rita_borgo
Minor	changes	to self-intersection

Revision 1.7 2005/02/07	15:44:31 rita_borgo
Fixed	Color	and	Volume

Revision 1.6 2005/02/01	17:37:53 rita_borgo
Fixed	Volume and Color

Revision 1.5 2005/01/18	16:33:12 rita_borgo
Added	OFF	file Option

Revision 1.4 2005/01/17	18:19:00 rita_borgo
Added	new	routines.
Self-intersection	first	release

Revision 1.2 2005/01/03	16:13:09 rita_borgo
Added	Standard comments

****************************************************************************/
#include <vector>
#include <string>
#include <stack>

using namespace std;

#include <vcg/complex/trimesh/base.h>
#include <vcg/complex/trimesh/update/topology.h>
#include <vcg/complex/trimesh/update/edges.h>
#include <vcg/complex/trimesh/update/bounding.h>
#include <vcg/complex/trimesh/update/flag.h>
#include <vcg/complex/trimesh/clean.h>
#include <vcg/space/intersection/triangle_triangle3.h>

#include <vcg/math/histogram.h>
#include <wrap/io_trimesh/import.h>
#include <wrap/io_trimesh/export.h>

#include <vcg/simplex/vertexplus/base.h>
#include <vcg/simplex/vertexplus/component.h>

#include <vcg/simplex/faceplus/base.h>
#include <vcg/simplex/faceplus/component.h>
#include <vcg/simplex/face/pos.h> 

#include "XMLTree.h"

#include <vcg/space/index/grid_static_ptr.h>
#include "defs.h"

using namespace std;
using namespace vcg;


class CFace;
class CEdge;
class CVertex  : public VertexSimp2< CVertex, CEdge, CFace, vert::Coord3f, vert::BitFlags, vert::Normal3f >{};
class CFace    : public FaceSimp2< CVertex, CEdge, CFace, face::FFAdj, face::VertexRef, face::BitFlags > {};
class CMesh    : public vcg::tri::TriMesh< vector<CVertex>, vector<CFace> > {};

typedef CMesh::VertexPointer VertexPointer;
typedef CMesh::VertexIterator VertexIterator;
typedef Point3<CMesh::ScalarType> Point3x;
typedef vector<Point3x> Hole;

void OpenMesh(const	char *filename,	CMesh &m)
{
	int	err	= tri::io::Importer<CMesh>::Open(m,filename);
	if(err)
	{
		printf("Error in reading %s: '%s'\n",filename, 
			tri::io::Importer<CMesh>::ErrorMsg(err));
		exit(-1);
	}
	printf("\t read mesh `%s'\n", filename);
}

typedef CMesh::VertexPointer VertexPointer;
typedef CMesh::VertexIterator VertexIterator;
typedef CMesh::FaceContainer FaceContainer;
typedef CMesh::ScalarType ScalarType;

struct MeshInfo
{
	string FileName;
	int vn,fn;
	bool Manifold;
	int count_e,boundary_e,count_fd,count_uv,numholes;
	int BEdges;
	float Volume;
	int numcomponents,Genus;
	bool Regular,Semiregular;
	bool Orientable,Oriented;
	int dv;
	bool SelfIntersect;
};

void initMeshInfo(MeshInfo &mi)
{
	mi.vn = mi.fn = mi.count_e = mi.boundary_e = mi.count_fd = mi.count_uv = 
		mi.numholes = mi.BEdges = mi.numcomponents = mi.Genus = mi.dv = 0;
	mi.Volume = 0;
}
void PrintAsciiInfo(MeshInfo &mi)
{
	printf("\t Mesh info:\n");
	printf(" \t M: '%s'\n\t Number	of vertices: %d	\n", mi.FileName.c_str(), mi.vn);
	printf("\t Number of faces: %d \n",	mi.fn);

	if (!mi.Manifold) printf("\t Manifold: NO\n");
	else
		printf("\t Manifold: YES\n");

	printf("\t Number of edges: %d \n", mi.count_e);
	printf("\t Number of internal edges: %d \n", mi.count_e-mi.boundary_e);
	printf("\t Number of boundary edges: %i \n", mi.boundary_e);
	printf("\t Number of degenerated faces: %d\n",	mi.count_fd);
	printf("\t Number of unreferenced	vertices: %d\n",mi.count_uv);
	printf("\t Number of holes/boundaries: %d \n", mi.numholes);
	if(mi.Volume) printf("\t Volume: %f \n", mi.Volume);
	else 
		printf("\t Volume: UNDEFINED, mesh is either non-manifold or has holes \n");

	printf("\t Number of connected components: %d\n",	mi.numcomponents);
	if (mi.Manifold)
		printf("\t Genus: %d \n", mi.Genus);
	else
		printf("\t Genus (n/a)\n");

	if (mi.Regular) 
		printf("\t Type of Mesh: REGULAR\n");
	else if (mi.Semiregular)
		printf("\t Type of Mesh: SEMIREGULAR\n");
	else
		printf("\t Type of Mesh: IRREGULAR\n");

	if (!mi.Manifold)
	{
		printf("\t Orientable Mesh: NO\n");
		printf("\t Oriented Mesh: NO\n");
	}
	else
	{
		printf("\t Orientable Mesh: %s\n",mi.Orientable?"Yes":"No");
		printf("\t Oriented Mesh: %s\n",mi.Oriented?"Yes":"No");
	}

	printf("\t Number of duplicated vertices found: %d\n",mi.dv);
	printf("\t Self Intersection: %s\n",mi.SelfIntersect?"Yes":"No");
}

void PrintXMLInfo(MeshInfo &mi)
{
	XMLTree	doc;
	doc.initializeMain();

	char s[256];
	sprintf(s,"%d",mi.vn);
	doc.addNode(s, VALUE_INTEGER, "Number of	Vertices");
	sprintf(s,"%d",mi.fn);
	doc.addNode(s, VALUE_INTEGER,	"Number	of Faces");

	if(mi.Manifold)
		doc.addNode("No", VALUE_BOOL,"Manifold");
	else
		doc.addNode("Yes", VALUE_BOOL,"Manifold");

	sprintf(s,"%d",mi.count_e);
	doc.addNode(s, VALUE_INTEGER,"Number of Edges");
	sprintf(s,"%d",mi.count_fd);
	doc.addNode(s, VALUE_INTEGER,"Number of Degenerated Faces");

	sprintf(s,"%d",mi.count_uv);
	doc.addNode(s, VALUE_INTEGER,"Number of unreferenced vertices");
	sprintf(s,"%d",mi.numholes);
	doc.addNode(s, VALUE_INTEGER,"Number of Holes");
	sprintf(s,"%d",mi.BEdges);
	doc.addNode(s, VALUE_INTEGER,"Number of Border Edges");
	sprintf(s,"%f",mi.Volume);
	doc.addNode(s, VALUE_FLOAT,"Volume");
	sprintf(s,"%d",mi.numcomponents);
	doc.addNode(s, VALUE_INTEGER,"Number of Connected Components");
	sprintf(s,"%d",mi.Genus);		doc.addNode(s, VALUE_INTEGER,"Genus");

	if (mi.Regular)
		doc.addNode("REGULAR", VALUE_STRING,"Type of Mesh");
	else if (mi.Semiregular)
		doc.addNode("SEMIREGULAR", VALUE_STRING,"Type of Mesh");
	else 
		doc.addNode("IRREGULAR",   VALUE_STRING,"Type of Mesh");
	
	if (!mi.Manifold) 
	{
		doc.addNode("NO", VALUE_STRING,"Orientable Mesh");
		doc.addNode("NO", VALUE_STRING,"Oriented Mesh");
	}
	else
	{
		doc.addNode(mi.Orientable?"Yes":"No", VALUE_STRING,"Orientable Mesh");
		doc.addNode(  mi.Oriented?"Yes":"No", VALUE_STRING,"Oriented Mesh");
	}

	sprintf(s,"%d",mi.dv);		doc.addNode(s, VALUE_INTEGER,"Duplicated Vertices");
	doc.addNode(  mi.SelfIntersect?"Yes":"No", VALUE_STRING,"Self	Intersection");

	doc.finalizeMain();
	doc.printXMLTree();
}

int main(int	argc,char	** argv)
{
	CMesh m;
	bool SaveFlag=false;
	bool AsciiFlag=true;
	bool XmlFlag=false;

	string SaveName;
	
	MeshInfo mi;
	initMeshInfo(mi);
	printf("-------------------------------\n"
		"    TriMeshInfo V.1.2 \n"
		"    http://vcg.isti.cnr.it\n"
		"    release date: "__DATE__"\n"
		"-------------------------------\n\n");


	// load input meshes.
	if (argc <= 1)
	{
		printf(MSG_ERR_N_ARGS);
		exit(-1);
	}
	mi.FileName=argv[1];

	for(int i=3; i < argc; i++)
	{
		if(argv[i][0]=='-')
		{
			switch(argv[i][1])
			{
			case 'o' :
				SaveFlag=true; SaveName=argv[i][2];
			break;
			case 'x' :
				if(argv[i][2]=='y')
				{
					XmlFlag = true;
					printf("Enable XML Printing\n");
				}
				else 
				{
					XmlFlag = false;
					printf("Disable XML Printing\n");
				}
			break;
			case 'a' :
				if(argv[i][2]=='y')
				{
					AsciiFlag = true;
					printf("Enable Ascii Printing\n");
				}
				else
				{
					AsciiFlag = false;
					printf("Disable Ascii Printing\n");
				}
			break;
			default:
				printf(MSG_ERR_INVALID_OPTION, argv[i]);
				exit(0);
			}
		i++;
	}

	OpenMesh(mi.FileName.c_str(),m);
	mi.vn=m.vn;
	mi.fn=m.fn;

	// DEGENERATED FACES
	mi.count_fd = tri::Clean<CMesh>::DegeneratedFaces(m);
	
	vcg::tri::UpdateTopology<CMesh>::FaceFace(m);

	// UNREFERENCED	VERTEX
	mi.count_uv = tri::Clean<CMesh>::DetectUnreferencedVertex(m);

	tri::UpdateFlags<CMesh>::Clear(m);

	// IS	MANIFOLD
	mi.Manifold	=	tri::Clean<CMesh>::IsComplexManifold(m);	

	// COUNT EDGES
	tri::Clean<CMesh>::CountEdges(m,	mi.count_e, mi.boundary_e);

	// HOLES COUNT
	if(mi.Manifold)
	{
		mi.numholes = tri::Clean<CMesh>::CountHoles(m);
		mi.BEdges = tri::Clean<CMesh>::BorderEdges(m, mi.numholes);
	}

	// Mesh	Volume
	if(mi.numholes==0) mi.Volume	=	m.Volume();
	
	// CONNECTED COMPONENTS
	mi.numcomponents = tri::Clean<CMesh>::ConnectedComponents(m);

	if(mi.Manifold) 
		mi.Genus = tri::Clean<CMesh>::MeshGenus(m,mi.count_uv, mi.numholes, mi.numcomponents, mi.count_e);
	
	// REGULARITY
	tri::Clean<CMesh>::IsRegularMesh(m, mi.Regular,	mi.Semiregular);
	
	// ORIENTABLE E ORIENTED MESH
	if (mi.Manifold) tri::Clean<CMesh>::IsOrientedMesh(m,	mi.Oriented,	mi.Orientable);

	mi.dv = tri::Clean<CMesh>::RemoveDuplicateVertex(m);

	// SELF INTERSECTION
	mi.SelfIntersect = tri::Clean<CMesh>::SelfIntersections(m);

	if(SaveFlag) tri::io::Exporter<CMesh>::Save(m,SaveName.c_str());

	if(AsciiFlag) PrintAsciiInfo(mi);
	if(XmlFlag) PrintXMLInfo(mi);
	return 0;
}

