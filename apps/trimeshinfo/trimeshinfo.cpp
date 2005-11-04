/****************************************************************************
*	VCGLib o o *
*	Visual and Computer	Graphics Library o o *
*	_	O	_	*
*	Copyright(C) 2004	\/)\/	*
*	Visual Computing Lab /\/|	*
*	ISTI - Italian National	Research Council | *
*	\	*
*	All	rights reserved. *
*	*
*	This program is	free software; you can redistribute	it and/or	modify *
*	it under the terms of	the	GNU	General	Public License as	published	by *
*	the	Free Software	Foundation;	either version 2 of	the	License, or	*
*	(at	your option) any later version.	*
*	*
*	This program is	distributed	in the hope	that it	will be	useful,	*
*	but	WITHOUT	ANY	WARRANTY;	without	even the implied warranty	of *
*	MERCHANTABILITY	or FITNESS FOR A PARTICULAR	PURPOSE. See the *
*	GNU	General	Public License (http://www.gnu.org/licenses/gpl.txt) *
*	for	more details.	*
*	*
****************************************************************************/
/****************************************************************************
History

$Log: not supported by cvs2svn $
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
using	namespace	std;


#include<vcg/simplex/face/with/afav.h>
#include<vcg/simplex/face/pos.h> 
#include<vcg/complex/trimesh/base.h>

#include<vcg/complex/trimesh/update/topology.h>
#include <vcg/complex/trimesh/update/edges.h>
#include <vcg/complex/trimesh/update/bounding.h>
#include <vcg/complex/trimesh/clean.h>
#include <vcg/space/intersection/triangle_triangle3.h>

#include <vcg/math/histogram.h>
#include <wrap/io_trimesh/import.h>
#include <wrap/io_trimesh/export_ply.h>
#include <wrap/io_trimesh/export_stl.h>
#include <wrap/io_trimesh/export_dxf.h>
#include "XMLTree.h"
#include <wrap/io_trimesh/export_off.h>
// loader
#include<wrap/io_trimesh/import_ply.h>

#include <vcg/space/index/grid_static_ptr.h>
#include "defs.h"
#include "trimeshtype.h"
string ans;

using	namespace	vcg;
using	namespace	tri;
using	namespace	face;






typedef	CMesh::VertexPointer	VertexPointer;
typedef	CMesh::VertexIterator VertexIterator;
typedef	Point3<CMesh::ScalarType> Point3x;
typedef	vector<Point3x>	Hole;

void OpenMesh(const	char *filename,	CMesh &m)
{
	int	err	=	tri::io::Importer<CMesh>::Open(m,filename);
	if(err){
		printf("Error	in reading %s: '%s'\n",filename,tri::io::Importer<CMesh>::ErrorMsg(err));
		exit(-1);
	}
	printf("\t read mesh `%s'\n", filename);
}

typedef	CMesh::VertexPointer	VertexPointer;
typedef	CMesh::VertexIterator VertexIterator;


typedef CMesh::FaceContainer  FaceContainer;

typedef CMesh::ScalarType	ScalarType;



void main(int	argc,char	** argv)
{
	CMesh m;
	XMLTree	doc;
	


	//load the mesh
	//argv[1]=(char*)"c:\\checkup\\debug\\column1m.ply";
	//argv[1]	=	"C:\\sf\\apps\\msvc\\trimeshinfo\\Release\\prism.off";
	//argv[1]	=	"C:\\sf\\apps\\msvc\\trimeshinfo\\Release\\prova1.ply";

	// print program info
	printf("-------------------------------\n"
		"	TriMeshInfo	V.1.01 \n"
		"	http://vcg.isti.cnr.it\n"
		"	release	date:	"__DATE__"\n"
		"-------------------------------\n\n");



	
	
		// load	input	meshes.
		if(argc	<= 1)
		{
			printf(MSG_ERR_N_ARGS);
			exit(-1);
		}
	



	OpenMesh(argv[1],m);

	doc.initializeMain();	

	printf("\t Mesh info:\n");
	printf(" \t M: '%s'\n\t Number	of vertices: %d	\n", argv[1],	m.vn);
	printf("\t Number of faces: %d \n",	m.fn);



	/*------------XML	file part	------------------*/
	char*	s	=new(char[25]);
	sprintf(s,"%d",m.vn);
	doc.addNode(s, VALUE_INTEGER,"Number of	Vertices");
	sprintf(s,"%d",m.fn);
	doc.addNode(s, VALUE_INTEGER,	"Number	of Faces");

	/*--------------------------------------------*/



	if(m.HasPerFaceColor()||m.HasPerVertexColor())
	{
		Color4b	Color=m.C();
		printf(	"\t	Object color(4b):	%f %f	%f \n",	Color[0],	Color[1],	Color[2]);

		/*------------XML	file part	------------------*/

		sprintf(s,"%f	%f %f	",Color[0],	Color[1],	Color[2]);
		doc.addNode(s, VALUE_FLOAT,"Colors");
		/*--------------------------------------------*/
	}


	vcg::tri::UpdateTopology<CMesh>::FaceFace(m);

	// IS	MANIFOLD


	tri::Clean<CMesh>::Initialize(m);



	bool Manifold	=	tri::Clean<CMesh>::IsComplexManifold(m);	

	if (!Manifold)
	{
		printf(	"\t Manifold: NO\n");

		/*------------XML	file part	------------------*/

		s	=	"No";
		doc.addNode(s, VALUE_BOOL,"Manifold");

		/*--------------------------------------------*/
	}
	else
	{

		printf(	"\t Manifold: YES\n");

		/*------------XML	file part	------------------*/
		s	=	"Yes";
		doc.addNode(s, VALUE_BOOL,"Manifold");
		/*--------------------------------------------*/
	}


	// COUNT EDGES
	int	count_e=0;
	int	boundary_e = 0;
	tri::Clean<CMesh>::CountEdges(m,	count_e, boundary_e);

	printf("\t Number of edges: %d \n", count_e);
	printf("\t Number of internal edges: %d \n", count_e-boundary_e);
	printf("\t Number of boundary edges: %i \n", boundary_e);

	/*------------XML	file part	------------------*/
	s	=	new(char[25]);
	sprintf(s,"%d",count_e);
	doc.addNode(s, VALUE_INTEGER,"Number of	Edges");

	/*--------------------------------------------*/



	// DEGENERATED FACES

	int	count_fd = tri::Clean<CMesh>::DegeneratedFaces(m);

	printf("\t Number of degenerated faces: %d\n",	count_fd);

	/*------------XML	file part	------------------*/

	sprintf(s,"%d",count_fd);
	doc.addNode(s, VALUE_INTEGER,"Number of	Degenerated	Faces");

	/*--------------------------------------------*/

	// UNREFERENCED	VERTEX

	int	count_uv = tri::Clean<CMesh>::DetectUnreferencedVertex(m);
	printf("\t Number of unreferenced	vertices: %d\n",count_uv);

	/*------------XML	file part	------------------*/

	sprintf(s,"%d",count_uv);
	doc.addNode(s, VALUE_INTEGER,"Number of	unreferenced vertices");

	/*--------------------------------------------*/


	// HOLES COUNT

	int	numholes =0	;
	if(Manifold)
	{
		numholes = tri::Clean<CMesh>::CountHoles(m);
		printf("\t Number of holes/boundaries: %d \n", numholes);

		/*------------XML	file part	------------------*/
		sprintf(s,"%d",numholes);
		doc.addNode(s, VALUE_INTEGER,"Number of	Holes");
		/*--------------------------------------------*/
		if(numholes>0)
		{
			printf("\t Edges per hole/boundary:\n\t	(");
			int	BEdges = tri::Clean<CMesh>::BorderEdges(m, numholes);

			/*------------XML	file part	------------------*/

			sprintf(s,"%d",BEdges);
			doc.addNode(s, VALUE_INTEGER,"Number of	Border Edges");

			/*--------------------------------------------*/
		}
	}
	else
		printf(	"\t Number of holes: UNDEFINED, mesh is non-manifold \n");


	// Mesh	Volume
	float	vol	=	m.Volume();
	int	nuh	=	numholes;
	if((m.Volume()!=0.)&&(Manifold)&&(numholes==0))
	{

		printf("\t Volume: %f \n", m.Volume());

		/*------------XML	file part	------------------*/

		sprintf(s,"%f",m.Volume());
		doc.addNode(s, VALUE_FLOAT,"Volume");
		/*--------------------------------------------*/
	}
	else
	{
		printf("\t Volume: UNDEFINED, mesh is either non-manifold or has holes \n");
	}


	// CONNECTED COMPONENTS


	int	numcomponents	=	tri::Clean<CMesh>::ConnectedComponents(m);
	printf("\t Number of connected components: %d\n",	numcomponents);

	/*------------XML	file part	------------------*/

	sprintf(s,"%d",numcomponents);
	doc.addNode(s, VALUE_INTEGER,"Number of	Connected	Components");
	/*--------------------------------------------*/

	//GENUS	-->	2( #components - genus ) = #vertices + #faces	-	#edge	-	#boundary_loops	=	eulernumber	-	#holes
	//eulero = (mesh.vn-count_uv)	-	(count_e)+mesh.fn;



	if(Manifold)
	{
		printf(	"\t Genus: %d \n", tri::Clean<CMesh>::MeshGenus(m,count_uv, numholes, numcomponents, count_e));

		/*------------XML	file part	------------------*/

		sprintf(s,"%d",tri::Clean<CMesh>::MeshGenus(m,count_uv, numholes, numcomponents, count_e));
		doc.addNode(s, VALUE_INTEGER,"Genus");
		/*--------------------------------------------*/
	}
	else//(!Manifold)
		printf(	"\t Genus: UNDEFINED, mesh is non-manifold \n");
	
	// REGULARITY
	bool Regular=true;
	bool Semiregular=true;

	tri::Clean<CMesh>::IsRegularMesh(m, Regular,	Semiregular);
	if (Regular)
	{

		printf("\t Type of Mesh: REGULAR\n");

		/*------------XML	file part	------------------*/
		s="REGULAR";
		doc.addNode(s, VALUE_STRING,"Type	of Mesh");
		/*--------------------------------------------*/
	}
	else if	(Semiregular)
	{

		printf("\t Type of Mesh: SEMIREGULAR\n");
		s="SEMIREGULAR";
		doc.addNode(s, VALUE_STRING,"Type	of Mesh");
	}
	else
	{

		printf("\t Type of Mesh: IRREGULAR\n");

		/*------------XML	file part	------------------*/

		s="IRREGULAR";
		doc.addNode(s, VALUE_STRING,"Type	of Mesh");

		/*--------------------------------------------*/
	}
	// ORIENTABLE	E	ORIENTED MESH

	bool Orientable=true;
	bool Oriented=true;
	if (!Manifold)
	{

		printf(	"\t Orientable Mesh: NO\n");

		/*------------XML	file part	------------------*/

		s="NO";
		doc.addNode(s, VALUE_STRING,"Orientable	Mesh");
		/*--------------------------------------------*/
	}
	else
	{
		tri::Clean<CMesh>::IsOrientedMesh(m,	Oriented,	Orientable);

		if (Orientable)
		{

			printf(	"\t Orientable Mesh: YES\n");

		/*------------XML	file part	------------------*/

			s="YES";
			doc.addNode(s, VALUE_STRING,"Orientable	Mesh");

		/*--------------------------------------------*/
		}
		else
		{

			printf(	"\t	Orientable Mesh: NO\n");

		/*------------XML	file part	------------------*/

			s="NO";
			doc.addNode(s, VALUE_STRING,"Orientable	Mesh");
		/*--------------------------------------------*/
		}
	}
	if (Oriented &&	Manifold)
	{

		printf(	"\t Oriented Mesh: YES\n");

		/*------------XML	file part	------------------*/

		s="YES";
		doc.addNode(s, VALUE_STRING,"Oriented	Mesh");
		/*--------------------------------------------*/
	}
	else
	{

		printf(	"\t Oriented Mesh: NO\n");

		/*------------XML	file part	------------------*/

		s="NO";
		doc.addNode(s, VALUE_STRING,"Oriented	Mesh");
		/*--------------------------------------------*/
	}
	int	dv = tri::Clean<CMesh>::RemoveDuplicateVertex(m);
	if(dv>0)
	{

		printf(	"\t Number of duplicated vertices found: %d\n",dv);

		/*------------XML	file part	------------------*/

		char*	s	=new(char[25]);
		sprintf(s,"%d",dv);
		doc.addNode(s, VALUE_INTEGER,"Duplicated Vertices");
		/*--------------------------------------------*/
	}
	else
	{

		printf(	"\t Duplicated vertices: NO\n");

		/*------------XML	file part	------------------*/

		s="NO";
		doc.addNode(s, VALUE_STRING,"Duplicated	Vertices");
		/*--------------------------------------------*/
	}
	// SELF	INTERSECTION


	printf("\t Model Bounding Box Diagonal: %f\n", m.bbox.Diag());
	if (tri::Clean<CMesh>::SelfIntersections(m))
	{

		printf(	"\t Self Intersection: YES\n");

		/*------------XML	file part	------------------*/

		s="YES";
		doc.addNode(s, VALUE_STRING,"Self	Intersection");

		/*--------------------------------------------*/
	}
	else
	{

		printf(	"\t Self Intersection: NO\n");

		/*------------XML	file part	------------------*/

		s="NO";
		doc.addNode(s, VALUE_STRING,"Self	Intersection");
		/*--------------------------------------------*/
	}


	
	


	string fs;

	cout<< "\t To save the file: [s/S]\n \t ";
	cin>>ans;


	if((ans	== "S")||(ans	== "s"))
	{
		cout<< "\t available formats: [ply, off, dxf, stl]"<<endl;
		cout<< "\t enter format \n \t ";
		cin>>ans;
		cout<<"\t enter filename \n \t ";
		cin>>fs;
		const	char*	filesave = fs.c_str();
		if(ans ==	"ply")
			tri::io::ExporterPLY<CMesh>::Save(m,	filesave);
		else if(ans	== "stl")
			tri::io::ExporterSTL<CMesh>::Save(m,filesave);
		else if(ans	== "dxf")
			tri::io::ExporterDXF<CMesh>::Save(m,filesave);

		else if(ans	== "off")
			tri::io::ExporterOFF<CMesh>::Save(m,filesave);
	}

	cout<<"\t create XML files? [y/Y|n/N] \n \t ";
	cin>>ans;
	if((ans=="Y")||(ans=="y"))
	{
		doc.finalizeMain();
		doc.printXMLTree();
	}
}

