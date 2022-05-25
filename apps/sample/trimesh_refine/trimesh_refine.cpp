/****************************************************************************
* VCGLib                                                            o o     *
* Visual and Computer Graphics Library                            o     o   *
*                                                                _   O  _   *
* Copyright(C) 2004-2016                                           \/)\/    *
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
#include<vcg/complex/complex.h>

#include <vcg/complex/algorithms/refine.h>
#include <vcg/complex/algorithms/refine_loop.h>
#include<vcg/complex/algorithms/refine_doosabin.h>


#include <vcg/complex/algorithms/bitquad_creation.h>

// input output
#include <wrap/io_trimesh/import.h>
#include <wrap/io_trimesh/export.h>

using namespace vcg;
using namespace std;


class MyEdge;    // dummy prototype never used
class MyFace;
class MyVertex;

struct MyUsedTypes : public UsedTypes<	Use<MyVertex>::AsVertexType, Use<MyFace>::AsFaceType>{};

class MyVertex  : public Vertex< MyUsedTypes, vertex::Coord3f, vertex::Normal3f, vertex::BitFlags  >{};
class MyFace    : public Face  < MyUsedTypes, face::InfoOcf, face::Normal3f, face::FFAdjOcf,  face::VertexRef, face::BitFlags > {};
class MyMesh    : public vcg::tri::TriMesh< vector<MyVertex>, face::vector_ocf<MyFace> > {};



class PFace;
class PVertex;

struct PUsedTypes: public UsedTypes<Use<PVertex>  ::AsVertexType,
									 Use<PFace>	::AsFaceType>{};

class PVertex:public Vertex<	PUsedTypes,
							  vertex::Coord3f,
							  vertex::Normal3f,
							  vertex::Mark,
							  vertex::Qualityf,
							  vertex::VFAdj,
							  vertex::BitFlags>{} ;

class PFace:public Face<
				  PUsedTypes
				  ,face::PolyInfo // this is necessary  if you use component in vcg/simplex/face/component_polygon.h
				  // It says "this class is a polygon and the memory for its components (e.g. pointer to its vertices
				  // will be allocated dynamically")
				  ,face::PFVAdj	 // Pointer to the vertices (just like FVAdj )
				  ,face::PFFAdj	 // Pointer to edge-adjacent face (just like FFAdj )
				  ,face::BitFlags // bit flags
				  ,face::Qualityf // quality
				  ,face::Normal3f // normal
				  
				  > {};

class PMesh: public
			  tri::TriMesh<
				  std::vector<PVertex>,	// the vector of vertices
				  std::vector<PFace >     // the vector of faces
				  >{};



#define FLAT	0
#define LOOP	1
#define CATMULL	2
#define BUTTERFLY 3
#define ONE_QUAD_X_EDGE 4
#define DOOSABIN 5

int  main(int argc, char **argv)
{
	if(argc<4)
	{
		printf(
			"\n                  PlyRefine (" __DATE__ ")\n"
			"						Visual Computing Group I.S.T.I. C.N.R.\n"
			"Usage: PlyRefine filein.ply fileout.[ply|off|obj|...] ref_step [opt] \n"
			"Commands: \n"
			" Refinement rules:\n"
			"     -m  use simple midpoint subdivision (default) \n"
			"     -b  use butterfly subdivision scheme \n"
			"     -l  use loop subdivision scheme \n"
			"     -o  use one-quad-per-edge schema (*) \n"
			"     -c  use Catmull-Clark (*) \n"
			"     -d  use DooSabin (*) \n"		 
			"     -e# refine only if the the edge is longer than #(default 0.0)\n"
			"Info:\n"
			"     (*) produces quad-only meshes, but updates topology only, \n"
			"         and leaves geometry unaffected \n"
			);
		exit(2);
	}
	
	int RefMode = FLAT	;
	int i=4; int n_steps; float length=0;
	while(i<argc)
	{
		if(argv[i][0]!='-')
		{printf("Error unable to parse option '%s'\n",argv[i]); exit(5);}
		switch(argv[i][1])
		{				
		case 'm' :	RefMode=FLAT; break;
		case 'b' :	RefMode=BUTTERFLY; break;
		case 'l' :	RefMode=LOOP; break;
		case 'c' :	RefMode=CATMULL; break;
		case 'o' :	RefMode=ONE_QUAD_X_EDGE; break;
		case 'd' :	RefMode=DOOSABIN; break;		
		case 'e' :	length=(float)atof(argv[i]+2); break;
		default : {printf("Error unable to parse option '%s'\n",argv[i]); exit(0);}
		}
		++i;
	}
	
	MyMesh m;
	
	if(tri::io::Importer<MyMesh>::Open(m,argv[1])!=0) {
		printf("Error reading file  %s\n",argv[1]);
		exit(1);
	}
	
	m.face.EnableFFAdjacency();
	tri::UpdateTopology<MyMesh>::FaceFace(m);
	tri::UpdateFlags<MyMesh>::FaceBorderFromFF(m);
	tri::UpdateNormal<MyMesh>::PerVertexNormalized(m);
	printf("Input mesh  vn:%i fn:%i\n",m.VN(),m.FN());
	
	n_steps=atoi(argv[3]);
	
	for(i=0;i < n_steps;++i)			
	{
		switch(RefMode)
		{
		case FLAT:
			tri::Refine<MyMesh, tri::MidPoint<MyMesh> >(m,tri::MidPoint<MyMesh>(&m),length);
			break;
		case LOOP:
			tri::RefineOddEven<MyMesh, tri::OddPointLoop<MyMesh>, tri::EvenPointLoop<MyMesh> >(m, tri::OddPointLoop<MyMesh>(m), tri::EvenPointLoop<MyMesh>(), length);
			break;
		case CATMULL:
			tri::BitQuadCreation<MyMesh>::MakePureByCatmullClark(m);
			tri::UpdateNormal<MyMesh>::PerBitQuadFaceNormalized(m);
			break;      
		case   ONE_QUAD_X_EDGE:
			tri::BitQuadCreation<MyMesh>::MakePureByRefine(m);
			tri::UpdateNormal<MyMesh>::PerBitQuadFaceNormalized(m);
			break;
		case BUTTERFLY:
			tri::Refine<MyMesh, tri::MidPointButterfly<MyMesh> >(m,tri::MidPointButterfly<MyMesh>(m),length);
			break;
		case DOOSABIN:
			PMesh base, refined;
			tri::UpdateTopology<MyMesh>::FaceFace(m);
			tri::UpdateNormal<MyMesh>::PerVertexNormalizedPerFaceNormalized(m);
			tri::PolygonSupport<MyMesh,PMesh>::MergeFlatFaces(m);
			tri::PolygonSupport<MyMesh,PMesh>::ImportFromTriMesh(base,m);
			tri::Clean<PMesh>::RemoveUnreferencedVertex(base);
			tri::Allocator<PMesh>::CompactEveryVector(base);
			
			tri::DooSabin<PMesh>::Refine(base, refined);
			tri::Allocator<PMesh>::CompactEveryVector(refined);
			tri::io::ExporterOBJ<PMesh>::Save(refined,"./refinedp.obj",tri::io::Mask::IOM_BITPOLYGONAL);
			tri::Allocator<PMesh>::CompactEveryVector(refined);
			tri::PolygonSupport<MyMesh,PMesh>::ImportFromPolyMesh(m,refined);
			PMesh testp;
			tri::UpdateTopology<MyMesh>::FaceFace(m);
			tri::PolygonSupport<MyMesh,PMesh>::ImportFromTriMesh(testp,m);
			tri::io::ExporterOBJ<PMesh>::Save(testp,"./refinedmp.obj",tri::io::Mask::IOM_BITPOLYGONAL);			
			tri::io::ExporterOBJ<MyMesh>::Save(m,"./refinedm.obj",tri::io::Mask::IOM_BITPOLYGONAL);
			break;
		}					
	}
	
	printf("Output mesh vn:%i fn:%i\n",m.VN(),m.FN());
	
	vcg::tri::io::Exporter<MyMesh>::Save(m,argv[2],vcg::tri::io::Mask::IOM_BITPOLYGONAL);
	return 0;
}
