/****************************************************************************
* VCGLib                                                            o o     *
* Visual and Computer Graphics Library                            o     o   *
*                                                                _   O  _   *
* Copyright(C) 2004                                                \/)\/    *
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


****************************************************************************/
#include <vector>
#include <string>  // anche questo
using namespace std;

#include<vcg/simplex/vertex/vertex.h>
#include<vcg/simplex/face/with/afav.h>
#include<vcg/simplex/face/pos.h>   // mi sembra di averlo aggiunto!


#include<vcg/complex/trimesh/base.h>
#include<vcg/complex/trimesh/update/topology.h>
#include <vcg/complex/trimesh/update/edges.h>
#include <vcg/complex/trimesh/update/bounding.h>
#include <vcg/math/histogram.h>
#include <vcg/complex/trimesh/clean.h>
#include <wrap/io_trimesh/import.h>
#include <wrap/io_trimesh/export_ply.h>

// loader 
#include<wrap/io_trimesh/import_ply.h>

using namespace vcg;


class MyFace;
class MyEdge;
class MyVertex:public Vertex<float,MyEdge,MyFace>{};
class MyFace :public FaceAFAV<MyVertex,MyEdge,MyFace>{};
class MyMesh: public tri::TriMesh< std::vector<MyVertex>, std::vector<MyFace > >{};

void OpenMesh(const char *filename, MyMesh &m)
{
  int err = tri::io::Importer<MyMesh>::Open(m,filename);
  if(err) {
      printf("Error in reading %s: '%s'\n",filename,tri::io::Importer<MyMesh>::ErrorMsg(err));
      exit(-1);
    }
  printf("read mesh `%s'\n", filename);		  
}



void main(int argc,char ** argv){

	MyMesh m;
	//load the mesh
	//argv[1]=(char*)"c:\\checkup\\debug\\column1m.ply";
	//argv[1] = "C:\\Documents and Settings\\Rita\\Desktop\\MeshReader\\prova2.ply";
	vcg::tri::io::ImporterPLY<MyMesh>::Open(m,argv[1]);
    FILE * index;
	index = fopen((string(argv[1])+string("2.html")).c_str(),"w");
	fprintf(index,"<p>Checkup: This is the check up result for %s </p>\n\n\n", argv[1]);
	
	fprintf(index,"<p>GENERAL INFO </p>\n\n");
	fprintf(index,"<p>Number of vertices: %d </p>\n", m.vn);
	fprintf(index,"<p>Number of faces: %d </p>\n", m.fn);
	if (m.Volume()!=0)
        fprintf(index,"<p>Volume: %d </p>\n", m.Volume());
	Color4b Color=m.C();
	fprintf(index, "<p>Object color(4b): %f %f %f </p>\n\n", Color[0], Color[1], Color[2]);

	
	/*fprintf(index,"BOUNDING BOX\n\n");
	if (m.bbox.IsEmpty())
		fprintf(index,"There's no info about bounding box/n/n");
	else
	{
		Point3f Center=m.bbox.Center();
		fprintf(index,"Bounding box center:     %f %f %f \n", Center.V(0),Center.V(1),Center.V(2));
		Point3f Dim=m.bbox.Dim();
		fprintf(index,"Bounding box dimensions: %f %f %f \n", abs(Dim.V(0)),abs(Dim.V(1)),abs(Dim.V(2)));
		fprintf(index,"Bounding box volume: %f \n", abs(m.bbox.Volume()));
		fprintf(index,"Bounding box diagonal: %f \n", m.bbox.Diag());
		fprintf(index,"Bounding box vertices are: %f %f %f \n", m.bbox.P(0).V(0), m.bbox.P(0).V(1), m.bbox.P(0).V(2));
		fprintf(index,"                           %f %f %f \n", m.bbox.P(1).V(0), m.bbox.P(1).V(1), m.bbox.P(1).V(2));
		fprintf(index,"                           %f %f %f \n", m.bbox.P(2).V(0), m.bbox.P(2).V(1), m.bbox.P(2).V(2));
		fprintf(index,"                           %f %f %f \n", m.bbox.P(3).V(0), m.bbox.P(3).V(1), m.bbox.P(3).V(2));
		fprintf(index,"                           %f %f %f \n", m.bbox.P(4).V(0), m.bbox.P(4).V(1), m.bbox.P(4).V(2));
		fprintf(index,"                           %f %f %f \n", m.bbox.P(5).V(0), m.bbox.P(5).V(1), m.bbox.P(5).V(2));
		fprintf(index,"                           %f %f %f \n", m.bbox.P(6).V(0), m.bbox.P(6).V(1), m.bbox.P(6).V(2));
		fprintf(index,"                           %f %f %f \n", m.bbox.P(7).V(0), m.bbox.P(7).V(1), m.bbox.P(7).V(2));

	}

	fprintf(index,"CAMERA INFO\n\n");
	if (m.camera.UberFlags()==NULL)
		fprintf(index, "No Camera Info\n\n");
	else
	{
	fprintf(index,"Flags: %s\n", m.camera.UberFlags());
	if (m.camera.IsOrtho())
		fprintf(index,"Ortogonal Camera:Yes\n");
	else
		fprintf(index,"Ortogonal Camera:No\n");
	fprintf(index,"Focal Distance: %f\n", m.camera.f);
	fprintf(index,"End of frustum: %f\n", m.camera.farend);
	fprintf(index,"Radial lens distortion coefficients: %f %f %f %f \n", m.camera.k[0],m.camera.k[1],m.camera.k[2],m.camera.k[3]);
	fprintf(index,"Pixel/size ratio: %f\n\n", m.camera.viewportM);
	}

	fprintf(index,"SHOT INFO\n\n");
	if (m.shot.IsValid())
		fprintf(index, "No Shot Info\n\n");

	fprintf(index,"VERTEX INFO\n\n");
	if (m.HasPerVertexNormal())
		fprintf(index, "Per Vertex Normal: YES\n");
	else
		fprintf(index, "Per Vertex Normal: NO\n");	
	if (m.HasPerVertexColor())
		fprintf(index, "Per Vertex Color: YES\n");
	else
		fprintf(index, "Per Vertex Color: NO\n");
	if (m.HasPerVertexMark())
		fprintf(index, "Per Vertex Mark: YES\n");
	else
		fprintf(index, "Per Vertex Mark: NO\n");
	if (m.HasPerVertexQuality())
		fprintf(index, "Per Vertex Quality: YES\n");
	else
		fprintf(index, "Per Vertex Quality: NO\n");
	if (m.HasPerVertexTexture())
		fprintf(index, "Per Vertex Texture: YES\n\n");
	else
		fprintf(index, "Per Vertex Texture: NO\n\n");

fprintf(index,"FACE INFO\n\n");
	if (m.HasPerFaceNormal())
		fprintf(index, "Per Face Normal: YES\n");
	else
		fprintf(index, "Per Face Normal: NO\n");	
	if (m.HasPerFaceColor())
		fprintf(index, "Per Face Color: YES\n");
	else
		fprintf(index, "Per Face Color: NO\n");
	if (m.HasPerFaceMark())
		fprintf(index, "Per Face Mark: YES\n");
	else
		fprintf(index, "Per Face Mark: NO\n");
	if (m.HasPerFaceQuality())
		fprintf(index, "Per Face Quality: YES\n\n");
	else
		fprintf(index, "Per Face Quality: NO\n\n");

fprintf(index,"WEDGE INFO\n\n");
	if (m.HasPerWedgeNormal())
		fprintf(index, "Per Wedge Normal: YES\n");
	else
		fprintf(index, "Per Wedge Normal: NO\n");	
	if (m.HasPerWedgeColor())
		fprintf(index, "Per Wedge Color: YES\n");
	else
		fprintf(index, "Per Wedge Color: NO\n");
	/*if (m.HasPerWedgeMark())
		fprintf(index, "Per Wedge Mark: YES\n");
	else
		fprintf(index, "Per Wedge Mark: NO\n");
	if (m.HasPerWedgeQuality())
		fprintf(index, "Per Wedge Quality: YES\n");
	else
		fprintf(index, "Per Wedge Quality: NO\n");
	if (m.HasPerWedgeTexture())
		fprintf(index, "Per Wedge Texture: YES\n\n");
	else
		fprintf(index, "Per Wedge Texture: NO\n\n");

	fprintf(index,"TOPOLOGY INFO\n\n");
	if (m.HasFFTopology())
		fprintf(index, "FFTopology: YES\n");
	else
	{
		fprintf(index, "FFTopology: NO\n");*/
		vcg::tri::UpdateTopology<MyMesh>::FaceFace(m);
	/*}
	if (m.HasVFTopology())
		fprintf(index, "VFTopology: YES\n\n");
	else
		fprintf(index, "VFTopology: NO\n\n");*/

	// COUNT EDGES

	MyMesh::FaceIterator fi;
	int count_e = 0;
	for(fi=m.face.begin();fi!=m.face.end();++fi)
		(*fi).ClearS();

	for(fi=m.face.begin();fi!=m.face.end();++fi)
		{
			(*fi).SetS();
			count_e +=3;
			for(int i=0; i<3; ++i)
				if((*fi).FFp(i)->IsS()) count_e--;
		}
	fprintf(index, "<p>Number of edges: %d </p>\n", count_e);

	// DA QUI IN POI!!!
	
	// DEGENERATED FACES

	int count_fd = 0;
	for(fi=m.face.begin(); fi!=m.face.end();++fi)
		if((*fi).Area() == 0)
			count_fd++;
	fprintf(index, "<p>Number of degenerated faces: %d </p>\n", count_fd);

	// UNREFERENCED VERTEX

	int count_uv = 0;
	MyMesh::FaceIterator f;
	MyMesh::VertexIterator v;
	int j;
	int deleted = 0;
	
	for(v=m.vert.begin();v!=m.vert.end();++v)
		(*v).ClearV();

	for(f=m.face.begin();f!=m.face.end();++f)
			for(j=0;j<3;++j)
					(*f).V(j)->SetV();

	for(v=m.vert.begin();v!=m.vert.end();++v)
		if( !(*v).IsV() )
			++count_uv;
	fprintf(index,"<p>Number of unreferenced vertices: %d</p>\n",count_uv);

// Holes count	

	for(f=m.face.begin();f!=m.face.end();++f)
		(*f).ClearS();
	MyMesh::FaceIterator g;
	g=m.face.begin(); f=g;
	int flag=0;
	int BEdges=0; int numholes=0;
	vcg::face::Pos<MyMesh::FaceType> hei;
	vcg::face::Pos<MyMesh::FaceType> he;

	while (f!=m.face.end())
	{
		flag=0;
		for(f=g;f!=m.face.end();++f)
		{
			for(j=0;j<3;j++)
			{
				if ((*f).IsBorder(j) && !(*f).IsS())
				{
					(*f).SetS();
					flag=1;
					g=f;
					break;
				}
			}
			if (flag==1)
				break;
		}
		numholes++;
		BEdges++;
		if (j>2) break;
		hei.Set(&(*g),j,g->V(j));
		he=hei;
	//		hei=he;
	//		do
	//			hei.Nextb()
	//		while(hei!=he);
		do
		{
			he.NextB();
			he.f->SetS();
			BEdges++;
		}
		while (he.f!=hei.f);
	}
	fprintf(index, "<p> Number of holes: %d </p> \n <p> Number of border edges: %d </p>", numholes, BEdges); 

}

