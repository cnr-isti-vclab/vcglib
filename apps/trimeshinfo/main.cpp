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
Revision 1.2  2005/01/03 16:13:09  rita_borgo
Added Standard comments



****************************************************************************/
#include <vector>
#include <string>  
#include <stack>
using namespace std;

#include<vcg/simplex/vertex/vertex.h>
#include<vcg/simplex/face/with/afav.h>
#include<vcg/simplex/face/topology.h>
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
using namespace face;

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
	argv[1] = "C:\\Documents and Settings\\Rita\\Desktop\\MeshReader\\trimeshinfo\\Debug\\prova0.ply";
	OpenMesh(argv[1],m);
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

	
	
		vcg::tri::UpdateTopology<MyMesh>::FaceFace(m);
	
	// IS MANIFOLD
	
	MyMesh::FaceIterator f;
	MyMesh::FaceIterator g;
	vcg::face::Pos<MyMesh::FaceType> he;
	vcg::face::Pos<MyMesh::FaceType> hei;
	int j;
	int man=0;
	bool Manifold = true;
	bool Manifold_lib = true;
	for(f=m.face.begin();f!=m.face.end();++f)
	{
		for (j=0;j<3;++j)
		{
			if(!IsManifold(*f,j))
			{
				Manifold_lib = false;
				f= m.face.end();
				break;
			}
		}
	}
	if (!Manifold_lib)
		fprintf(index, "<p> Manifold from lib gives: NO </p>"); 
	else
		fprintf(index, "<p> Manifold from lib gives: YES </p>"); 

	for(f=m.face.begin();f!=m.face.end();++f)
	{
		for (j=0;j<3;++j)
		{
      if ((*f).IsBorder(j))
			{}
			else if (&(*f) == (*f).FFp(j)->FFp((*f).FFi(j)))
			{}
			else
			{	
					hei.Set(&(*f), j , f->V(j));
					he=hei;
					he.NextF();
					while (he.f!=hei.f)
					{
						man++;
						he.NextF();
					}
					Manifold=false;
			}
		}
	}
	if (!Manifold)
		fprintf(index, "<p> Manifold from Matteo gives: NO </p>"); 
	else
		fprintf(index, "<p> Manifold from Matteo gives: YES </p>"); 


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
	MyMesh::VertexIterator v;
	
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
	g=m.face.begin(); f=g;
	
	int BEdges=0; int numholes=0;
	
	for(f=g;f!=m.face.end();++f)
	{
		if(!(*f).IsS())
		{
			for(j=0;j<3;j++)
			{
				if ((*f).IsBorder(j))
				{
					BEdges++;
					if(!(IsManifold(*f,j)))
					{
						(*f).SetS();
						hei.Set(&(*f),j,f->V(j));
						he=hei;
						do
						{
							he.NextB();
							he.f->SetS();
					//		BEdges++;
						}
						while (he.f!=hei.f);
						numholes++;
					}
				}
			}
		}
	}
	fprintf(index, "<p> Number of holes: %d </p> \n <p> Number of border edges: %d </p>", numholes, BEdges); 

	// CONNECTED COMPONENTS


	for(f=m.face.begin();f!=m.face.end();++f)
		(*f).ClearS();
	g=m.face.begin(); f=g;
	int CountComp=0; int CountOrient=0;
	stack<MyMesh::FaceIterator> sf;	
	MyMesh::FaceType *l;
	for(f=m.face.begin();f!=m.face.end();++f)
	{
		if (!(*f).IsS())
		{
			(*f).SetS();
			sf.push(f);
			while (!sf.empty())
			{
				g=sf.top();
				he.Set(&(*g),0,g->V(0));
				sf.pop();
				for(j=0;j<3;++j)
						if( !(*g).IsBorder(j) )
							{
								l=he.f->FFp(j);
								if( !(*l).IsS() )
									{
										(*l).SetS();
										sf.push(l);
									}
							}
			}
		CountComp++;
		}
	}
	fprintf(index, "<p> Number of connected components: %d </p>", CountComp); 
	
// ORIENTABLE E ORIENTED MESH

	int flag=0;
	bool Oriented=true;
	if (!Manifold)
		fprintf(index, "<p> Orientable Mesh: NO</p>"); 
	else
	{
		for(f=m.face.begin();f!=m.face.end();++f)
		{
			(*f).ClearS();
		//	(*f).ClearR();
		}
		g=m.face.begin(); f=g; 
		for(f=m.face.begin();f!=m.face.end();++f)
		{
			if (!(*f).IsS())
			{
				(*f).SetS();
				sf.push(f);
				
				while (!sf.empty())
				{
					g=sf.top();
					sf.pop();
					for(j=0;j<3;++j)
					{
						int prova = (*g).IsR();
						if( !(*g).IsBorder(j) )
						{
							he.Set(&(*g),0,g->V(0));
							l=he.f->FFp(j);
							if( !(*l).IsS() )
							{
								(*l).SetS();
								sf.push(l);
							}
							he.Set(&(*g),j,g->V(j));								
							hei.Set(he.f->FFp(j),he.f->FFi(j), (he.f->FFp(j))->V(he.f->FFi(j)));
							if (he.v!=hei.v)    // bene
							{
								if ((*l).IsS())
								{}
								else
								{
									(*l).SetS();
									sf.push(l);
								}
							}	
						}
						else if ((*l).IsS() && !(*l).IsR())
						{
							flag=1;
							break;
						}
						else
						{
							Oriented=false;
							(*l).SetS();
							(*l).SetR();
							sf.push(l);
						}
					}
				}
			}
			if (flag==1)
				break;
		}
		if (flag==0)
				fprintf(index, "<p> Orientable Mesh: YES</p>"); 
		else
				fprintf(index, "<p> Orientable Mesh: NO</p>"); 
	}
	if (Oriented && Manifold)
			fprintf(index, "<p> Oriented Mesh: YES</p>"); 
	else
			fprintf(index, "<p> Oriented Mesh: NO</p>"); 

	fclose(index);
}

