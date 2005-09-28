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
Revision 1.2  2005/09/28 10:46:04  rita_borgo
Added possibility of saving File in OFF format

Revision 1.1  2005/09/20 10:15:27  rita_borgo
Changed file name to uniform with other solution projects,
 before was main.cpp

Revision 1.8  2005/02/15 12:26:06  rita_borgo
Minor changes to self-intersection

Revision 1.7  2005/02/07 15:44:31  rita_borgo
Fixed Color and Volume

Revision 1.6  2005/02/01 17:37:53  rita_borgo
Fixed Volume and Color

Revision 1.5  2005/01/18 16:33:12  rita_borgo
Added OFF file Option

Revision 1.4  2005/01/17 18:19:00  rita_borgo
Added new routines.
Self-intersection first release

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

#include "defs.h"

using namespace vcg;
using namespace tri;
using namespace face;

class MyFace;
class MyEdge;
class MyVertex:public Vertex<float,MyEdge,MyFace>{};
class MyFace :public FaceAFAV<MyVertex,MyEdge,MyFace>{};
class MyMesh: public tri::TriMesh< std::vector<MyVertex>, std::vector<MyFace > >{};


typedef MyMesh::VertexPointer  VertexPointer;
typedef MyMesh::VertexIterator  VertexIterator;

typedef Point3<MyMesh::ScalarType> Point3x;
typedef vector<Point3x> Hole;

string ans;

void OpenMesh(const char *filename, MyMesh &m)
{
  int err = tri::io::Importer<MyMesh>::Open(m,filename);
  if(err) {
      printf("Error in reading %s: '%s'\n",filename,tri::io::Importer<MyMesh>::ErrorMsg(err));
      exit(-1);
    }
  printf("read mesh `%s'\n", filename);		  
}


inline char* GetExtension(char* filename)
{
    for(int i=strlen(filename)-1; i >= 0; i--)
        if(filename[i] == '.')
            break;
    if(i > 0)
        return &(filename[i+1]);
    else
        return NULL;
}


typedef MyMesh::VertexPointer  VertexPointer;
typedef MyMesh::VertexIterator  VertexIterator;

/* classe di confronto per l'algoritmo di individuazione vertici duplicati*/
template <class VertexIterator>
class DuplicateVert_Compare{
public:
	inline bool operator() (VertexIterator a, VertexIterator b)
		{
			return *a < *b;
		}
};

static int DuplicateVertex( MyMesh & m )    // V1.0
{
	
	if(m.vert.size()==0 || m.vn==0)
		return 0;
	std::map<VertexPointer, VertexPointer> mp;
	int i,j;
	VertexIterator vi; 
	int deleted=0;
	int k=0;
	int num_vert = m.vert.size();
	vector<VertexPointer> perm(num_vert);
	for(vi=m.vert.begin(); vi!=m.vert.end(); ++vi, ++k)
		perm[k] = &(*vi);

	DuplicateVert_Compare<VertexPointer> c_obj;

	std::sort(perm.begin(),perm.end(),c_obj);
	j = 0;
  i = j;
  mp[perm[i]] = perm[j];
  ++i;
  for(;i!=num_vert;)
	{
		if( (! (*perm[i]).IsD()) && 
        (! (*perm[j]).IsD()) && 
				(*perm[i]).P() == (*perm[j]).cP() )
		{
			if(deleted ==0)
			{
				cout<< "\t Found Duplicated Vertices"<< endl;
				cout<< "\t do you want to remove them? [y/Y| n/N]\n";
				cin>>ans;
			}
			if((ans == "s") ||(ans=="S"))
			{
				VertexPointer t = perm[i];
				mp[perm[i]] = perm[j];
				(*t).SetD();
				m.vn--;
			}
			++i;
			deleted++;// per evitare che su risposta no rifaccia la domanda ad ogni vertice che incontra;
		}
		else
		{
			j = i;
	    ++i;
		}
	}
	ans.clear();
	return deleted;
}
void main(int argc,char ** argv){

	
	MyMesh m;
	bool DEBUG = false;
	

/*------------XML file part ------------------*/

	static char* XML_SCHEMA_NAME = "protegekb";
	XMLTree	doc;
	MainNode* mn = new MainNode;

/*--------------------------------------------*/

	//load the mesh
	//argv[1]=(char*)"c:\\checkup\\debug\\column1m.ply";
	//argv[1] = "C:\\sf\\apps\\msvc\\trimeshinfo\\Release\\prism.off";
//argv[1] = "C:\\sf\\apps\\msvc\\trimeshinfo\\Release\\prova1.ply";

    // print program info
    printf("-------------------------------\n"
           "         TriMeshInfo V.1.01 \n"
           "     http://vcg.isti.cnr.it\n"
           "   release date: "__DATE__"\n"
           "-------------------------------\n\n");


 if(DEBUG)
	argv[1] = "C:\\sf\\apps\\msvc\\trimeshinfo\\cube.ply";
 
 else
 {
 // load input meshes.
  if(argc <= 1)
  {
     printf(MSG_ERR_N_ARGS);
     exit(-1);
  }
 }


  OpenMesh(argv[1],m);

	

	
	doc.initializeMain(XML_SCHEMA_NAME);
	char* s1 = "http://www.w3.org/2001/XMLSchema-instance";
	char* s2 = new(char[50]);
	sprintf(s2,"\"%s\"",s1);	
	doc.addHeaders(" xmlns:xsi=", s2);
	s1 = "http://protege.stanford.edu/plugins/xmlbackend/protege_xml_backend.xsd";
	s2 = new(char[100]);
	sprintf(s2,"\"%s\"",s1);	
	doc.addHeaders(" xsi:noNamespaceSchemaLocation=", s2);

  FILE * index;
	index = fopen((string(argv[1])+string("2.html")).c_str(),"w");
	fprintf(index,"<p>Mesh info: %s </p>\n\n\n", argv[1]);
	
	fprintf(index,"<p>GENERAL INFO </p>\n\n");
	fprintf(index,"<p>Number of vertices: %d </p>\n", m.vn);
	fprintf(index,"<p>Number of faces: %d </p>\n", m.fn);
	printf("Mesh info:\n");
	printf("	M: '%s'\n\t Number of vertices: %d \n", argv[1], m.vn);
	printf("\t Number of faces: %d \n", m.fn);

	//XML
	
	NodeGroup* ng = new NodeGroup;

	char* s =new(char[25]);
	sprintf(s,"%d",m.vn);	
	ValueNode* vn = new ValueNode;
	EntryNode* en = new EntryNode;
	vn->setValue(s);
	en->addValue(*vn);
	en->type = "Integer";
	OwnSlotNode* osn = new OwnSlotNode;
	osn->setName("Number of Vertices");
	osn->addEntry(*en);
	SlotNode* sn = new SlotNode;
	sn->addOwnSlot(osn);
	ng->addNode(osn);
	
	s = new(char[25]);
	vn = new ValueNode;
	en = new EntryNode;
	osn = new OwnSlotNode;
	sprintf(s,"%d",m.fn);	
	vn->setValue(s);
	en->addValue(*vn);
	en->type = "Integer";

	osn->setName("Number of Faces");
	osn->addEntry(*en);
	sn->addOwnSlot(osn);
	ng->addNode(osn);

	

	if(m.HasPerFaceColor()||m.HasPerVertexColor())
	{
		Color4b Color=m.C();
		fprintf(index, "<p>Object color(4b): %f %f %f </p>\n\n", Color[0], Color[1], Color[2]);
		printf( "\t Object color(4b): %f %f %f \n", Color[0], Color[1], Color[2]);

		s = new(char[25]);
		vn = new ValueNode;
		en = new EntryNode;
		osn = new OwnSlotNode;
		sprintf(s,"%f %f %f ",Color[0], Color[1], Color[2]);	
		vn->setValue(s);
		en->addValue(*vn);
		en->type = "Float";

		osn->setName("Colors");
		osn->addEntry(*en);
		sn->addOwnSlot(osn);
		ng->addNode(osn);
	}
	
	


	
		vcg::tri::UpdateTopology<MyMesh>::FaceFace(m);

	// IS MANIFOLD
	MyMesh::FaceIterator fi;
	MyMesh::FaceIterator gi;
	vcg::face::Pos<MyMesh::FaceType> he;
	vcg::face::Pos<MyMesh::FaceType> hei;
	int j;
	for(fi=m.face.begin();fi!=m.face.end();++fi)
		(*fi).ClearS();
	
	int man=0;
	bool Manifold = true;
	
	
	for(fi=m.face.begin();fi!=m.face.end();++fi)
	{
		for (j=0;j<3;++j)
		{
			if(!IsManifold(*fi,j))
			{
				Manifold = false;
				fi= m.face.end();
				--fi;
				j=3;
			}
		}
	}
	if (!Manifold)
	{
		fprintf(index, "<p> Manifold: NO </p>"); 
	  printf( "\t Manifold: NO\n"); 
		s = new(char[25]);
		vn = new ValueNode;
		en = new EntryNode;
		osn = new OwnSlotNode;
		s = "No";
		vn->setValue(s);
		en->addValue(*vn);
		en->type = "Bool";

		osn->setName("Manifold");
		osn->addEntry(*en);
		sn->addOwnSlot(osn);
		ng->addNode(osn);
	}
	else
	{
		fprintf(index, "<p> Manifold: YES </p>"); 
	  printf( "\t Manifold: YES\n "); 
		s = new(char[25]);
		vn = new ValueNode;
		en = new EntryNode;
		osn = new OwnSlotNode;
		s = "Yes";
		vn->setValue(s);
		en->addValue(*vn);
		en->type = "Bool";

		osn->setName("Manifold");
		osn->addEntry(*en);
		sn->addOwnSlot(osn);
		ng->addNode(osn);
	}

	
	// COUNT EDGES
	int count_e = 0;
	int boundary_e = 0;
	bool counted =false;
	for(fi=m.face.begin();fi!=m.face.end();fi++)
		{
			(*fi).SetS();
			count_e +=3;								//assume that we have to increase the number of edges with three
			for(int j=0; j<3; j++)
			{
				if (fi->IsBorder(j))			//If this edge is a border edge
					boundary_e++;						//  then increase the number of boundary edges
				else if (IsManifold(*fi,j))		//If this edge is manifold
				{
					if((*fi).FFp(j)->IsS()) //If the face on the other side of the edge is already selected
						count_e--;						//  we counted one edge twice
				}
				else											//We have a non-manifold edge
				{
					hei.Set(&(*fi), j , fi->V(j));
					he=hei;
					he.NextF();
					while (he.f!=hei.f)			//	so we have to iterated all faces that are connected to this edge
					{
						if (he.f->IsS())			//  if one of the other faces was already visited than this edge was counted already.
						{
							counted=true;
							break;
						}
						else 
						{
							he.NextF();
						}
					}
					if (counted)
					{
						count_e--;
						counted=false;
					}
				}
			}
		}	
		fprintf(index, "<p>Number of edges: %d </p>\n", count_e);
		fprintf(index, "<p>Number of internal edges: %d </p>\n", count_e-boundary_e);
		fprintf(index, "<p>Number of boundary edges: %i </p>\n", boundary_e);
		printf("\t Number of edges: %d \n", count_e);
		printf("\t Number of internal edges: %d \n", count_e-boundary_e);
		printf("\t Number of boundary edges: %i \n", boundary_e);

		s = new(char[25]);
		vn = new ValueNode;
		en = new EntryNode;
		osn = new OwnSlotNode;
		sprintf(s,"%d",count_e);	
		vn->setValue(s);
		en->addValue(*vn);
		en->type = "Integer";

		osn->setName("Number of Edges");
		osn->addEntry(*en);
		sn->addOwnSlot(osn);
		ng->addNode(osn);


	
	// DEGENERATED FACES

	int count_fd = 0;
	for(fi=m.face.begin(); fi!=m.face.end();++fi)
		if((*fi).Area() == 0)
			count_fd++;
	fprintf(index, "<p>Number of degenerated faces: %d </p>\n", count_fd);
  printf("\t Number of degenerated faces: %d \n", count_fd);
		s = new(char[25]);
		vn = new ValueNode;
		en = new EntryNode;
		osn = new OwnSlotNode;
		sprintf(s,"%d",count_fd);	
		vn->setValue(s);
		en->addValue(*vn);
		en->type = "Integer";

		osn->setName("Number of Degenerated Faces");
		osn->addEntry(*en);
		sn->addOwnSlot(osn);
		ng->addNode(osn);

	// UNREFERENCED VERTEX

	int count_uv = 0;
	MyMesh::VertexIterator v;
	
	
	
	for(v=m.vert.begin();v!=m.vert.end();++v)
		(*v).ClearV();

	for(fi=m.face.begin();fi!=m.face.end();++fi)
			for(j=0;j<3;++j)
					(*fi).V(j)->SetV();

	for(v=m.vert.begin();v!=m.vert.end();++v)
		if( !(*v).IsV() )
			++count_uv;
	fprintf(index,"<p>Number of unreferenced vertices: %d</p>\n",count_uv);
  printf("\t Number of unreferenced vertices: %d\n",count_uv);
		s = new(char[25]);
		vn = new ValueNode;
		en = new EntryNode;
		osn = new OwnSlotNode;
		sprintf(s,"%d",count_uv);	
		vn->setValue(s);
		en->addValue(*vn);
		en->type = "Integer";

		osn->setName("Number of unreferenced vertices");
		osn->addEntry(*en);
		sn->addOwnSlot(osn);
		ng->addNode(osn);	


// HOLES COUNT	

	int numholes=0;
	int numholev=0;
	int BEdges=0; 
	vector<vector<Point3x> > holes; //indices of vertices

	for(fi=m.face.begin();fi!=m.face.end();++fi)
		(*fi).ClearS();
	gi=m.face.begin(); fi=gi;
	
	
		if (Manifold)
	{
    for(fi=m.face.begin();fi!=m.face.end();fi++)					//for all faces do
		{
				for(j=0;j<3;j++)										//for all edges
				{
					if(fi->V(j)->IsS()) continue;

					if((*fi).IsBorder(j))							//found an unvisited border edge
					{
						he.Set(&(*fi),j,fi->V(j));			//set the face-face iterator to the current face, edge and vertex
						vector<Point3x> hole;						//start of a new hole
						hole.push_back(fi->P(j));				//  including the first vertex
						numholev++;
						he.v->SetS();										//set the current vertex as selected
						he.NextB();											//go to the next boundary edge
						
						
						while(fi->V(j) != he.v)					//will we do not encounter the first boundary edge.
						{
							Point3x newpoint = he.v->P();		//select its vertex.
							if(he.v->IsS())									//check if this vertex was selected already, because then we have an additional hole.
							{
								//cut and paste the additional hole.
								vector<Point3x> hole2;				
								int index = find(hole.begin(),hole.end(),newpoint) - hole.begin();
								for(int i=index; i<hole.size(); i++)
									hole2.push_back(hole[i]);

								hole.resize(index);
								if(hole2.size()!=0)						//annoying in degenerate cases
									holes.push_back(hole2);							
							}
							hole.push_back(newpoint);
							numholev++;
							he.v->SetS();										//set the current vertex as selected
							he.NextB();											//go to the next boundary edge
						}
						holes.push_back(hole);
					}
				}
		}
		numholes = holes.size();
    fprintf(index,"<p>Number of holes/boundaries: %d </p>\n", numholes); 
    printf("\t Number of holes/boundaries: %d \n", numholes); 
		s = new(char[25]);
		vn = new ValueNode;
		en = new EntryNode;
		osn = new OwnSlotNode;
		sprintf(s,"%d",numholes);	
		vn->setValue(s);
		en->addValue(*vn);
		en->type = "Integer";

		osn->setName("Number of Holes");
		osn->addEntry(*en);
		sn->addOwnSlot(osn);
		ng->addNode(osn);	

		

		if(numholes>0)
		{
			int BEdges = 0;
			printf("\t Edges per hole/boundary:\n\t (");
			for(int i=0; i<numholes; i++)
			{
				if(i==numholes-1){ printf("%i)\n",holes[i].size()); BEdges++;}
				else{ printf("%i, ",holes[i].size()); BEdges++;}
			}
			s = new(char[25]);
			vn = new ValueNode;
			en = new EntryNode;
			osn = new OwnSlotNode;
			sprintf(s,"%d",BEdges);	
			vn->setValue(s);
			en->addValue(*vn);
			en->type = "Integer";
			osn->setName("Number of Border Edges");
			osn->addEntry(*en);
			sn->addOwnSlot(osn);
			ng->addNode(osn);		
		}
	}
	else
		printf( "\t Number of holes: UNDEFINED, mesh is non-manifold \n");

	
	// Mesh Volume
	float vol = m.Volume();
	int nuh = numholes;
	if((m.Volume()>0.)&&(Manifold)&&(numholes==0))
	{
        fprintf(index,"<p>Volume: %d </p>\n", m.Volume());
        printf("\t Volume: %f \n", m.Volume());
		s = new(char[25]);
		vn = new ValueNode;
		en = new EntryNode;
		osn = new OwnSlotNode;
		sprintf(s,"%f",m.Volume());	
		vn->setValue(s);
		en->addValue(*vn);
		en->type = "Float";

		osn->setName("Volume");
		osn->addEntry(*en);
		sn->addOwnSlot(osn);
		ng->addNode(osn);		
	}
	else 
	{
		printf("\t Volume: UNDEFINED, mesh is either non-manifold or has holes \n");
    fprintf(index,"Volume: UNDEFINED, mesh is either non-manifold or has holes \n");
	}


	// CONNECTED COMPONENTS

	vector<int> nrfaces;
	for(fi=m.face.begin();fi!=m.face.end();++fi)
		(*fi).ClearS();
	gi=m.face.begin(); fi=gi;
	int Compindex=0;
	stack<MyMesh::FaceIterator> sf;	
	MyMesh::FaceType *l;
	for(fi=m.face.begin();fi!=m.face.end();++fi)
	{
		if (!(*fi).IsS())
		{
			(*fi).SetS();
			(*fi).Q()=Compindex;
			nrfaces.push_back(1);
			sf.push(fi);
			while (!sf.empty())
			{
				gi=sf.top();
				he.Set(&(*gi),0,gi->V(0));
				sf.pop();
				for(j=0;j<3;++j)
				{
						if( !(*gi).IsBorder(j) )
							{
								l=he.f->FFp(j);
								if( !(*l).IsS() )
									{
										(*l).SetS();
										sf.push(l);
									}
							}
				}
			}
		Compindex++;
		}
	}
	int numcomponents = nrfaces.size();
	fprintf(index, "<p> Number of connected components: %d </p>", numcomponents); 
  printf("\t Number of connected components: %d\n", numcomponents); 
	s = new(char[25]);
		vn = new ValueNode;
		en = new EntryNode;
		osn = new OwnSlotNode;
		sprintf(s,"%d",numcomponents);	
		vn->setValue(s);
		en->addValue(*vn);
		en->type = "Integer";

		osn->setName("Number of Connected Components");
		osn->addEntry(*en);
		sn->addOwnSlot(osn);
		ng->addNode(osn);

	//GENUS  --> 2( #components - genus ) = #vertices + #faces - #edge - #boundary_loops = eulernumber - #holes
	//eulero = (mesh.vn-count_uv) - (count_e)+mesh.fn;

	int eulernumber = (m.vn-count_uv) + m.fn - count_e;
	if(Manifold)
	{
		int genus = -( 0.5 * (eulernumber - numholes) - numcomponents );
		fprintf(index, "<p> Genus: %d </p> \n ", genus); 
		printf( "\t Genus: %d \n", genus);
		s = new(char[25]);
		vn = new ValueNode;
		en = new EntryNode;
		osn = new OwnSlotNode;
		sprintf(s,"%d",genus);	
		vn->setValue(s);
		en->addValue(*vn);
		en->type = "Integer";

		osn->setName("Genus");
		osn->addEntry(*en);
		sn->addOwnSlot(osn);
		ng->addNode(osn);			
	}
	else //(!Manifold) 
	{
		fprintf( index,"<p>Genus: UNDEFINED, mesh is non-manifold </p>\n");
		printf( "Genus: UNDEFINED, mesh is non-manifold \n");
	}
// REGULARITY

	bool Regular=true;
	bool Semiregular=true;
	int inc=0;
	for(v=m.vert.begin();v!=m.vert.end();++v)
		(*v).ClearS();
	for(fi=m.face.begin();fi!=m.face.end();++fi)
	{
		for (j=0; j<3; j++)
		{
			he.Set(&(*fi),j,fi->V(j));
			if (!(*fi).IsBorder(j) && !(*fi).IsBorder((j+2)%3) && !fi->V(j)->IsS())
			{
				hei=he;
				inc=1;
				he.FlipE();
				he.NextF();
				while (he.f!=hei.f)
				{
					he.FlipE();
					if (he.IsBorder())
					{
						inc=6;
						break;
					}
					he.NextF();
					inc++;
				}
				if (inc!=6)
					Regular=false;
				if (inc!=6 && inc!=5)
					Semiregular=false;
				fi->V(j)->SetS();

			}
			else
				fi->V(j)->SetS();
		}
		if (Semiregular==false)
			break;

	}

	if (Regular)
	{
			fprintf(index, "<p> Type of Mesh: REGULAR</p>"); 
		  printf("\t Type of Mesh: REGULAR\n"); 
			s = new(char[25]);
		vn = new ValueNode;
		en = new EntryNode;
		osn = new OwnSlotNode;
		s="REGULAR";	
		vn->setValue(s);
		en->addValue(*vn);
		en->type = "String";

		osn->setName("Type of Mesh");
		osn->addEntry(*en);
		sn->addOwnSlot(osn);
		ng->addNode(osn);		
	}
	else if (Semiregular)
	{
			fprintf(index, "<p> Type of Mesh: SEMIREGULAR</p>");
		  printf("\t Type of Mesh: SEMIREGULAR\n");
			s = new(char[25]);
		vn = new ValueNode;
		en = new EntryNode;
		osn = new OwnSlotNode;
		s="SEMIREGULAR";	
		vn->setValue(s);
		en->addValue(*vn);
		en->type = "String";

		osn->setName("Type of Mesh");
		osn->addEntry(*en);
		sn->addOwnSlot(osn);
		ng->addNode(osn);		
	}
	else 
	{
		fprintf(index, "<p> Type of Mesh: IRREGULAR</p>"); 
	  printf("\t Type of Mesh: IRREGULAR\n"); 
			s = new(char[25]);
		vn = new ValueNode;
		en = new EntryNode;
		osn = new OwnSlotNode;
		s="IRREGULAR";	
		vn->setValue(s);
		en->addValue(*vn);
		en->type = "String";

		osn->setName("Type of Mesh");
		osn->addEntry(*en);
		sn->addOwnSlot(osn);
		ng->addNode(osn);		
	}
// ORIENTABLE E ORIENTED MESH

	bool Orientable=true;
	bool Oriented=true;
	if (!Manifold)
	{
		fprintf(index, "<p> Orientable Mesh: NO</p>"); 
	  printf( "\t Orientable Mesh: NO\n"); 
			s = new(char[25]);
		vn = new ValueNode;
		en = new EntryNode;
		osn = new OwnSlotNode;
		s="NO";	
		vn->setValue(s);
		en->addValue(*vn);
		en->type = "String";

		osn->setName("Orientable Mesh");
		osn->addEntry(*en);
		sn->addOwnSlot(osn);
		ng->addNode(osn);		
	}
	else
	{
		for(fi=m.face.begin();fi!=m.face.end();++fi)
		{
			(*fi).ClearS();
			(*fi).ClearUserBit(0);
		}
		gi=m.face.begin(); fi=gi; 
		for(fi=m.face.begin();fi!=m.face.end();++fi)
		{
			if (!(*fi).IsS())
			{
				(*fi).SetS();
				sf.push(fi);
				
				while (!sf.empty())
				{
					gi=sf.top();
					sf.pop();
					for(j=0;j<3;++j)
					{
						if( !(*gi).IsBorder(j) )
						{
							he.Set(&(*gi),0,gi->V(0));
							l=he.f->FFp(j);
							he.Set(&(*gi),j,gi->V(j));								
							hei.Set(he.f->FFp(j),he.f->FFi(j), (he.f->FFp(j))->V(he.f->FFi(j)));
							if( !(*gi).IsUserBit(0) )
							{
								if (he.v!=hei.v)    // bene
								{
									if ((*l).IsS() && (*l).IsUserBit(0))
									{
										Orientable=false;
										break;
									}
									else if (!(*l).IsS())
									{
										(*l).SetS();
										sf.push(l);
									}
								}	
								else if (!(*l).IsS())
								{
									Oriented=false;
									(*l).SetS();
									(*l).SetUserBit(0);
									sf.push(l);
								}
								else if ((*l).IsS() && !(*l).IsUserBit(0))
								{
									Orientable=false;
									break;
								}
							}
							else if (he.v==hei.v)    // bene
							{
								if ((*l).IsS() && (*l).IsUserBit(0))
								{
									Orientable=false;
									break;
								}
								else if (!(*l).IsS())
								{
									(*l).SetS();
									sf.push(l);
								}
							}	
							else if (!(*l).IsS())
							{
								Oriented=false;
								(*l).SetS();
								(*l).SetUserBit(0);
								sf.push(l);
							}
							else if ((*l).IsS() && !(*l).IsUserBit(0))
							{
								Orientable=false;
								break;
							}
						}
					}
				}
			}
			if (!Orientable)
				break;
		}
		if (Orientable)
		{
				fprintf(index, "<p> Orientable Mesh: YES</p>"); 
			  printf( "\t Orientable Mesh: YES\n");
							s = new(char[25]);
		vn = new ValueNode;
		en = new EntryNode;
		osn = new OwnSlotNode;
		s="YES";	
		vn->setValue(s);
		en->addValue(*vn);
		en->type = "String";

		osn->setName("Orientable Mesh");
		osn->addEntry(*en);
		sn->addOwnSlot(osn);
		ng->addNode(osn);		
		}
		else
		{
				fprintf(index, "<p> Orientable Mesh: NO</p>"); 
			  printf( "\t Orientable Mesh: NO\n"); 
							s = new(char[25]);
		vn = new ValueNode;
		en = new EntryNode;
		osn = new OwnSlotNode;
		s="NO";
		vn->setValue(s);
		en->addValue(*vn);
		en->type = "String";

		osn->setName("Orientable Mesh");
		osn->addEntry(*en);
		sn->addOwnSlot(osn);
		ng->addNode(osn);		
		}
	}
	if (Oriented && Manifold)
	{
			fprintf(index, "<p> Oriented Mesh: YES</p>"); 
		  printf( "\t Oriented Mesh: YES\n"); 
			s = new(char[25]);
		vn = new ValueNode;
		en = new EntryNode;
		osn = new OwnSlotNode;
		s="YES";	
		vn->setValue(s);
		en->addValue(*vn);
		en->type = "String";

		osn->setName("Oriented Mesh");
		osn->addEntry(*en);
		sn->addOwnSlot(osn);
		ng->addNode(osn);			
	}
	else
	{
			fprintf(index, "<p> Oriented Mesh: NO</p>"); 
		  printf( "\t Oriented Mesh: NO\n"); 
			s = new(char[25]);
		vn = new ValueNode;
		en = new EntryNode;
		osn = new OwnSlotNode;
		s="NO";	
		vn->setValue(s);
		en->addValue(*vn);
		en->type = "String";

		osn->setName("Oriented Mesh");
		osn->addEntry(*en);
		sn->addOwnSlot(osn);
		ng->addNode(osn);			
	}
	int dv = DuplicateVertex(m);
	if(dv>0)
	{
		fprintf(index, "<p> Duplicated vertices: %d</p>", dv); 
		printf( "\t Number of duplicated vertices found: %d\n",dv);		
		s = new(char[25]);
		vn = new ValueNode;
		en = new EntryNode;
		osn = new OwnSlotNode;
		sprintf(s,"%d",dv);	
		vn->setValue(s);
		en->addValue(*vn);
		en->type = "Int";

		osn->setName("Duplicated Vertices");
		osn->addEntry(*en);
		sn->addOwnSlot(osn);
		ng->addNode(osn);		
	}
	else
	{
		fprintf(index, "<p> Duplicated vertices: NO</p>"); 
		printf( "\t Duplicated vertices: NO\n");

			s = new(char[25]);
		vn = new ValueNode;
		en = new EntryNode;
		osn = new OwnSlotNode;
		s="NO";	
		vn->setValue(s);
		en->addValue(*vn);
		en->type = "String";

		osn->setName("Duplicated Vertices");
		osn->addEntry(*en);
		sn->addOwnSlot(osn);
		ng->addNode(osn);	
	}
	// SELF INTERSECTION
	if (m.fn<300000)
	{
		bool SelfInt=false;
		for(fi=m.face.begin();fi!=m.face.end();++fi)
		{
			for(gi=++fi , fi--;gi!=m.face.end();++gi)
			{
				if ((*fi).FFp(0)!=&(*gi) && (*fi).FFp(1)!=&(*gi) && (*fi).FFp(2)!=&(*gi) &&
					fi->V(0)!=gi->V(0) && fi->V(0)!=gi->V(1) && fi->V(0)!=gi->V(2) &&
					fi->V(1)!=gi->V(0) && fi->V(1)!=gi->V(1) && fi->V(1)!=gi->V(2) &&
					fi->V(2)!=gi->V(0) && fi->V(2)!=gi->V(1) && fi->V(2)!=gi->V(2))
				{
					if (NoDivTriTriIsect(fi->V(0)->P(), fi->V(1)->P(), fi->V(2)->P(),gi->V(0)->P(), gi->V(1)->P(), gi->V(2)->P()) )
						SelfInt=true;
				}
			}
			if (SelfInt)
				break;			
		}
		if (SelfInt)
		{
			fprintf(index, "<p> Self Intersection: YES</p>"); 
			printf( "\t Self Intersection: YES\n");
			s = new(char[25]);
		vn = new ValueNode;
		en = new EntryNode;
		osn = new OwnSlotNode;
		s="YES";	
		vn->setValue(s);
		en->addValue(*vn);
		en->type = "String";

		osn->setName("Self Intersection");
		osn->addEntry(*en);
		sn->addOwnSlot(osn);
		ng->addNode(osn);	
		}
		else
		{
			fprintf(index, "<p> Self Intersection: NO</p>"); 
			 printf( "\t Self Intersection: NO\n"); 
			 s = new(char[25]);
		vn = new ValueNode;
		en = new EntryNode;
		osn = new OwnSlotNode;
		s="NO";	
		vn->setValue(s);
		en->addValue(*vn);
		en->type = "String";

		osn->setName("Self Intersection");
		osn->addEntry(*en);
		sn->addOwnSlot(osn);
		ng->addNode(osn);	
		}
	}



	string fs;
	
	cout<< "\t To save the file: [s/S]"<< endl;
	cin>>ans;

		
		if((ans == "S")||(ans == "s"))
		{
			cout<< "\t available formats: [ply, off, stl] "<<endl;
			cout<< "\t enter format"<<endl;
			cin>>ans;
			cout<<"\t enter filename"<<endl;
			cin>>fs;
			const char* filesave = fs.c_str();
			if(ans == "ply")
				tri::io::ExporterPLY<MyMesh>::Save(m, filesave);
			else if(ans == "stl")
				tri::io::ExporterSTL<MyMesh>::Save(m,filesave);
			else if(ans == "dxf")
				tri::io::ExporterDXF<MyMesh>::Save(m,filesave);
	
			else if(ans == "off")
				tri::io::ExporterOFF<MyMesh>::Save(m,filesave);
		}

	doc.addSlots(sn);
	OwnSlotsNode* ossn = new OwnSlotsNode;
	ossn->addOwnSlot(ng);
	ClassNode* cn = new ClassNode;
	cn->addOwnSlots(ossn);
	doc.addClasses(cn);

	InstanceNode*	in = new InstanceNode;
	in->addOwnSlots(ossn);
	doc.addInstances(in);
	doc.finalizeMain("/",XML_SCHEMA_NAME);
	doc.printXMLTree();
	fclose(index);
}

