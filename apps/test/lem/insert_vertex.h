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
****************************************************************************/

#ifndef __VCG_TRI_INSERT_VERTEX
#define __VCG_TRI_INSERT_VERTEX


#include<vcg\simplex\face\pos.h>
#include<vcg\simplex\face\topology.h>
#include<vcg\complex\trimesh\update\topology.h>
#include<vcg\complex\trimesh\allocate.h>

/// This Class is used for insertiong a vertex in a face
namespace vcg{
  namespace tri{	
	
  template <class FACE_TYPE>
  void FFAttach(FACE_TYPE *f0,int z0,FACE_TYPE *f1,int z1)
  {
	f0->FFp(z0)=f1;
	f0->FFi(z0)=z1;
	f1->FFp(z1)=f0;
	f1->FFi(z1)=z0;
  }


   template <class MESH_TYPE> 
  ///insert a vertex iside a face and re-triangolarize v will be pointer to inserted vertex
  void InsertVertEdge(MESH_TYPE &m,typename MESH_TYPE::FaceType** f,int edge,typename MESH_TYPE::VertexType *&v)
  {
	std::vector<MESH_TYPE::FaceType **> local_var;
    local_var.push_back(f);
	MESH_TYPE::VertexIterator Vi=vcg::tri::Allocator<MESH_TYPE>::AddVertices(m,1);
	MESH_TYPE::FaceIterator Finit=vcg::tri::Allocator<MESH_TYPE>::AddFaces(m,4,local_var);
	
	if (MESH_TYPE::HasVFTopology())
			Vi->VFp()=0;
		
	MESH_TYPE::FaceIterator Fi=Finit;
	

	//take new faces added
	MESH_TYPE::FaceType *Fr0=&(*Fi);
	Fi++;
	MESH_TYPE::FaceType *Fr1=&(*Fi);
	Fi++;
	MESH_TYPE::FaceType *Fl0=&(*Fi);
	Fi++;
	MESH_TYPE::FaceType *Fl1=&(*Fi);

	//old faces
	MESH_TYPE::FaceType *fdl=(*f);
	MESH_TYPE::FaceType *fdr=(*f)->FFp(edge);

	int edgel=edge;
	int edger=(*f)->FFi(edge);

	//opposite vertex used to build new triangles
	MESH_TYPE::VertexType *voppl=fdl->V((edgel+2)%3);
	MESH_TYPE::VertexType *voppr=fdr->V((edger+2)%3);
	
	//inizialize the for new faces
	Fl0->V(0)=voppl;
	Fl0->V(1)=fdl->V(edgel);
	Fl0->V(2)=&(*Vi);
	Fl1->V(0)=voppl;
	Fl1->V(1)=&(*Vi);
	Fl1->V(2)=fdl->V((edgel+1)%3);

	Fr0->V(0)=voppr;
	Fr0->V(1)=&(*Vi);
	Fr0->V(2)=fdr->V((edger+1)%3);
	Fr1->V(0)=voppr;
	Fr1->V(1)=fdr->V(edger);
	Fr1->V(2)=&(*Vi);

	//VFTopology setting
	if (MESH_TYPE::HasVFTopology())
	{	
		for (int i=0;i<3;i++)
		{
			//initial settings
			Fl0->VFp(i)=0;
			Fl1->VFp(i)=0;
			Fr0->VFp(i)=0;
			Fr1->VFp(i)=0;

			//append new faces to VF topology
			vcg::face::VFAppend<MESH_TYPE::FaceType>(Fl0,i);
			vcg::face::VFAppend<MESH_TYPE::FaceType>(Fl1,i);
			vcg::face::VFAppend<MESH_TYPE::FaceType>(Fr0,i);
			vcg::face::VFAppend<MESH_TYPE::FaceType>(Fr1,i);

			//erase old faces from VF topology
			vcg::face::VFDetach<MESH_TYPE::FaceType>((*fdl),i);		
			vcg::face::VFDetach<MESH_TYPE::FaceType>((*fdr),i);		
		
		}
	}

	//FFTopology setting
	if (MESH_TYPE::HasFFTopology())
		{
			//attach old faces that was attached by 
			//the one that should be substituted
			FFAttach(Fl0,0,fdl->FFp((edgel+2)%3),fdl->FFi((edgel+2)%3));
			FFAttach(Fl1,2,fdl->FFp((edgel+1)%3),fdl->FFi((edgel+1)%3));
		
			FFAttach(Fr0,2,fdr->FFp((edger+1)%3),fdr->FFi((edger+1)%3));
			FFAttach(Fr1,0,fdr->FFp((edger+2)%3),fdr->FFi((edger+2)%3));

			//then connect between thenselfes
			
			FFAttach(Fl0,2,Fl1,0);
			FFAttach(Fl0,1,Fr0,1);

			FFAttach(Fl1,1,Fr1,1);
			FFAttach(Fr0,0,Fr1,2);
			


		}

	//finally set as deleted the old faces
	fdl->SetD();
	fdr->SetD();
	v=&(*Vi);
  }


  template <class MESH_TYPE> 
  ///insert a vertex iside a face and re-triangolarize v will be pointer to inserted vertex
  void InsertVertFace(MESH_TYPE &m,typename MESH_TYPE::FaceType** f,typename MESH_TYPE::VertexType *&v)
  {
    assert(!(*f)->IsD());
	std::vector<MESH_TYPE::FaceType **> local_var;
    local_var.push_back(f);
	MESH_TYPE::VertexIterator Vi=vcg::tri::Allocator<MESH_TYPE>::AddVertices(m,1);
	MESH_TYPE::FaceIterator Finit=vcg::tri::Allocator<MESH_TYPE>::AddFaces(m,3,local_var);
	
	if (MESH_TYPE::HasVFTopology())
			Vi->VFp()=0;
		
	MESH_TYPE::FaceIterator Fi=Finit;
	MESH_TYPE::FaceType *F;
	MESH_TYPE::FaceType *fd=(*f);
	//set vertex pointer of new face
	for (int i=0;i<3;i++)
	{
		F=&(*Fi);

		assert(!fd->V(i)->IsD());
		assert(!Vi->IsD());

		F->V(0)=fd->V(i);
		F->V(1)=fd->V((i+1)%3);
		F->V(2)=&(*Vi);

		

		//assign topology in substitution of the old one
		if (MESH_TYPE::HasFFTopology())
		{
			FFAttach(F,0,fd->FFp(i),fd->FFi(i));
		}

		if (MESH_TYPE::HasVFTopology())
		{
			F->VFp(0)=0;
			F->VFp(1)=0;
			F->VFp(2)=0;
			//put new faces on list of the old vertex and new one
			vcg::face::VFAppend<MESH_TYPE::FaceType>(F,0);
			vcg::face::VFAppend<MESH_TYPE::FaceType>(F,1);
			vcg::face::VFAppend<MESH_TYPE::FaceType>(F,2);

			vcg::face::VFDetach<MESH_TYPE::FaceType>((*fd),i);
		}
		Fi++;
	}
	//then attach the faces between themselfes
	Fi=Finit;	
	MESH_TYPE::FaceIterator Fsucc=Fi;
	Fsucc++;

	MESH_TYPE::FaceType *F0=&(*Fi);
	MESH_TYPE::FaceType *F1=&(*Fsucc);

	FFAttach<MESH_TYPE::FaceType>(F0,1,F1,2);

	Fi++;
	Fsucc++;
	F0=&(*Fi);
	F1=&(*Fsucc);

	FFAttach<MESH_TYPE::FaceType>(F0,1,F1,2);

	Fsucc=Finit;
	Fi++;
	F0=&(*Fi);
	F1=&(*Fsucc);

	FFAttach<MESH_TYPE::FaceType>(F0,1,F1,2);

	//at the end set as deleted the old face that was substituted
	fd->SetD();
    v=&(*Vi);

  }


}
}
#endif 