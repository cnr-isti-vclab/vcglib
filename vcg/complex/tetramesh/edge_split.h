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
#ifndef __VCG_TETRA_EDGE_SPLIT
#define __VCG_TETRA_EDGE_SPLIT

#include <vcg/simplex/tetrahedron/pos.h>
#include <vcg/complex/tetramesh/allocate.h>
#include <vcg/complex/tetramesh/update/topology.h>
#include <vcg/space/tetra.h>
namespace vcg{
namespace tetra{

/** \addtogroup tetramesh */
/*@{*/
/// This Class is used for split the edges

template <class TETRA_MESH_TYPE> 
class EdgeSplit
	{
	public:
  /// The tetrahedral mesh type
  typedef	typename TETRA_MESH_TYPE TetraMeshType;
  /// The tetrahedron type
  typedef	typename TetraMeshType::TetraType TetraType;
	/// The vertex type
	typedef	typename TetraType::VertexType VertexType;
  /// The vertex iterator type
  typedef	typename TetraMeshType::VertexIterator VertexIterator;
  /// The tetra iterator type
  typedef	typename TetraMeshType::TetraIterator TetraIterator;
  /// The coordinate type
	typedef	typename TetraType::VertexType::CoordType CoordType;
  ///the container of tetrahedron type
  typedef typename TetraMeshType::TetraContainer TetraContainer;
  ///the container of vertex type
  typedef typename TetraMeshType::VertexContainer VertexContainer;
  /// The HEdgePos type
  typedef PosLoop<TetraType> PosType;
  /// The topology updater type
  typedef vcg::tetra::UpdateTetraTopology<VertexContainer,TetraContainer> Topology;
  /// The allocator type
  typedef  vcg::tetra::Allocator<TetraMeshType> TetraAllocator;
  /// Default Constructor
	EdgeSplit()
		{
		};

  ~EdgeSplit()
		{
		};

private:
///the tetrahedron that must mark as Deleted after substitution
TetraType* _toDel[30];
///The number of tetrahedrons that are substituted
int _nT;

Topology _Topo;

///add a vertex into the edge at distance alfa ( 0<alfa<1) to the first vertex of the edge

VertexType* _AddVertexEdge(TetraMeshType &tm,const TetraType &t,const int &edge,const double &alfa)
{
  VertexType *v0=(VertexType*)t.V(Tetra::VofE(edge,0));
	VertexType *v1=(VertexType*)t.V(Tetra::VofE(edge,1));
	Allocator<TetraMeshType> All= Allocator<TetraMeshType>();
	VertexIterator vn=All.AddVertices(tm,1);
	vn->Flags()=0;
	vn->VTb()=NULL;
	vn->VTi()=-1;
  vn->P()=(v0->P()*alfa)+(v1->P()*(1.0f-alfa));
  return (&(*vn));
}

///set the default v-t topology
void _SetDefultVTTopology(TetraType *t)
{
unsigned int j;
for (j=0;j<4;j++)
{
	t->TVp(j) = NULL;
	t->TVi(j) = -1;
}
}

/// Transform the vertex index according to rotation of the tetraedron 
/// that trasform it in the basic case
static int _GetMapVertEdgeRot(const int &indexE,const int &indexV)
	{	
		static int mapvertedgerot[12][4]={
						  {0,3,1,2},
						  {0,1,2,3},
						  {0,2,3,1},
						  {1,3,2,0},
						  {1,0,3,2},
						  {2,1,3,0},

						  {1,2,0,3},
						  {2,3,0,1},
						  {3,1,0,2},
						  {2,0,1,3},
						  {3,2,1,0},
						  {3,0,2,1},
						};
		assert ((indexE<12)&&(indexV<4));
		return mapvertedgerot[indexE][indexV];
	}	

/// Transform the face index according to rotation of the tetraedron 
/// that trasform it in the basic case
  static int _GetMapFaceEdgeRot(const int &indexE,const int &indexF)
	{	
		static int mapfaceedgerot[12][4]={
						  {1,2,0,3},
						  {0,1,2,3},
						  {2,0,1,3},
						  {3,1,0,2},
						  {1,0,3,2},
						  {3,0,2,1},

						  {0,3,1,2},
						  {2,3,0,1},
						  {1,3,2,0},
						  {0,2,3,1},
						  {3,2,1,0},
						  {2,1,3,0},
						};
		assert ((indexE<12)&&(indexF<4));
		return mapfaceedgerot[indexE][indexF];
	}	

/// Returns the rotation sense during the loop on the edge to divide
/// according to rotation of the tetraedron that trasform it in the basic case
 static int _GetDirRot(int indexE,int indexF)
	{	
    static int mapfaceedgerot[12][4]={
						  {2,0,-1,-1},
						  {0,-1,2,-1},
						  {-1,2,0,-1},
						  {2,-1,-1,0},
						  {-1,0,-1,2},
						  {-1,-1,2,0},
               
              {0,2,-1,-1},
						  {2,-1,0,-1},
						  {-1,0,2,-1},
						  {0,-1,-1,2},
						  {-1,2,-1,0},
						  {-1,-1,0,2},
						};

		assert ((indexE<12)&&(indexF<4));
		return mapfaceedgerot[indexE][indexF];
	}	

///Built an Half edge on tetrahedron t using edge edge
PosType _FindPos(TetraType *t,int edge)
{	
	int face0=Tetra::FofE(edge,0);
	int ve0=Tetra::VofE(edge,0);
	PosType pos(t,face0,edge,ve0);
	return pos;
}

///Assert the right order of vertex that compose the tetrahedron
void _AssertOrder(TetraType *t,VertexType *v0,VertexType *v1,VertexType *v2,VertexType *v3)
{
  assert(t->V(0)==v0);
	assert(t->V(1)==v1);
	assert(t->V(2)==v2);
	assert(t->V(3)==v3);
}

///Connect trought Tetrahedron-Tetrahedron Topology t0 and t1 with faces i0 and i1
void _ConnectTTTopology(TetraType *t0,int i0,TetraType *t1,int i1) 
{
  assert((i0>=0)&&(i0<4));
  assert((i1>=0)&&(i1<4));
  assert((!t0->IsD())&&(!t1->IsD()));
  t0->TTp(i0)=t1;
  t0->TTi(i0)=i1;
  t1->TTp(i1)=t0;
  t1->TTi(i1)=i0;
  assert( (((t0->TTp(i0))->TTp(t0->TTi(i0)))==t0));
  assert( (((t1->TTp(i1))->TTp(t1->TTi(i1)))==t1));
}

///Divide the tetrahadron in pos in two tetrahedrons using the new vertex vnew
void _Divide(PosType pos,VertexType *vnew,TetraType *newtp0,TetraType *newtp1,bool invert)
{
      int curredge=pos.E();
     
			//control if the edge vertices arein the right order for the table
			if (invert)
			curredge+=6;

			//find the new position to vertex according
			int ie0=_GetMapVertEdgeRot(curredge,0);
			int ie1=_GetMapVertEdgeRot(curredge,2);
			int in0=_GetMapVertEdgeRot(curredge,1);
			int in1=_GetMapVertEdgeRot(curredge,3);
		
			//as first the ones that appartain to the selected curredge
			VertexType *ve0=pos.T()->V(ie0);
			VertexType *ve1=pos.T()->V(ie1);

			//and after the others that will be in the cutting plane
			VertexType *vn0=pos.T()->V(in0);
			VertexType *vn1=pos.T()->V(in1);
		
			newtp0->V(0)=ve0;
      newtp0->V(1)=vn0;
      newtp0->V(2)=vnew;
      newtp0->V(3)=vn1;

      newtp1->V(0)=vnew;
      newtp1->V(1)=vn0;
      newtp1->V(2)=ve1;
      newtp1->V(3)=vn1;
     
		  //right order of the vertices
#ifdef _DEBUG
      _AssertOrder(newtp0,ve0,vn0,vnew,vn1);
			_AssertOrder(newtp1,vnew,vn0,ve1,vn1);
			//end asserts
#endif
}

bool _InvertRotation(PosType pos)
{
  return (pos.V()!=Tetra::VofE(pos.E(),0));
}

///substitute the told tetrahedon on VT topology with newtp0 and newtp1 as created
void _SubstituteVTTopology(TetraType *told,TetraType *newtp0,TetraType *newtp1)
{
  _SetDefultVTTopology(newtp0);
	_SetDefultVTTopology(newtp1);
			
  //detach the old tetrahedron from VTtopology
  _Topo.DetachVTTopology(told);
	//tetrahedron 0
	_Topo.InsertVTTopology(newtp0);
	//tetrahedron 1
	_Topo.InsertVTTopology(newtp1);		
}

///control if the connections between tetrahedron created have the right shared vertices
void _ControlConnection(TetraType *oldtp0,TetraType *newtp0)
{
  VertexType *v00=oldtp0->V(0);
  VertexType *v01=oldtp0->V(1);
  VertexType *v02=oldtp0->V(2);
  VertexType *v03=oldtp0->V(3);

  VertexType *v10=newtp0->V(0);
  VertexType *v11=newtp0->V(1);
  VertexType *v12=newtp0->V(2);
  VertexType *v13=newtp0->V(3);
  
  assert(((v00==v10)&&(v02==v12))||((v00==v12)&&(v02==v10)));
  assert(((v01==v13)&&(v03!=v11))||((v01!=v13)&&(v03==v11)));
}

///set as extern the 4 faces of the tetrahedron
void _SetDefaultTTExtern(TetraType *t)
{
  for (int y=0;y<4;y++)
			{
				t->TTp(y)=t;
				t->TTi(y)=y;
      }
}

///substitute in Tetra Tetra Topology the tetrahedron old_t with new_t in according 
///to face and edge 

void _SubstituteTTTopology(TetraType *old_t,TetraType *new_t,int edgerot,int face)
{
   int indexface=_GetMapFaceEdgeRot(edgerot,face);

			if (old_t->IsBorderF(indexface))
			{
	        new_t->TTp(face)=new_t;
			  	new_t->TTi(face)=face;
      }
      else
      {
       	  TetraType *tetrad=old_t->TTp(indexface);
			    int fad=old_t->TTi(indexface);
          _ConnectTTTopology(new_t,face,tetrad,fad);
          assert (!tetrad->IsD());
      }
}

/// sobstitute the old tetrahedrons that share the edge in pos with new ones
/// that share the vertex vnew that divide the old edge

void _AddNewTetrahedrons(TetraMeshType &tm,PosType pos,VertexType *vnew)
{

	TetraType *oldtp0=NULL;
	TetraType *oldtp1=NULL;

	TetraType *newtp0;
  TetraType *newtp1;

  TetraType *firsttp0=NULL;
  TetraType *firsttp1=NULL;


  int curredge;
  int direction=-1;
  bool invert=false;

TetraAllocator All=TetraAllocator();
pos.Reset();
_nT=0;

while (!pos.LoopEnd())
  {	
			assert(!pos.T()->IsD());
      
      invert=_InvertRotation(pos);

			//CREATE THE NEW TETRAHEDRONS

			//create the new ones putting the veritices in the right order
      TetraIterator ti=All.AddTetra(tm,2);		
      newtp0 = &(*ti);
			ti++;
      newtp1 = &(*ti);	
      

      _Divide(pos,vnew,newtp0,newtp1,invert);

#ifdef _DEBUG
      if ((oldtp0!=NULL)&&(!pos.Jump()))
      _ControlConnection(oldtp0,newtp0);
      if ((oldtp1!=NULL)&&(!pos.Jump()))
      _ControlConnection(oldtp1,newtp1);
#endif
      
      // SUBSTITUTE NEW TETRAHEDRONS ON VT TOPOLOGY
      if (tm.HasVTTopology())
      _SubstituteVTTopology(pos.T(),newtp0,newtp1);
			

			//THEN SET THE T-T TOPOLOGY

      _SetDefaultTTExtern(newtp0);
      _SetDefaultTTExtern(newtp1);
			
      curredge=pos.E();
      if (invert)
			  curredge+=6;

      //face3 
      _SubstituteTTTopology(pos.T(),newtp1,curredge,3);
		
			//face1 
      _SubstituteTTTopology(pos.T(),newtp0,curredge,1);
      
			//now I set t-t topology between themselfes

			_ConnectTTTopology(newtp0,3,newtp1,1);

      if (pos.Jump())
      {
        vnew->SetB();
        oldtp0=NULL;
        oldtp1=NULL;
      }
      
      direction=_GetDirRot(curredge,pos.F());
      assert(direction!=-1);
       //control the direction of moving
      if ((oldtp0!=NULL)&&(oldtp1!=NULL))
      {
        //direction=_GetDirRot(oldtp0,newtp0);

        //find direction of moving
        if (direction==0)
		  	{
         _ConnectTTTopology(oldtp0,0,newtp0,2);
         _ConnectTTTopology(oldtp1,0,newtp1,2);
		  	}
			  else 
        if (direction==2)
			  {
          _ConnectTTTopology(oldtp0,2,newtp0,0);
          _ConnectTTTopology(oldtp1,2,newtp1,0);
			  }
      }
      //assign if it is the first one
      if (firsttp0==NULL)
        firsttp0=newtp0;
      if (firsttp1==NULL)
        firsttp1=newtp1;

      oldtp0=newtp0;
      oldtp1=newtp1;

      _toDel[_nT]=pos.T();
      _nT++;
      pos.NextT();
      }
  
  //at the end I finish the connections
  if (!(pos.Jump())&&(direction==0)&&(firsttp0!=NULL)&&(firsttp1!=NULL)&&(oldtp0!=NULL)&&(oldtp1!=NULL))
   {
     _ConnectTTTopology(oldtp0,0,firsttp0,2);
     _ConnectTTTopology(oldtp1,0,firsttp1,2);
   }
  else if (!(pos.Jump())&&(direction==2)&&(firsttp0!=NULL)&&(firsttp1!=NULL)&&(oldtp0!=NULL)&&(oldtp1!=NULL))
   {
      _ConnectTTTopology(oldtp0,2,firsttp0,0);
      _ConnectTTTopology(oldtp1,2,firsttp1,0);
   }
  else if (pos.Jump())
    vnew->SetB();

}

///Mark as deleted the tetrahedron that must be substituted
void _DeleteOldTetra()
{
  for (int i=0;i<_nT;i++)
  _toDel[i]->SetD();
}

//=========================================================================

public:

/// Split the edge with local remeshing 
/// Tetrahedron-Tetrahedron topology is required
VertexType* DoSplit(TetraMeshType &tm,TetraType *t,int edge,double alfa)
{
	assert(!t->IsD());
  assert(tm.HasTTTopology());
	assert((alfa>0)&&(alfa<1));
	assert((edge>=0)&&(edge<6));
	VertexType *vnew=_AddVertexEdge(tm,*t,edge,alfa);
	_AddNewTetrahedrons(tm,_FindPos(t,edge),vnew);
  _DeleteOldTetra();
	return(vnew);
}

};//end class

}//end namespace tetra
}//end namespace vcg
#endif