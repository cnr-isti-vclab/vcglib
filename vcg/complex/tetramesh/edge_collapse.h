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

#ifndef __VCG_TETRA_EDGE_COLLAPSE
#define __VCG_TETRA_EDGE_COLLAPSE

#include <vcg/space/tetra3.h>
#include <vcg/complex/tetramesh/update/topology.h>
#include <vcg/complex/tetramesh/update/normal.h>

//#include <vcg/complex/trimesh/edge_collapse.h>
//
////used to verify the edge collapse on the surface of the mesh
//#include <vcg/simplex/face/with/afav.h>
//#include <vcg/complex/trimesh/base.h>
//#include <vcg/complex/trimesh/update/topology.h>

#include <vcg/complex/tetramesh/update/triconvert.h>


namespace vcg{
namespace tetra{

/** \addtogroup tetramesh */
/*@{*/
/// This Class is used for the edge collapse

template <class TETRA_MESH_TYPE> 
class EdgeCollapse
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
  /// The scalar type
  typedef	typename TetraMeshType::VertexType::ScalarType ScalarType;
  ///the container of tetrahedron type
  typedef typename TetraMeshType::TetraContainer TetraContainer;
  ///the container of vertex type
  typedef typename TetraMeshType::VertexContainer VertexContainer;
  /// The HEdgePos type
  typedef Pos<TetraType> PosType;
  /// The HEdgePos Loop type
  typedef PosLoop<TetraType> PosLType;
  /// The topology updater type
  typedef vcg::tetra::UpdateTetraTopology<VertexContainer,TetraContainer> Topology;
  ///the normal updater type
  typedef  vcg::tetra::UpdateNormals<TetraMeshType> UpdateNormals;


  /// Default Constructor
	EdgeCollapse()
		{
      c0=0;
      flip=0;
      lkv=0;
      lke=0;

		};

  ~EdgeCollapse()
		{
		};

  private:
  vector<TetraType*> To_Del;
  Topology _Topo;
  UpdateNormals _UN;
  typedef pair <VertexType*,VertexType*> VertPair;
  typedef pair <int,int> FacePair;

public:
  int c0;
  int flip;
  int lkv;
  int lke;

private:
///select the 2 faces that does not share the edge
FacePair _FindNoEdgeFace(TetraType *t,int edge)
{	
				//as first I find the 2 faces on the opposite sides of the egde
        int fa0=Tetra::FofE(edge,0);
				int fa1=Tetra::FofE(edge,1);

				//then find the faces that remain
				int fa2=(fa0+1)%4;
				while ((fa2==fa0)||(fa2==fa1))
				{
					fa2=(fa2+1)%4;
				}
				int fa3=(fa2+1)%4;
				while ((fa3==fa0)||(fa3==fa1)||(fa3==fa2))
				{
					fa3=(fa3+1)%4;
				}
        return FacePair(fa2,fa3);
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

///Calculate the volume on the vertex resulting after collapse...

ScalarType _VolumeSimulateCollapse(PosType Pos,ScalarType alfa)
{
  VertexType *Vrem=(Pos.T()->V(Tetra::VofE(Pos.E(),0)));
  VertexType *Vdel=(Pos.T()->V(Tetra::VofE(Pos.E(),1)));

  if (Vrem!=Pos.T()->V(Pos.V()))
    swap<VertexType*>(Vdel,Vrem);

  ScalarType vol=0;
  CoordType oldpos = Vrem->P();
  CoordType newpos=((Vrem->P()*alfa)+(Vdel->P()*(1.f-alfa)));

//move vertex that remain in the new position
  Vrem->P() = newpos;

  vector< TetraType *>::iterator ti=_Sets.no_E.begin(); 

  while (ti!=_Sets.no_E.end())
  {
        Tetra3<ScalarType> T=Tetra3<ScalarType>();
        T.P0(0)=(*ti)->V(0)->cP();
        T.P1(0)=(*ti)->V(1)->cP();
        T.P2(0)=(*ti)->V(2)->cP();
        T.P3(0)=(*ti)->V(3)->cP();
        vol+=T.ComputeVolume();
        ti++;
  }
  Vrem->P()=oldpos;
  return vol;
}

///return the sum of volumes of the union of  stars on vertices
ScalarType _VolumeUnion()
{
  vector< TetraType *>::iterator ti=_Sets.v0_U_v1.begin(); 
  ScalarType vol=0;
  while (ti!=_Sets.v0_U_v1.end())
  {
    vol+=(*ti)->Volume();
    ti++;
  }
  return vol;
}

#ifdef _DEBUG
void _AssertingVolume(TetraType *t)
{
		assert(t->ComputeVolume() >0);
}
#endif

///collpse de edge specified by pos (the first vertex on edge remain)

void _Collapse(PosType p,ScalarType alfa)
{
      VertexType *Vrem=(p.T()->V(Tetra::VofE(p.E(),0)));
      VertexType *Vdel=(p.T()->V(Tetra::VofE(p.E(),1)));
      Vrem->P()=(Vrem->P()*alfa)+(Vdel->P()*(1.f-alfa));
      PosLType pos(p.T(),p.F(),p.E(),p.V());
      pos.Reset();
      To_Del.reserve(40);
      To_Del.clear();
			while (!pos.LoopEnd())
			{
				//get the two faces that doesn't share the edge
        FacePair fp=_FindNoEdgeFace(pos.T(),pos.E());
        int fa0=fp.first;
        int fa1=fp.second;

				//now set the T-T topology on that faces
				TetraType *tleft=pos.T()->TTp(fa0);
				TetraType *tright=pos.T()->TTp(fa1);
				int ileft=pos.T()->TTi(fa0);
				int iright=pos.T()->TTi(fa1);
				
        //in this case I cannot do the collapse
				assert (!((pos.T()==tleft)&&(pos.T()==tright)));
				
        //case no one is extern face
				if ((!pos.T()->IsBorderF(fa0))&&(!pos.T()->IsBorderF(fa1)))
          //connect the 2 tetrahedrons
          _ConnectTTTopology(tleft,ileft,tright,iright);	
				else
				  //case f2 is an extern face
				if (pos.T()->IsBorderF(fa0))
				{
					tright->TTp(iright)=tright;
					tright->TTi(iright)=iright;
				}
        
        else //case fa1 is an extern face
        //if ((pos.T()->IsBorderF(fa3))
				{
					tleft->TTp(ileft)=tleft;
					tleft->TTi(ileft)=ileft;	
				}
				
				//end setting T-T topology

				//setting the V-T topology

				//i remove the tetrahedrons that have the edge
				// to collapse
				_Topo.DetachVTTopology(pos.T());
				//end setting the V-T topology	
        To_Del.push_back(pos.T());
				pos.NextT();
			  
        tm.tn--;
      }

      //delting old tetrahedrons
      vector<TetraType*>::iterator ti;
      for (ti=To_Del.begin();ti<To_Del.end();ti++)
        (*ti)->SetD();

			//now I cycle on the tetrahedron that had the old vertex
			//reassegning the new one.

      VTIterator< TetraType> VTi(Vdel->VTb(),Vdel->VTi());
			while (!VTi.End())
			{	
        TetraType *T_Change=VTi.Vt();
        int index=VTi.Vi();
        //VTi++;
				//assegning the vertex that remain
				T_Change->V(index)=Vrem;
        _Topo.DetachVTTopology(Vdel,T_Change);
				_Topo.InsertVTTopology(Vrem,index,T_Change);
				//that's cause i restart everytime in the chain 
				//from the vertex
        VTi.Vt()=Vdel->VTb();
				VTi.Vi()=Vdel->VTi();
#ifdef _DEBUG
       _AssertingVolume(T_Change);
#endif
       
			}	
      if (Vdel->IsB())
        Vrem->SetB();
      //set as deleted the vertex
			Vdel->SetD();
			tm.vn--;

}

struct Face
{
			VertexType* v[3];

      Face(	VertexType* a, VertexType* b,VertexType* c)
           {
						assert((a!=b)&&(b!=c)&&(a!=c));
            v[0]=a;
            v[1]=b;
            v[2]=c;
            sort(v,v+3);
						}

      const bool operator <(const Face & f) const 
      {
				return ((v[0]==f.v[0])?((v[1]==f.v[1])?(v[2]<f.v[2]):(v[1]<f.v[1])):(v[0]<f.v[0]));
			}

			const bool operator ==(const Face & f) const {
			return ((v[0]==f.v[0])&&(v[1]==f.v[1])&&(v[2]==f.v[2]));
			}

      };

struct Edge{
			VertexType* v0;
      VertexType* v1;
      Edge(	VertexType* a, VertexType* b){
						assert(a!=b);
						if(a<b) 
							{v0=a;v1=b;}
						else
						{v1=a;v0=b;}
      }
			
			const bool operator <(const Edge & e) const {
				return (v0==e.v0)?(v1<e.v1):(v0<e.v0);
			}

			const bool operator ==(const Edge & e) const {
			return (v0==e.v0)&&(v1==e.v1);
			}

};

map<Edge,char> EdgeMark;
map<Face,char> FaceMark;

VertexType _dummyV;

void orMarkE(Edge E,char M)
{
 map<Edge,char>::iterator EI;
 EI=EdgeMark.find(E);
 if (EI==EdgeMark.end())
  EdgeMark.insert (pair<Edge,char>(E,M));
 else
   (*EI).second|=M;
}

bool isMarkedE(Edge E,char M)
{
 map<Edge,char>::iterator EI;
 EI=EdgeMark.find(E);
 if (EI==EdgeMark.end())
    return false;
 else return ((*EI).second & M);
}

void orMarkF(Face F,char M)
{
 map<Face,char>::iterator FI;
 FI=FaceMark.find(F);
 if (FI==FaceMark.end())
  FaceMark.insert (pair<Face,char>(F,M));
 else
   (*FI).second|=M;
}

bool isMarkedF(Face F,char M)
{
 map<Face,char>::iterator FI;
 FI=FaceMark.find(F);
 if (FI==FaceMark.end())
    return false;
 else return ((*FI).second & M);
}

///this structure is used to find the sets that are used in link conditions and more
struct TetraSets
{
  std::vector <TetraType*> v0;
  std::vector <TetraType*> v1;
  std::vector <TetraType*> v0_U_v1;
  std::vector <TetraType*> no_E;
  std::vector <TetraType*> E;
  std::vector <char> indexE;
  std::vector <char> indexv0;
  std::vector <char> indexv1;
};

TetraSets _Sets;

/////verify if the collapse can done looking to the edges 
//bool _LinkConditionsE(PosType pos) 
//{ 	
//    const int LINK_V0 = 0x00000001;
//    const int LINK_EE = 0x00000002;
//   
//    EdgeMark.clear();
//
//		// Mark edges of ve0 
//    vector< TetraType *>::iterator ti=_Sets.v0.begin(); 
//    vector< char >::iterator en=_Sets.indexv0.begin();
//    while (ti!=_Sets.v0.end())
//    {
//      for(int i=0;i<6;i++)
//        //put the edge of each tetrahedron on the map
//        orMarkE(Edge((*ti)->V(Tetra::VofE(i,0)),(*ti)->V(Tetra::VofE(i,1))),LINK_V0);
//
//     
//      //put dummy edge
//      for (int f=0;f<3;f++)
//      {
//       int f_test=Tetra::FofV((*en),f);
//        if  ((*ti)->IsBorderF(f_test))
//        {
//          orMarkE(Edge((*ti)->V(Tetra::VofF(f_test,0)),&_dummyV),LINK_V0);
//          orMarkE(Edge((*ti)->V(Tetra::VofF(f_test,1)),&_dummyV),LINK_V0);
//          orMarkE(Edge((*ti)->V(Tetra::VofF(f_test,2)),&_dummyV),LINK_V0);
//        }
//      }
//     ti++;
//     en++;
//    }
//    
//    ti=_Sets.E.begin(); 
//    en=_Sets.indexE.begin();
//    //mark them as intersection
//    while (ti!=_Sets.E.end())
//    {
//      for(int i=0;i<6;i++)
//        //put the edge of each tetrahedron on the map
//        orMarkE(Edge((*ti)->V(Tetra::VofE(i,0)),(*ti)->V(Tetra::VofE(i,1))),LINK_EE);
//      //  //edges with dummy vertex
//
//      //faces on the edge
//      int f0=Tetra::FofE((*en),0);
//      int f1=Tetra::FofE((*en),1);
//
//      if  ((*ti)->IsBorderF(f0))
//      {
//          orMarkE(Edge((*ti)->V(Tetra::VofF(f0,0)),&_dummyV),LINK_EE);
//          orMarkE(Edge((*ti)->V(Tetra::VofF(f0,1)),&_dummyV),LINK_EE);
//          orMarkE(Edge((*ti)->V(Tetra::VofF(f0,2)),&_dummyV),LINK_EE);
//      }
//
//      if  ((*ti)->IsBorderF(f1))
//      {
//          orMarkE(Edge((*ti)->V(Tetra::VofF(f1,0)),&_dummyV),LINK_EE);
//          orMarkE(Edge((*ti)->V(Tetra::VofF(f1,1)),&_dummyV),LINK_EE);
//          orMarkE(Edge((*ti)->V(Tetra::VofF(f1,2)),&_dummyV),LINK_EE);
//      }
//
//      ti++;
//      en++;
//    }
//   
//    //and at the end I verify if the intersection is equal to the star of the edge
//    ti=_Sets.v1.begin(); 
//    en=_Sets.indexv1.begin();
//    while (ti!=_Sets.v1.end())
//    {
//      for(int i=0;i<6;i++)
//      {
//        Edge e_test=Edge((*ti)->V(Tetra::VofE(i,0)),(*ti)->V(Tetra::VofE(i,1)));
//        if ((isMarkedE(e_test,LINK_V0))&&(!isMarkedE(e_test,LINK_EE)))
//        {        
//          lke++;
//          return false;
//        }
//      }
//
//      //dummy edges control
//      //put dummy edge
//      for (int f=0;f<3;f++)
//      {
//       int f_test=Tetra::FofV((*en),f);
//        if  ((*ti)->IsBorderF(f_test))
//        {
//          //control all the 3 edges
//          Edge e_test0=Edge((*ti)->V(Tetra::VofF(f_test,0)),&_dummyV);
//          Edge e_test1=Edge((*ti)->V(Tetra::VofF(f_test,1)),&_dummyV);
//          Edge e_test2=Edge((*ti)->V(Tetra::VofF(f_test,2)),&_dummyV);
//          if (((isMarkedE(e_test0,LINK_V0))&&(!isMarkedE(e_test0,LINK_EE)))||
//            ((isMarkedE(e_test1,LINK_V0))&&(!isMarkedE(e_test1,LINK_EE)))||
//            ((isMarkedE(e_test2,LINK_V0))&&(!isMarkedE(e_test2,LINK_EE))))
//          {        
//            lke++;
//            return false;
//          }
//        }
//      }
//      ti++;
//      en++;
//    }
//    return true;
//}


///verify if the collapse can done looking to the edges 
bool _LinkConditionsF(PosType pos) 
{ 	
    const int LINK_V0 = 0x00000001;
    const int LINK_EE = 0x00000002;
   
    EdgeMark.clear();

		// Mark edges of ve0 
    vector< TetraType *>::iterator ti=_Sets.v0.begin(); 
    vector< char >::iterator en=_Sets.indexv0.begin();
    while (ti!=_Sets.v0.end())
    {
      //put dummy face
      for (int f=0;f<3;f++)
      {
       int f_test=Tetra::FofV((*en),f);
        if  ((*ti)->IsBorderF(f_test))
        {
          orMarkF(Face((*ti)->V(Tetra::VofF(f_test,0)),(*ti)->V(Tetra::VofF(f_test,1)),&_dummyV),LINK_V0);
          orMarkF(Face((*ti)->V(Tetra::VofF(f_test,1)),(*ti)->V(Tetra::VofF(f_test,2)),&_dummyV),LINK_V0);
          orMarkF(Face((*ti)->V(Tetra::VofF(f_test,2)),(*ti)->V(Tetra::VofF(f_test,0)),&_dummyV),LINK_V0);
        }
      }
     ti++;
     en++;
    }
    
    ti=_Sets.E.begin(); 
    en=_Sets.indexE.begin();
    //mark them as intersection
    while (ti!=_Sets.E.end())
    {
      //faces on the edge
      int f0=Tetra::FofE((*en),0);
      int f1=Tetra::FofE((*en),1);

      if  ((*ti)->IsBorderF(f0))
      {
          orMarkF(Face((*ti)->V(Tetra::VofF(f0,0)),(*ti)->V(Tetra::VofF(f0,1)),&_dummyV),LINK_EE);
          orMarkF(Face((*ti)->V(Tetra::VofF(f0,1)),(*ti)->V(Tetra::VofF(f0,2)),&_dummyV),LINK_EE);
          orMarkF(Face((*ti)->V(Tetra::VofF(f0,2)),(*ti)->V(Tetra::VofF(f0,0)),&_dummyV),LINK_EE);
      }

      if  ((*ti)->IsBorderF(f1))
      {
          orMarkF(Face((*ti)->V(Tetra::VofF(f1,0)),(*ti)->V(Tetra::VofF(f1,1)),&_dummyV),LINK_EE);
          orMarkF(Face((*ti)->V(Tetra::VofF(f1,1)),(*ti)->V(Tetra::VofF(f1,2)),&_dummyV),LINK_EE);
          orMarkF(Face((*ti)->V(Tetra::VofF(f1,2)),(*ti)->V(Tetra::VofF(f1,0)),&_dummyV),LINK_EE);
      }

      ti++;
      en++;
    }
   
    //and at the end I verify if the intersection is equal to the star of the edge
    ti=_Sets.v1.begin(); 
    en=_Sets.indexv1.begin();
    while (ti!=_Sets.v1.end())
    {

      //dummy edges control
      for (int f=0;f<3;f++)
      {
       int f_test=Tetra::FofV((*en),f);
        if  ((*ti)->IsBorderF(f_test))
        {
          //control all the 3 edges
          Face f_test0=Face((*ti)->V(Tetra::VofF(f_test,0)),(*ti)->V(Tetra::VofF(f_test,1)),&_dummyV);
          Face f_test1=Face((*ti)->V(Tetra::VofF(f_test,1)),(*ti)->V(Tetra::VofF(f_test,2)),&_dummyV);
          Face f_test2=Face((*ti)->V(Tetra::VofF(f_test,2)),(*ti)->V(Tetra::VofF(f_test,0)),&_dummyV);
          if (((isMarkedF(f_test0,LINK_V0))&&(!isMarkedF(f_test0,LINK_EE)))||
            ((isMarkedF(f_test1,LINK_V0))&&(!isMarkedF(f_test1,LINK_EE)))||
            ((isMarkedF(f_test2,LINK_V0))&&(!isMarkedF(f_test2,LINK_EE))))
          {        
            lke++;
            return false;
          }
        }
      }
      ti++;
      en++;
    }
    return true;
}

///verify if the collapse can done looking to the edges 
bool _LinkConditionsE(PosType pos) 
{ 	
    const int LINK_V0 = 0x00000001;
    const int LINK_EE = 0x00000002;
   
    FaceMark.clear();

		// Mark edges of ve0 
    vector< TetraType *>::iterator ti=_Sets.v0.begin(); 
    vector< char >::iterator en=_Sets.indexv0.begin();
    while (ti!=_Sets.v0.end())
    {
      //put dummy edge
      for (int f=0;f<3;f++)
      {
       int f_test=Tetra::FofV((*en),f);
        if  ((*ti)->IsBorderF(f_test))
        {
          orMarkE(Edge((*ti)->V(Tetra::VofF(f_test,0)),&_dummyV),LINK_V0);
          orMarkE(Edge((*ti)->V(Tetra::VofF(f_test,1)),&_dummyV),LINK_V0);
          orMarkE(Edge((*ti)->V(Tetra::VofF(f_test,2)),&_dummyV),LINK_V0);
        }
      }
     ti++;
     en++;
    }
    
    ti=_Sets.E.begin(); 
    en=_Sets.indexE.begin();
    //mark them as intersection
    while (ti!=_Sets.E.end())
    {
      //faces on the edge
      int f0=Tetra::FofE((*en),0);
      int f1=Tetra::FofE((*en),1);

      if  ((*ti)->IsBorderF(f0))
      {
          orMarkE(Edge((*ti)->V(Tetra::VofF(f0,0)),&_dummyV),LINK_EE);
          orMarkE(Edge((*ti)->V(Tetra::VofF(f0,1)),&_dummyV),LINK_EE);
          orMarkE(Edge((*ti)->V(Tetra::VofF(f0,2)),&_dummyV),LINK_EE);
      }

      if  ((*ti)->IsBorderF(f1))
      {
          orMarkE(Edge((*ti)->V(Tetra::VofF(f1,0)),&_dummyV),LINK_EE);
          orMarkE(Edge((*ti)->V(Tetra::VofF(f1,1)),&_dummyV),LINK_EE);
          orMarkE(Edge((*ti)->V(Tetra::VofF(f1,2)),&_dummyV),LINK_EE);
      }

      ti++;
      en++;
    }
   
    //and at the end I verify if the intersection is equal to the star of the edge
    ti=_Sets.v1.begin(); 
    en=_Sets.indexv1.begin();
    while (ti!=_Sets.v1.end())
    {

      //dummy edges control
      for (int f=0;f<3;f++)
      {
       int f_test=Tetra::FofV((*en),f);
        if  ((*ti)->IsBorderF(f_test))
        {
          //control all the 3 edges
          Edge e_test0=Edge((*ti)->V(Tetra::VofF(f_test,0)),&_dummyV);
          Edge e_test1=Edge((*ti)->V(Tetra::VofF(f_test,1)),&_dummyV);
          Edge e_test2=Edge((*ti)->V(Tetra::VofF(f_test,2)),&_dummyV);
          if (((isMarkedE(e_test0,LINK_V0))&&(!isMarkedE(e_test0,LINK_EE)))||
            ((isMarkedE(e_test1,LINK_V0))&&(!isMarkedE(e_test1,LINK_EE)))||
            ((isMarkedE(e_test2,LINK_V0))&&(!isMarkedE(e_test2,LINK_EE))))
          {        
            lke++;
            return false;
          }
        }
      }
      ti++;
      en++;
    }
    return true;
}

///verify the link conditions for a collapse using vertices

bool _LinkConditionsV() 
{ 	
    const int LINK_V0 = VertexType::NewUserBit();
		const int LINK_V1 = VertexType::NewUserBit();
    const int LINK_EE = VertexType::NewUserBit();
    
		const int NOT_LINKED = ~(LINK_V0 | LINK_V1 | LINK_EE);
    _dummyV.Flags() &= NOT_LINKED;

    VertexType *vt0;
    VertexType *vt1;
    VertexType *vt2;
    VertexType *vt3;
		

   vector< TetraType *>::iterator ti=_Sets.v0_U_v1.begin(); 
    
    //reset all link flags
    while (ti!=_Sets.v0_U_v1.end())
    {
      for(int i=0;i<4;i++)
        (*ti)->V(i)->Flags() &= NOT_LINKED;
      ti++;
    }


    //also in the ones that appartain to the edge
    vector< char >::iterator en;
    ti=_Sets.E.begin(); 
    en=_Sets.indexE.begin();
    //reset all link flags for intersection and in the same
    //time mark them as intersection
    while (ti!=_Sets.E.end())
    {
      for(int i=0;i<4;i++)
      {
        (*ti)->V(i)->Flags() &= NOT_LINKED;
        (*ti)->V(i)->Flags() |= LINK_EE;
      }

      //dummy vertex

      //faces on the edge
      int f0=Tetra::FofE((*en),0);
      int f1=Tetra::FofE((*en),1);
      
      if  (((*ti)->IsBorderF(f0))||((*ti)->IsBorderF(f1)))
         _dummyV.Flags() |= LINK_EE;

      ti++;
      en++;
    }


		// Mark vertices of ve0 
    ti=_Sets.v0.begin();
    en=_Sets.indexv0.begin();

    while (ti!=_Sets.v0.end())
    {
      for(int i=0;i<4;i++)
        (*ti)->V(i)->Flags() |= LINK_V0;
      
      //dummy faces on the vertex
      int f0=Tetra::FofV((*en),0);
      int f1=Tetra::FofV((*en),1);
      int f2=Tetra::FofV((*en),2);

      if  (((*ti)->IsBorderF(f0))||((*ti)->IsBorderF(f1))||((*ti)->IsBorderF(f2)))
        _dummyV.Flags() |= LINK_V0;
         
      ti++;
      en++;
    }

    //and at the end I verify if the intersection is equal to the star of the edge
    bool correct=true;
    ti=_Sets.v1.begin(); 
    en=_Sets.indexv1.begin();

    while (ti!=_Sets.v1.end())
    {
      vt0=(*ti)->V(0);
      vt1=(*ti)->V(1);
      vt2=(*ti)->V(2);
      vt3=(*ti)->V(3);

      if ((vt0->Flags()& LINK_V0)&&(!(vt0->Flags()& LINK_EE)))
        correct=false;
      else
      if ((vt1->Flags()& LINK_V0)&&(!(vt1->Flags()& LINK_EE)))
        correct=false;
      else
      if ((vt2->Flags()& LINK_V0)&&(!(vt2->Flags()& LINK_EE)))
        correct=false;
      else
      if ((vt3->Flags()& LINK_V0)&&(!(vt3->Flags()& LINK_EE)))
        correct=false;

       //dummy vertex control
      int f0=Tetra::FofV((*en),0);
      int f1=Tetra::FofV((*en),1);
      int f2=Tetra::FofV((*en),2);

      if  (((*ti)->IsBorderF(f0))||((*ti)->IsBorderF(f1))||((*ti)->IsBorderF(f2)))
        if ((_dummyV.Flags()& LINK_V0)&&(!(_dummyV.Flags()& LINK_EE)))
          correct=false;

      if (!correct)
      {
        VertexType::DeleteUserBit(LINK_EE);
        VertexType::DeleteUserBit(LINK_V1);
        VertexType::DeleteUserBit(LINK_V0);  
        lkv++;
        return (false);
      }
      en++;
      ti++;
    }
    VertexType::DeleteUserBit(LINK_EE);
    VertexType::DeleteUserBit(LINK_V1);
    VertexType::DeleteUserBit(LINK_V0);
    return true;
}

///verify the flip conditions for a collapse
bool _FlipCondition(PosType pos,ScalarType alfa)
{	
  int edge=pos.E();
  VertexType *ve0=pos.T()->V(Tetra::VofE(edge,0));
	VertexType *ve1=pos.T()->V(Tetra::VofE(edge,1));
	CoordType oldpos0;
  CoordType oldpos1;
  CoordType newpos=((ve0->P()*alfa)+(ve1->P()*(1.f-alfa)));

  vector< TetraType *>::iterator ti=_Sets.no_E.begin();

  //verification
  oldpos0 = ve0->P();
	oldpos1 = ve1->P();

  //assegning new position
  ve0->P() =newpos;
	ve1->P() =newpos;

	while (ti!=_Sets.no_E.end())
    {
			assert(!(*ti)->IsD());
      assert((((*ti)->V(0)==ve0)||((*ti)->V(1)==ve0)||((*ti)->V(2)==ve0)||((*ti)->V(3)==ve0))^
            (((*ti)->V(0)==ve1)||((*ti)->V(1)==ve1)||((*ti)->V(2)==ve1)||((*ti)->V(3)==ve1)));

      Tetra3<ScalarType> T=Tetra3<ScalarType>();
      T.P0(0)=(*ti)->V(0)->cP();
      T.P1(0)=(*ti)->V(1)->cP();
      T.P2(0)=(*ti)->V(2)->cP();
      T.P3(0)=(*ti)->V(3)->cP();
      //flip comes if volume is less or equal to zero
			if (T.ComputeVolume()<=0)
			{	
        ve0->P()=oldpos0;
				ve1->P()=oldpos1;
        flip++;
  			return false;
			}
			ti++;	
		  }

  //reset initial value
	ve0->P()=oldpos0;
	ve1->P()=oldpos1;

	return true;
}

///update the normal of the modified tetrahedrons ond the normal of the vertex that remain after collapse
void _SetNormal(VertexType* v)
{
 if (TetraType::HasTetraNormal())
 {
  VTIterator<TetraType> VTi=VTIterator<TetraType>(v->VTb(),v->VTi());
  while (!VTi.End())
  {
    VTi.Vt()->ComputeNormal();
    VTi++;
  }
 }
  if (VertexType::HasNormal())
   _UN.PerVertex(v);
}
public:


///Return the aspect Ratio media of the tetrahedrons
///that share the adge to collapse
ScalarType AspectRatioCollapsed(PosType p)
{
  PosL pos=PosL(p.T(),p.F(),p.E(),p.V());
  pos.Reset();
  int num=0;
  ScalarType ratio_media=0.f;
  while(!pos.end())
  {
    ratio_media+=pos.T()->AspectRatio();
    pos.NextT();
    num++;
  }
  ratio_media=ratio_media/num;
  return (ratio_media);
}


///check the link conditions for the collapse indicated by pos
bool  CheckPreconditions(PosType pos,ScalarType alfa)
{	
  VertexType *v0=pos.T()->V(Tetra::VofE(pos.E(),0));
	VertexType *v1=pos.T()->V(Tetra::VofE(pos.E(),1));
  //if the two vertices are of border and the edge is not a border edge
  //we can do it.
  bool border0=v0->IsB();
  bool border1=v1->IsB();
  bool bordere=_Topo.IsExternEdge(pos.T(),pos.E());

  //first case vertex external and edge internal
	if ((border0 && border1)&&(!bordere))
  {
    c0++;
			return false;
  }
	else
  //if both vertex are internal so is enougth to verify flip conditions
	if ((!border0) && (!border1))
			return (_FlipCondition(pos,alfa));
  else
  //if the edge is internal is enougth to verify link condition on vertex
  if (!bordere)
      return((_FlipCondition(pos,alfa))&&(_LinkConditionsV()));
  else
  //at the end if trh edge is on the border we must verify also with the complete test
		return ((_FlipCondition(pos,alfa))&&(_LinkConditionsV())&&(_LinkConditionsE(pos))&&(_LinkConditionsF(pos)));
   //return false;
}


///Modify pos and alfa to obtain the collapse that minimize the error
void BestCollapse(PosType &pos,ScalarType &alfa,int nsteps)
{
  bool ext_v0=(pos.T()->V(Tetra::VofE(pos.E(),0)))->IsB();
  bool ext_v1=(pos.T()->V(Tetra::VofE(pos.E(),1)))->IsB();

   if ((ext_v0)&&(!ext_v1))
      alfa=1.f;
   else
   if ((!ext_v0)&&(ext_v1))
      alfa=0.f;
   else
   if ((!ext_v0)&&(!ext_v1))
     alfa=0.5f;
   else
   if ((ext_v0)&&(ext_v1))//both are external vertex
   {
    /*alfa=1.f;*/
    ScalarType step=1.f/(nsteps-1);
    ScalarType best_error=1000000.f;
    ScalarType Vol_Original=_VolumeUnion();
    for (int i=0;i<nsteps;i++)
    {
      ScalarType alfatemp=step*((double)i);
      //the error is the absolute value of difference of volumes
      ScalarType error=fabs(Vol_Original-_VolumeSimulateCollapse(pos,alfa));
      if(error<best_error)
      {
       alfa=alfatemp;
       best_error=error;
      }
    }
   }
}

///finds sets used for all test in edge collapse
void FindSets(vcg::tetra::Pos<TetraType> pos)
{
 
  _Sets.v0.clear();
  _Sets.indexv0.clear();
  _Sets.v1.clear();
  _Sets.indexv1.clear();
  _Sets.v0_U_v1.clear();
  _Sets.E.clear();
  _Sets.no_E.clear();
  _Sets.indexE.clear();
  int size=40;
  _Sets.v0.reserve(size);
  _Sets.indexv0.reserve(size);
  _Sets.v1.reserve(size);
  _Sets.indexv1.reserve(size);
  _Sets.v0_U_v1.reserve(size*2);
  _Sets.no_E.reserve(size*2);
  _Sets.E.reserve(size);
  _Sets.indexE.reserve(size);
  
  int edge =pos.E();

  VertexType *ve0=pos.T()->V(Tetra::VofE(edge,0));
  VertexType *ve1=pos.T()->V(Tetra::VofE(edge,1));
    	
	// put all tetrahedrons in the first one vector and in the union
  VTIterator<TetraType> vf0(ve0->VTb(),ve0->VTi());
  while (!vf0.End())
  {
    //set of ve0
    _Sets.v0.push_back(vf0.Vt());
    _Sets.indexv0.push_back(vf0.Vi());
    //set of union
    _Sets.v0_U_v1.push_back(vf0.Vt());
    //set of union minus intersection
    if ((vf0.Vt()->V(0)!=ve1)&&(vf0.Vt()->V(1)!=ve1)&&(vf0.Vt()->V(2)!=ve1)&&(vf0.Vt()->V(3)!=ve1))
		  _Sets.no_E.push_back(vf0.Vt());
    vf0++;
	}

  //second vertex iteration
  vf0.Vt()=ve1->VTb();
  vf0.Vi()=ve1->VTi();
  
  while (!vf0.End())
  {
    //set of ve1
    _Sets.v1.push_back(vf0.Vt());
    _Sets.indexv1.push_back(vf0.Vi());
    //set of union
    _Sets.v0_U_v1.push_back(vf0.Vt());
    //set of union minus intersection 
    if ((vf0.Vt()->V(0)!=ve0)&&(vf0.Vt()->V(1)!=ve0)&&(vf0.Vt()->V(2)!=ve0)&&(vf0.Vt()->V(3)!=ve0))
		  _Sets.no_E.push_back(vf0.Vt());
    vf0++;
	}

  //erase duplicated tetrahedrons from the union set
  sort(_Sets.v0_U_v1.begin(),_Sets.v0_U_v1.end());
  unique(_Sets.v0_U_v1.begin(),_Sets.v0_U_v1.end());

  //now compute the intersection
  PosLType PL(pos.T(),pos.F(),pos.E(),pos.V());

    //mark the vertex on the edge
    while (!PL.LoopEnd())
    {
      _Sets.E.push_back(PL.T());
      _Sets.indexE.push_back(PL.E());
      PL.NextT();
    }
}

///do the collapse on the edge in postype p
void DoCollapse(PosType p,ScalarType alfa)
{
  VertexType *v=p.T()->V(p.V());
  assert(p.T()->HasVTAdjacency());
  assert((alfa>=0)&&(alfa<=1.f));
  _Collapse(p,alfa);
  _SetNormal(v);
}


};
}//end namespace
}//end namespace
#endif