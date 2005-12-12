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
  typedef	TETRA_MESH_TYPE TetraMeshType;
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
  typedef typename vcg::tetra::UpdateTetraTopology<VertexContainer,TetraContainer> Topology;
  ///the normal updater type
  typedef  typename vcg::tetra::UpdateNormals<TetraMeshType> UpdateNormals;


  /// Default Constructor
	EdgeCollapse()
		{
		};

  ~EdgeCollapse()
		{
		};

  private:

	typedef pair <int,int> FacePair;
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

  void clear()
  {
    v0.clear();
    v1.clear();
    v0_U_v1.clear();
    no_E.clear();
    E.clear();
	indexE.clear();
	indexv0.clear();
	indexv1.clear();
  }
};





static map<Edge,char> & _EdgeMark(){
	static map<Edge,char>  em;
	return em;
};

static map<Face,char> & _FaceMark(){
	static map<Face,char> fm;
	return fm;
}

static VertexType &_DummyV(){
	static VertexType _dv;
	return _dv;
}

static TetraSets &_Sets(){
	static TetraSets _s;
	return _s;
}



///select the 2 faces that does not share the edge
static FacePair _FindNoEdgeFace(TetraType *t,int edge)
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

#ifdef _DEBUG
static void _AssertingVolume(TetraType *t)
{
		//assert(t->ComputeVolume() >0);
		assert(vcg::ComputeVolume<TetraType>(*t)>0);
}
#endif


///collpse de edge specified by pos (the first vertex on edge remain)
static int _Collapse(PosType p,CoordType NewP)
{
	  int n_deleted=0;
	  vector<TetraType*> To_Del;
      VertexType *Vrem=(p.T()->V(Tetra::VofE(p.E(),0)));
      VertexType *Vdel=(p.T()->V(Tetra::VofE(p.E(),1)));
      //Vrem->P()=(Vrem->P()*alfa)+(Vdel->P()*(1.f-alfa));
      Vrem->P()=NewP;
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
				Topology::_AttachTTTopology(tleft,ileft,tright,iright);	
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
				Topology::DetachVTTopology(pos.T());

				//end setting the V-T topology	
				To_Del.push_back(pos.T());
				pos.NextT();
				n_deleted++;
//        tm.tn--;
      }

      //delting old tetrahedrons
      typename vector<TetraType*>::iterator ti;
      for (ti=To_Del.begin();ti<To_Del.end();ti++)
        (*ti)->SetD();

			//now I cycle on the tetrahedron that had the old vertex
			//reassegning the new one.

			VTIterator<TetraType> VTi(Vdel->VTb(),Vdel->VTi());
			while (!VTi.End())
			{	
				TetraType *T_Change=VTi.Vt();
				int index=VTi.Vi();
				//VTi++;
				//assegning the vertex that remain
				T_Change->V(index)=Vrem;
				Topology::DetachVTTopology(Vdel,T_Change);
				Topology::InsertVTTopology(Vrem,index,T_Change);
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
			return n_deleted;
}


static void orMarkE(Edge E,char M)
{
 typename map<Edge,char>::iterator EI;
 EI=_EdgeMark().find(E);
 if (EI==_EdgeMark().end())
  _EdgeMark().insert (pair<Edge,char>(E,M));
 else
   (*EI).second|=M;
}

static bool isMarkedE(Edge E,char M)
{
 typename map<Edge,char>::iterator EI;
 EI=_EdgeMark().find(E);
 if (EI==_EdgeMark().end())
    return false;
 else return (((*EI).second & M)!=0);
}

static void orMarkF(Face F,char M)
{
 typename map< Face,char>::iterator FI;
 FI=_FaceMark().find(F);
 if (FI==_FaceMark().end())
  _FaceMark().insert (pair<Face,char>(F,M));
 else
   (*FI).second|=M;
}

static bool isMarkedF(Face F,char M)
{
 typename map<Face,char>::iterator FI;
 FI=_FaceMark().find(F);
 if (FI==_FaceMark().end())
    return false;
 else return (((*FI).second & M)!=0);
}


///verify the link conditions on faces
static bool _LinkConditionsF(PosType pos) 
{ 	
    const int LINK_V0 = 0x00000001;
    const int LINK_EE = 0x00000002;
   
    _EdgeMark().clear();

		// Mark edges of ve0 
    typename vector< TetraType *>::iterator ti=_Sets().v0.begin(); 
    typename vector< char >::iterator en=_Sets().indexv0.begin();
	VertexType *v0=(*ti)->V(*en);
    while (ti!=_Sets().v0.end())
    {
	  assert(v0==(*ti)->V(*en));
      //put dummy face
      for (int f=0;f<3;f++)
      {
       int f_test=Tetra::FofV((*en),f);
        if  ((*ti)->IsBorderF(f_test))
        {
          orMarkF(Face((*ti)->V(Tetra::VofF(f_test,0)),(*ti)->V(Tetra::VofF(f_test,1)),&_DummyV()),LINK_V0);
          orMarkF(Face((*ti)->V(Tetra::VofF(f_test,1)),(*ti)->V(Tetra::VofF(f_test,2)),&_DummyV()),LINK_V0);
          orMarkF(Face((*ti)->V(Tetra::VofF(f_test,2)),(*ti)->V(Tetra::VofF(f_test,0)),&_DummyV()),LINK_V0);
        }
      }
     ti++;
     en++;
    }
    
    ti=_Sets().E.begin(); 
    en=_Sets().indexE.begin();
    //mark them as intersection
    while (ti!=_Sets().E.end())
    {
      //faces on the edge
      int f0=Tetra::FofE((*en),0);
      int f1=Tetra::FofE((*en),1);

      if  ((*ti)->IsBorderF(f0))
      {
          orMarkF(Face((*ti)->V(Tetra::VofF(f0,0)),(*ti)->V(Tetra::VofF(f0,1)),&_DummyV()),LINK_EE);
          orMarkF(Face((*ti)->V(Tetra::VofF(f0,1)),(*ti)->V(Tetra::VofF(f0,2)),&_DummyV()),LINK_EE);
          orMarkF(Face((*ti)->V(Tetra::VofF(f0,2)),(*ti)->V(Tetra::VofF(f0,0)),&_DummyV()),LINK_EE);
      }

      if  ((*ti)->IsBorderF(f1))
      {
          orMarkF(Face((*ti)->V(Tetra::VofF(f1,0)),(*ti)->V(Tetra::VofF(f1,1)),&_DummyV()),LINK_EE);
          orMarkF(Face((*ti)->V(Tetra::VofF(f1,1)),(*ti)->V(Tetra::VofF(f1,2)),&_DummyV()),LINK_EE);
          orMarkF(Face((*ti)->V(Tetra::VofF(f1,2)),(*ti)->V(Tetra::VofF(f1,0)),&_DummyV()),LINK_EE);
      }

      ti++;
      en++;
    }
   
    //and at the end I verify if the intersection is equal to the star of the edge
    ti=_Sets().v1.begin(); 
    en=_Sets().indexv1.begin();
	VertexType *v1=(*ti)->V(*en);
    while (ti!=_Sets().v1.end())
    {
	  assert(v1==(*ti)->V(*en));
      //dummy edges control
      for (int f=0;f<3;f++)
      {
       int f_test=Tetra::FofV((*en),f);
        if  ((*ti)->IsBorderF(f_test))
        {
          //control all the 3 edges
          Face f_test0=Face((*ti)->V(Tetra::VofF(f_test,0)),(*ti)->V(Tetra::VofF(f_test,1)),&_DummyV());
          Face f_test1=Face((*ti)->V(Tetra::VofF(f_test,1)),(*ti)->V(Tetra::VofF(f_test,2)),&_DummyV());
          Face f_test2=Face((*ti)->V(Tetra::VofF(f_test,2)),(*ti)->V(Tetra::VofF(f_test,0)),&_DummyV());
          if (((isMarkedF(f_test0,LINK_V0))&&(!isMarkedF(f_test0,LINK_EE)))||
            ((isMarkedF(f_test1,LINK_V0))&&(!isMarkedF(f_test1,LINK_EE)))||
            ((isMarkedF(f_test2,LINK_V0))&&(!isMarkedF(f_test2,LINK_EE))))
					{
//						FAIL::LKF();
            return false;
					}
        }
      }
      ti++;
      en++;
    }
    return true;
}


///verify the link conditions on edges
static bool _LinkConditionsE(PosType pos) 
{ 	
    const int LINK_V0 = 0x00000001;
    const int LINK_EE = 0x00000002;
   
    _FaceMark().clear();

		// Mark edges of ve0 
    typename vector< TetraType *>::iterator ti=_Sets().v0.begin(); 
    typename vector< char >::iterator en=_Sets().indexv0.begin();
    while (ti!=_Sets().v0.end())
    {
      //put dummy edge
      for (int f=0;f<3;f++)
      {
       int f_test=Tetra::FofV((*en),f);
        if  ((*ti)->IsBorderF(f_test))
        {
          orMarkE(Edge((*ti)->V(Tetra::VofF(f_test,0)),&_DummyV()),LINK_V0);
          orMarkE(Edge((*ti)->V(Tetra::VofF(f_test,1)),&_DummyV()),LINK_V0);
          orMarkE(Edge((*ti)->V(Tetra::VofF(f_test,2)),&_DummyV()),LINK_V0);
        }
      }
     ti++;
     en++;
    }
    
    ti=_Sets().E.begin(); 
    en=_Sets().indexE.begin();
    //mark them as intersection
    while (ti!=_Sets().E.end())
    {
      //faces on the edge
      int f0=Tetra::FofE((*en),0);
      int f1=Tetra::FofE((*en),1);

      if  ((*ti)->IsBorderF(f0))
      {
          orMarkE(Edge((*ti)->V(Tetra::VofF(f0,0)),&_DummyV()),LINK_EE);
          orMarkE(Edge((*ti)->V(Tetra::VofF(f0,1)),&_DummyV()),LINK_EE);
          orMarkE(Edge((*ti)->V(Tetra::VofF(f0,2)),&_DummyV()),LINK_EE);
      }

      if  ((*ti)->IsBorderF(f1))
      {
          orMarkE(Edge((*ti)->V(Tetra::VofF(f1,0)),&_DummyV()),LINK_EE);
          orMarkE(Edge((*ti)->V(Tetra::VofF(f1,1)),&_DummyV()),LINK_EE);
          orMarkE(Edge((*ti)->V(Tetra::VofF(f1,2)),&_DummyV()),LINK_EE);
      }

      ti++;
      en++;
    }
   
    //and at the end I verify if the intersection is equal to the star of the edge
    ti=_Sets().v1.begin(); 
    en=_Sets().indexv1.begin();
    while (ti!=_Sets().v1.end())
    {

      //dummy edges control
      for (int f=0;f<3;f++)
      {
       int f_test=Tetra::FofV((*en),f);
        if  ((*ti)->IsBorderF(f_test))
        {
          //control all the 3 edges
          Edge e_test0=Edge((*ti)->V(Tetra::VofF(f_test,0)),&_DummyV());
          Edge e_test1=Edge((*ti)->V(Tetra::VofF(f_test,1)),&_DummyV());
          Edge e_test2=Edge((*ti)->V(Tetra::VofF(f_test,2)),&_DummyV());
          if (((isMarkedE(e_test0,LINK_V0))&&(!isMarkedE(e_test0,LINK_EE)))||
            ((isMarkedE(e_test1,LINK_V0))&&(!isMarkedE(e_test1,LINK_EE)))||
            ((isMarkedE(e_test2,LINK_V0))&&(!isMarkedE(e_test2,LINK_EE))))
					{
//						FAIL::LKE();
            return false;
					}
        }
      }
      ti++;
      en++;
    }
    return true;
}

static bool _QuickConditions(PosType pos)
{	
  VertexType *v0=pos.T()->V(Tetra::VofE(pos.E(),0));
  VertexType *v1=pos.T()->V(Tetra::VofE(pos.E(),1));

  //if the two vertices are of border and the edge is not a border edge
  //we can do it.

  bool border0=v0->IsB();
  bool border1=v1->IsB();
  bool bordere=Topology::IsExternEdge(pos.T(),pos.E());

  //first case vertex external and edge internal
	if ((border0 && border1)&&(!bordere))
	{
		return false;
	}
 	else /// look if the 2 other faces that don't share the vertex are external on not
	{
	
	typename vector< TetraType *>::iterator ti=_Sets().E.begin(); 
	typename vector< char >::iterator en=_Sets().indexE.begin();
		//mark them as intersection
    while (ti!=_Sets().E.end())
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
		if (((pos.T()==tleft)&&(pos.T()==tright)))
		{
			return false;
		}
		ti++;
		en++;
	}		
	
	}
		return true;
}

///verify the link conditions on vertices
static bool _LinkConditionsV() 
{ 	
    const int LINK_V0 = VertexType::NewBitFlag();
	const int LINK_V1 = VertexType::NewBitFlag();
    const int LINK_EE = VertexType::NewBitFlag();
    
	const int NOT_LINKED = ~(LINK_V0 | LINK_V1 | LINK_EE);
    _DummyV().Flags() &= NOT_LINKED;

    VertexType *vt0;
    VertexType *vt1;
    VertexType *vt2;
    VertexType *vt3;
		

   typename vector< TetraType *>::iterator ti=_Sets().v0_U_v1.begin(); 
    
    //reset all link flags
    while (ti!=_Sets().v0_U_v1.end())
    {
      for(int i=0;i<4;i++)
        (*ti)->V(i)->Flags() &= NOT_LINKED;
      ti++;
    }


    //also in the ones that appartain to the edge
    typename  vector< char >::iterator en;
    ti=_Sets().E.begin(); 
    en=_Sets().indexE.begin();
    //reset all link flags for intersection and in the same
    //time mark them as intersection
    while (ti!=_Sets().E.end())
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
         _DummyV().Flags() |= LINK_EE;

      ti++;
      en++;
    }


		// Mark vertices of ve0 
    ti=_Sets().v0.begin();
    en=_Sets().indexv0.begin();

    while (ti!=_Sets().v0.end())
    {
      for(int i=0;i<4;i++)
        (*ti)->V(i)->Flags() |= LINK_V0;
      
      //dummy faces on the vertex
      int f0=Tetra::FofV((*en),0);
      int f1=Tetra::FofV((*en),1);
      int f2=Tetra::FofV((*en),2);

      if  (((*ti)->IsBorderF(f0))||((*ti)->IsBorderF(f1))||((*ti)->IsBorderF(f2)))
        _DummyV().Flags() |= LINK_V0;
         
      ti++;
      en++;
    }

    //and at the end I verify if the intersection is equal to the star of the edge
    bool correct=true;
    ti=_Sets().v1.begin(); 
    en=_Sets().indexv1.begin();

    while (ti!=_Sets().v1.end())
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
        if ((_DummyV().Flags()& LINK_V0)&&(!(_DummyV().Flags()& LINK_EE)))
          correct=false;

      if (!correct)
      {
        VertexType::DeleteBitFlag(LINK_EE);
        VertexType::DeleteBitFlag(LINK_V1);
        VertexType::DeleteBitFlag(LINK_V0); 
//				FAIL::LKV();
        return (false);
      }
      en++;
      ti++;
    }
    VertexType::DeleteBitFlag(LINK_EE);
    VertexType::DeleteBitFlag(LINK_V1);
    VertexType::DeleteBitFlag(LINK_V0);
    return true;
}

///verify the flip condition
static bool _FlipCondition(PosType pos,CoordType NewP)
{	
  int edge=pos.E();
  VertexType *ve0=pos.T()->V(Tetra::VofE(edge,0));
  VertexType *ve1=pos.T()->V(Tetra::VofE(edge,1));
  CoordType oldpos0;
  CoordType oldpos1;

  typename  vector< TetraType *>::iterator ti=_Sets().no_E.begin();

  //saving old position
  oldpos0 = ve0->P();
  oldpos1 = ve1->P();

  //assegning new position
  ve0->P() =NewP;
  ve1->P() =NewP;

	while (ti!=_Sets().no_E.end())
    {
			assert(!(*ti)->IsD());
			assert((((*ti)->V(0)==ve0)||((*ti)->V(1)==ve0)||((*ti)->V(2)==ve0)||((*ti)->V(3)==ve0))^
            (((*ti)->V(0)==ve1)||((*ti)->V(1)==ve1)||((*ti)->V(2)==ve1)||((*ti)->V(3)==ve1)));
			if (vcg::ComputeVolume<TetraType>(**ti)<=0)
			{	
//				FAIL::VOL();
				ve0->P()=oldpos0;
				ve1->P()=oldpos1;
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
static void _InitTetrahedronValues(VertexType* v)
{
 
  VTIterator<TetraType> VTi= VTIterator<TetraType>(v->VTb(),v->VTi());
  while (!VTi.End())
  {
    if (TetraType::HasTetraQuality())
    {
      VTi.Vt()->ComputeAspectRatio();
    }

    if (TetraType::HasTetraNormal())
    {
      VTi.Vt()->ComputeNormal();
    }

    ++VTi;
  }

  VTi.Vt()=v->VTb();
  VTi.Vi()=v->VTi();
  while (!VTi.End())
  {
    for (int i=0;i<4;i++)
    {
      if (VTi.Vt()->V(i)->IsB())
      {
        if (VertexType::HasNormal())
        UpdateNormals::PerVertex(VTi.Vt()->V(i));
      }
      
    }
      ++VTi;
  }

}

public:

/// clean everything
static void Reset(){
	_EdgeMark().clear();
	_FaceMark().clear();
	_Sets().clear();
	_DummyV().ClearFlags();
}
///Return the aspect Ratio media of the tetrahedrons
///that share the adge to collapse
static ScalarType AspectRatioCollapsed(PosType p)
{
  //PosL pos=PosL(p.T(),p.F(),p.E(),p.V());
   PosLoop<TetraType> pos=PosLoop<TetraType>(p.T(),p.F(),p.E(),p.V());
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


///check the topologycal preserving  conditions for the collapse indicated by pos
static bool  CheckPreconditions(PosType pos,CoordType NewP)
{	
   VertexType *v0=pos.T()->V(Tetra::VofE(pos.E(),0));
   VertexType *v1=pos.T()->V(Tetra::VofE(pos.E(),1));
  //if the two vertices are of border and the edge is not a border edge
  //we can do it.
  bool border0=v0->IsB();
  bool border1=v1->IsB();
  bool bordere=Topology::IsExternEdge(pos.T(),pos.E());
	if (!_QuickConditions(pos))
	{
		//FAIL::BOR();
		return false;
	}
 // //first case vertex external and edge internal
	//if ((border0 && border1)&&(!bordere))
	//{
	//	//FAIL::BOR();
	//	return false;
	//}
 	else
  //if both vertex are internal so is enougth to verify flip conditions
	if ((!border0) && (!border1))
	   return (_FlipCondition(pos,NewP));
	else
	//if the edge is internal is enougth to verify link condition on vertex
	if (!bordere)
      return((_FlipCondition(pos,NewP))&&(_LinkConditionsV()));
	else
	//at the end if trh edge is on the border we must verify also with the complete test
	  return ((_FlipCondition(pos,NewP))&&(_LinkConditionsV())&&(_LinkConditionsE(pos))&&(_LinkConditionsF(pos)));
   //return false;
}

///return the sum of volumes of the union of  stars on vertices (the original volume of tetrahedrons)
static ScalarType VolumeOriginal()
{
  typename  vector< TetraType *>::iterator ti=_Sets().v0_U_v1.begin(); 
  ScalarType vol=0;
  while (ti!=_Sets().v0_U_v1.end())
  {
    vol+=(*ti)->Volume();
    ti++;
  }
  return vol;
}

///Calculate the volume on the vertex resulting after collapse...
static ScalarType VolumeSimulateCollapse(PosType Pos,CoordType &newP)
{
  VertexType *Vrem=(Pos.T()->V(Tetra::VofE(Pos.E(),0)));
  VertexType *Vdel=(Pos.T()->V(Tetra::VofE(Pos.E(),1)));

  if (Vrem!=Pos.T()->V(Pos.V()))
    swap<VertexType*>(Vdel,Vrem);

  ScalarType vol=0;
  CoordType oldpos = Vrem->P();
 
//move vertex that remain in the new position
  Vrem->P() = newP;

  typename vector< TetraType *>::iterator ti=_Sets().no_E.begin(); 

  while (ti!=_Sets().no_E.end())
  {
 /*       Tetra3<ScalarType> T=Tetra3<ScalarType>();
        T.P0(0)=(*ti)->V(0)->cP();
        T.P1(0)=(*ti)->V(1)->cP();
        T.P2(0)=(*ti)->V(2)->cP();
        T.P3(0)=(*ti)->V(3)->cP();
       
		vol+=T.ComputeVolume(); */
//		vol+= vcg::ComputeVolume<TetraType>(*((Tetra3<ScalarType>*)&*ti));
					
  		vol+= vcg::ComputeVolume(**ti);
        ti++;
  }
  Vrem->P()=oldpos;
  return vol;
}

///finds sets used for all test in edge collapse
static void FindSets(vcg::tetra::Pos<TetraType> pos)
{
 
  _Sets().clear();
  int size=40;
  _Sets().v0.reserve(size);
  _Sets().indexv0.reserve(size);
  _Sets().v1.reserve(size);
  _Sets().indexv1.reserve(size);
  _Sets().v0_U_v1.reserve(size*2);
  _Sets().no_E.reserve(size*2);
  _Sets().E.reserve(size);
  _Sets().indexE.reserve(size);
  
  int edge =pos.E();

  VertexType *ve0=pos.T()->V(Tetra::VofE(edge,0));
  VertexType *ve1=pos.T()->V(Tetra::VofE(edge,1));
    	
	// put all tetrahedrons in the first one vector and in the union
  VTIterator<TetraType> vf0(ve0->VTb(),ve0->VTi());
  while (!vf0.End())
  {
    //set of ve0
    _Sets().v0.push_back(vf0.Vt());
    _Sets().indexv0.push_back(vf0.Vi());
    //set of union
    _Sets().v0_U_v1.push_back(vf0.Vt());
    //set of union minus intersection
    if ((vf0.Vt()->V(0)!=ve1)&&(vf0.Vt()->V(1)!=ve1)&&(vf0.Vt()->V(2)!=ve1)&&(vf0.Vt()->V(3)!=ve1))
		  _Sets().no_E.push_back(vf0.Vt());
    ++vf0;
	}

  //second vertex iteration
  vf0.Vt()=ve1->VTb();
  vf0.Vi()=ve1->VTi();
  
  while (!vf0.End())
  {
    //set of ve1
    _Sets().v1.push_back(vf0.Vt());
    _Sets().indexv1.push_back(vf0.Vi());
    //set of union
    _Sets().v0_U_v1.push_back(vf0.Vt());
    //set of union minus intersection 
    if ((vf0.Vt()->V(0)!=ve0)&&(vf0.Vt()->V(1)!=ve0)&&(vf0.Vt()->V(2)!=ve0)&&(vf0.Vt()->V(3)!=ve0))
		  _Sets().no_E.push_back(vf0.Vt());
    ++vf0;
	}

  //erase duplicated tetrahedrons from the union set
  sort(_Sets().v0_U_v1.begin(),_Sets().v0_U_v1.end());
  unique(_Sets().v0_U_v1.begin(),_Sets().v0_U_v1.end());

  //now compute the intersection
  PosLType PL(pos.T(),pos.F(),pos.E(),pos.V());

    //mark the vertex on the edge
    while (!PL.LoopEnd())
    {
      _Sets().E.push_back(PL.T());
      _Sets().indexE.push_back(PL.E());
      PL.NextT();
    }

}

///do the collapse on the edge in postype p
static int DoCollapse(PosType p,CoordType newP)
{
  VertexType *v=p.T()->V(p.V());
  assert(p.T()->HasVTAdjacency());
  int n_del=_Collapse(p,newP);
  _InitTetrahedronValues(v);
  return n_del;
}


};
}//end namespace
}//end namespace
#endif
