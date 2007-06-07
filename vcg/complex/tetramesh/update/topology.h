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
Revision 1.17  2006/12/03 14:56:30  ganovelli
*** empty log message ***

Revision 1.16  2006/06/29 13:07:33  ganovelli
Aggiunta superclasse UpdateTopology templated sui container e con funzioni sui container


Revision 1.1  2004/16/04 14:32  pietroni
Initial commit


****************************************************************************/
#ifndef __VCG_TETRA_UPDATE_TOPOLOGY
#define __VCG_TETRA_UPDATE_TOPOLOGY
#include <algorithm>
#include <vector>
#include <map>
#include <vcg\simplex\tetrahedron\pos.h>
using namespace std;
namespace vcg {
namespace tetra {
/** Class Facet.
    This is class for definition of a face of tethahedron
		@param STL_VERT_CONT (Template Parameter) Specifies the type of the vertices container any the vertex type.
 */


template < class VERT_TYPE , class TETRA_TYPE>
class Facet{

public:

	/// The vertex type 
	typedef VERT_TYPE MVTYPE;
	typedef TETRA_TYPE MTTYPE;
	

private:
	MTTYPE *Tr;
	int numface;
	MVTYPE * vertex[3];	

public:

	Facet(MVTYPE *v0,MVTYPE *v1,MVTYPE *v2,TETRA_TYPE * t,int index)
	{	
		vertex[0]=v0;
		vertex[1]=v1;
		vertex[2]=v2;
		if(vertex[0] > vertex[1])	std::swap(vertex[0], vertex[1]);
		if(vertex[1] > vertex[2])
		{
			std::swap(vertex[1], vertex[2]);
			if(vertex[0] > vertex[1])
			{
				std::swap(vertex[0], vertex[1]);
			}
		}
		Tr = t;
		numface = index;
	}

	inline const MVTYPE  * V(int index) const
	{
		return vertex[index];
	}

	TETRA_TYPE *getTetrahedron()
	{
		return Tr;
	}

	void setTetrahedron(TETRA_TYPE * t)
	{
		Tr=t;
	}

	inline bool operator == ( Facet const & f) const
	{ 
		return ((vertex[0]==f.V(0))&&(vertex[1]==f.V(1))&&(vertex[2]==f.V(2)));
	}

	inline bool operator != ( Facet const & f) const
	{ 
		return !((*this) == f);
	}

	inline bool operator > ( Facet const & f) const
	{ 
		
			if (vertex[0]!=f.V(0))
			{
				if (vertex[0]>f.V(0))
					return true;
				else 
					return false;
			}
			else
			if (vertex[1]!=f.V(1))
			{
				if (vertex[1]>f.V(1))
					return true;
				else 
					return false;
			}
			else
			if (vertex[2]!=f.V(2))
			{
				if (vertex[2]>f.V(2))
					return true;
				else 
					return false;
			}else 
				return false;
		
	}	

	inline bool operator < ( Facet const & f) const
	{ 
		return (!((*this)>f)&&((*this)!=f));
	}

	inline bool operator <= ( Facet const & f) const
	{ 
		return (((*this)<f)||((*this)==f));
	}

	inline bool operator >= ( Facet const & f) const
	{ 
		return (((*this)>f)||((*this)==f));
	}

	int getFaceIndex()const
	{
		return numface;
	}
	
};//end class
/** \addtogroup tetramesh */
/*@{*/

 /**  Class UpdateTopology.
 This is class for Topology of a tetrahedralmesh.
		@param STL_VERT_CONT (Template Parameter) Specifies the type of the vertices container any the vertex type.
				@param STL_TETRA_CONT (Template Parameter) Specifies the type of the tetrahedrons container any the tetrahedrons type.
 */
template  < class STL_VERT_CONT ,class STL_TETRA_CONT >
class UpdateTopologyBase
{

public:

	/// The vertex container
	typedef STL_VERT_CONT VertexContainer;

	/// The tethaedhron container
	typedef STL_TETRA_CONT TetraContainer;

	/// The vertex type 
	typedef typename STL_VERT_CONT::value_type VertexType;
	
	/// The tetrahedron type 
	typedef typename STL_TETRA_CONT::value_type TetraType;

	/// The type of vertex iterator
	typedef typename STL_VERT_CONT::iterator VertexIterator;

	/// The type of tetra iterator
	typedef typename STL_TETRA_CONT::iterator TetraIterator;

	/// The type of constant vertex iterator
	typedef typename STL_VERT_CONT::const_iterator const_VertexIterator;

	/// The type of constant face iterator
	typedef typename STL_TETRA_CONT::const_iterator const_TetraIterator;

public:
  
/***********************************************/
/** @Vertex-Tetrahedron Topology Funtions
**/
//@{


/// Create the VT topology for tetrahedrons that are into containers.
static void VTTopology( VertexContainer & vert, TetraContainer & tetra )
{
	ClearVTTopology( vert, tetra );
	for( TetraIterator t = tetra.begin(); t != tetra.end(); ++t )
		if( !(*t).IsD() )
			for( int j = 0; j < 4; ++j )
			{
				(*t).VTp(j) = (*t).V(j)->VTp();
				(*t).VTi(j) = (*t).V(j)->VTi();
				(*t).V(j)->VTp() = &(*t);
				(*t).V(j)->VTi() = j;
			}
}


/// Clear the vertex-tetra (VT) topology.
static void ClearVTTopology( VertexContainer & vert, TetraContainer & tetra )
{
	for( VertexIterator v = vert.begin(); v != vert.end(); ++v ) { v->VTp() = 0; v->VTi() = 0; }
	for( TetraIterator t = tetra.begin(); t != tetra.end(); ++t )
		if( ! (*t).IsD() )
			for( int j = 0; j < 4; ++j ) { (*t).VTp(j) = 0; (*t).VTi(j) = 0; }
}


/// Erase one tetrahedron from VTTopology of all his vertices.
static void DetachVTTopology( TetraType *t )
{
	if( ! (*t).IsD() )
		for( int i = 0; i < 4; i++ ) DetachVTTopology( t->V(i), t );
}


/// Erase one tetrahedron from VTTopology of one specified vertex.
static void DetachVTTopology( VertexType *v, TetraType *t )
{
	TetraType *lastt;
	int lastz;
	VTIterator<TetraType> Et( v->VTb(), v->VTi() );
	if( Et.Vt() == t )
	{
		v->VTb() = (TetraType *) t->VTp( v->VTi() );
		v->VTi() = t->VTi( v->VTi() );
	}
	else
	{
		lastz = Et.Vi();
		while( ( Et.Vt() != t ) && ( !Et.End() ) )
		{
			lastz = Et.Vi();
			lastt = Et.Vt();
			++Et;
		}
		/// In the list of the vertex v must be present the tetrahedron that you want to detach
		assert( Et.Vt() != NULL );
		lastt->VTp(lastz) = Et.Vt()->VTp( Et.Vi() );
		lastt->VTi(lastz) = Et.Vt()->VTi( Et.Vi() );
	}
}


/// Insert the tetrahedron t in VT topology for vertex v of index z.
static void InsertVTTopology( VertexType *v, int z, TetraType *t )
{
	if( ! (*t).IsD() )
	{
		t->VTp(z) = v->VTb();
		t->VTi(z) = v->VTi();
		v->VTb() = &(*t);
		v->VTi() = z;
	}
}


/// Insert the tetrahedron t in VT topology for all his vertices.
static void InsertVTTopology( TetraType *t )
{
	assert( !( t->IsD() ) );
	for( int k = 0; k < 4; k++ )
	{
		assert( !( t->V(k)->IsD() ) );
		InsertVTTopology( t->V(k), k, t );
	}		
}


/// Test the Tetrahedron-Tetrahedron (TT) topology (by face).
static void TestVTTopology( VertexContainer & vert, TetraContainer & tetra )
{
	int i;
	for( VertexIterator vi = vert.begin(); vi != vert.end(); vi++ )
		if( !(*vi).IsD() )
        {
			TetraType *nextT = vi->VTb();
			int nextI = vi->VTi();
            int oldI;
            while( nextT != NULL )
            {
				assert( ( nextT->V(nextI) == &(*vi) ) );
				oldI = nextI;
				nextI = nextT->VTi(nextI);
				nextT = nextT->VTp(oldI);
            }
		}
}


/*@}*/
/***********************************************/
/** @Tetrahedron-Tetrahedron Topology Funtions
**/
//@{
///Build the Tetrahedron-Tetrahedron Topology (by Face)
static void TTTopology(const VertexContainer &vert,TetraContainer &tetra)
	{
		vector <Facet<VertexType,TetraType> > VF;
		VertexType* v0;
		VertexType* v1;
		VertexType* v2;
		
		for (TetraIterator ti=tetra.begin();ti!=tetra.end();ti++)
 		 if (!(*ti).IsD())
     {
			(*ti).TTi(0)=0;
			(*ti).TTi(1)=1;
			(*ti).TTi(2)=2;
			(*ti).TTi(3)=3;
		 	(*ti).TTp(0)=(&(*ti));
			(*ti).TTp(1)=(&(*ti));
			(*ti).TTp(2)=(&(*ti));
			(*ti).TTp(3)=(&(*ti));
			
      v0=(*ti).V(Tetra::VofF(0,0));
      v1=(*ti).V(Tetra::VofF(0,1));
      v2=(*ti).V(Tetra::VofF(0,2));
		
			VF.push_back(Facet<VertexType,TetraType>(v0,v1,v2,&(*ti),0));

      v0=(*ti).V(Tetra::VofF(1,0));
      v1=(*ti).V(Tetra::VofF(1,1));
      v2=(*ti).V(Tetra::VofF(1,2));

			VF.push_back(Facet<VertexType,TetraType>(v0,v1,v2,&(*ti),1));

      v0=(*ti).V(Tetra::VofF(2,0));
      v1=(*ti).V(Tetra::VofF(2,1));
      v2=(*ti).V(Tetra::VofF(2,2));
			
			VF.push_back(Facet<VertexType,TetraType>(v0,v1,v2,&(*ti),2));

      v0=(*ti).V(Tetra::VofF(3,0));
      v1=(*ti).V(Tetra::VofF(3,1));
      v2=(*ti).V(Tetra::VofF(3,2));

			VF.push_back(Facet<VertexType,TetraType>(v0,v1,v2,&(*ti),3));
  	}
	sort(VF.begin(),VF.end());
	
	TetraType *t0;
	TetraType *t1;
	int faceindex0;
	int faceindex1;
	int j;
	unsigned int i;
		for (i=0;i<VF.size()-1;i++)
		{	
			j=i+1;
			if (VF[i]==VF[j])
				{	
					t0=VF[i].getTetrahedron();
					t1=VF[j].getTetrahedron();
					faceindex0=VF[i].getFaceIndex();
					faceindex1=VF[j].getFaceIndex();
					t0->TTp(faceindex0)=(t1);
					t1->TTp(faceindex1)=(t0);
					t0->TTi(faceindex0)=(faceindex1);
					t1->TTi(faceindex1)=(faceindex0);
					i++;
				}
				
		}
  }



///Connect trought Tetrahedron-Tetrahedron Topology t0 and t1 with faces i0 and i1
static void _AttachTTTopology(TetraType *t0,int i0,TetraType *t1,int i1) 
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

///Detach Tetrahedron-Tetrahedron Topology 
static void DetachTTTopology(TetraType *t) 
{
  assert(!t->IsD());
	int i;
	for(i=0; i < 4; ++i)
		  t->TTp(i)->TTp(t->TTi(i)) = t->TTp(i);
}


///Test the Tetrahedron-Tetrahedron Topology (by Face)
static void TestTTTopology(VertexContainer &vert,TetraContainer &tetra)
	{
    int i;
		for (TetraIterator ti=tetra.begin();ti!=tetra.end();ti++)
		if ((!(*ti).IsD()))
			for (i=0;i<4;i++)
			{
					{	
					  assert( ((((*ti).TTp(i))->TTp((*ti).TTi(i)))==&(*ti)));
            
            VertexType	*v0=(*ti).V(Tetra::VofF(i,0));
            VertexType	*v1=(*ti).V(Tetra::VofF(i,1));
            VertexType	*v2=(*ti).V(Tetra::VofF(i,2));
						
						TetraType *t1=(TetraType*)(*ti).TTp(i);
            assert (!t1->IsD());
						int z1=(*ti).TTi(i);
            
            VertexType	*vo0=(*t1).V(Tetra::VofF(z1,0));
            VertexType	*vo1=(*t1).V(Tetra::VofF(z1,1));
            VertexType	*vo2=(*t1).V(Tetra::VofF(z1,2));

						assert((v0!=v1)&&(v0!=v2)&&(v1!=v2));
						assert((vo0!=vo1)&&(vo0!=vo2)&&(vo1!=vo2));

						assert ((v0==vo0)||(v0==vo1)||(v0==vo2));
						assert ((v1==vo0)||(v1==vo1)||(v1==vo2));
						assert ((v2==vo0)||(v2==vo1)||(v2==vo2));
					}
			}
	}

///test if all and only the exernal vertex are set of border
static void TestExternalVertex(VertexContainer &vert,TetraContainer &tetra)
{
  TetraIterator ti;
  VertexIterator vi;

 typedef pair <VertexType*, bool> VertBoolPair;
  map<VertexType*, bool> Inserted;
	typename map<VertexType*, bool>::iterator MapIte;

  for (ti=tetra.begin();ti<tetra.end();ti++)
  {
    int i;
    if (!ti->IsD())
    {
    for (i=0;i<4;i++)
      if (ti->IsBorderF(i))
      {
        VertexType *v0=ti->V(Tetra::VofF(i,0));
        VertexType *v1=ti->V(Tetra::VofF(i,1));
        VertexType *v2=ti->V(Tetra::VofF(i,2));

        MapIte = Inserted.find(v0);
        if ( MapIte == Inserted.end( ) )
          Inserted.insert (VertBoolPair(v0,true));

        MapIte = Inserted.find(v1);
        if ( MapIte == Inserted.end( ) )
          Inserted.insert (VertBoolPair(v1,true));

        MapIte = Inserted.find(v2);
        if ( MapIte == Inserted.end( ) )
          Inserted.insert (VertBoolPair(v2,true));

        assert(!((v0->IsD())||(v1->IsD())||(v2->IsD())));
        assert ((v0->IsB())&&(v1->IsB())&&(v2->IsB()));
      }
    }
  }

  for (vi=vert.begin();vi<vert.end();vi++)
  {
    if (!vi->IsD())
    {
      if (vi->IsB())
      {
        MapIte = Inserted.find(&(*vi));
        //control if the extrenal vertex appartain to an external face
        assert ( MapIte != Inserted.end( ) );
      } 
    }
  }
}

///set the external vertex according to Tetra-Tetra topology
static void setExternalVertices(VertexContainer &vert,TetraContainer &tetra)
{	
    
		TetraIterator tt;
    VertexIterator vi;
		int i;
    for (vi=vert.begin();vi<vert.end();++vi)
        vi->ClearB();
    for (tt=tetra.begin();tt<tetra.end();++tt)
			if(!(*tt).IsD())
    {
			for(i=0;i<4;i++)
			{
				if ((*tt).IsBorderF(i))
				{
          (*tt).V(Tetra::VofF(i,0))->SetB();
					(*tt).V(Tetra::VofF(i,1))->SetB();
					(*tt).V(Tetra::VofF(i,2))->SetB();
				}
        
      }
       
    }
	}


/*@}*/

private:

struct _triV
{
VertexType *v[3];

 _triV(VertexType *v0,VertexType *v1,VertexType *v2)
	{	
    v[0]=v0;
    v[1]=v1;
    v[2]=v2;
    sort(v,v+3);
	}

	inline const VertexType  * V(int index) const
	{
		return v[index];
	}

	inline bool operator == ( _triV const & tv) const
	{ 
		return ((v[0]==tv.V(0))&&(v[1]==tv.V(1))&&(v[2]==tv.V(2)));
	}

	inline bool operator != ( _triV const & tv) const
	{ 
		return !((*this) == tv);
	}

	inline bool operator > ( _triV const & tv ) const
	{ 
		
			if (v[0]!=tv.V(0))
			{
				if (v[0]>tv.V(0))
					return true;
				else 
					return false;
			}
			else
			if (v[1]!=tv.V(1))
			{
				if (v[1]>tv.V(1))
					return true;
				else 
					return false;
			}
			else
			if (v[2]!=tv.V(2))
			{
				if (v[2]>tv.V(2))
					return true;
				else 
					return false;
			}else 
				return false;
		
	}	

	inline bool operator < (_triV const & tv) const
	{ 
		return !(((*this)>tv)&&((*this)!=tv));
	}

	inline bool operator <= (_triV const & tv) const
	{ 
		return (((*this)<tv)||((*this)==tv));
	}

	inline bool operator >= ( _triV const & tv) const
	{ 
		return (((*this)>tv)||((*this)==tv));
	}
};


public:
///this function is used to test if an edge is extern
static bool IsExternEdge(TetraType *t,int edge)
{
	std::vector < _triV > Faces;

  assert((t->HasTTAdjacency())||(t->HasVTAdjacency()));
  if ((!t->V(Tetra::VofE(edge,0))->IsB())||(!t->V(Tetra::VofE(edge,1))->IsB()))
    return (false);

  if (t->HasTTAdjacency())
  {
    PosLoop<TetraType> pl(t,Tetra::FofE(edge,0),edge,Tetra::VofE(edge,0));
    pl.Reset();
    //stops if one of faces incident to the edge is an extern face
    while ((!pl.LoopEnd())&&(!pl.T()->IsBorderF(Tetra::FofE(pl.E(),0)))&&(!pl.T()->IsBorderF(Tetra::FofE(pl.E(),1))))
      pl.NextT();
    if (pl.LoopEnd())
      return false;
    else 
      return true;
  }
  else
  { //using vt adiacency
    VertexType *v0=t->V(Tetra::VofE(edge,0));
    VertexType *v1=t->V(Tetra::VofE(edge,1));
    assert(v0!=v1);
    VTIterator<TetraType> Vti(v0->VTb(),v0->VTi());
    int num=0;
    Faces.clear();
    Faces.reserve(40);
    while (!Vti.End())
    {   
        //take the three faces incident on one vertex
        int f0=Tetra::FofV(Vti.Vi(),0);
        int f1=Tetra::FofV(Vti.Vi(),1);
        int f2=Tetra::FofV(Vti.Vi(),2);
        VertexType *vf0=Vti.Vt()->V(Tetra::VofF(f0,0));
        VertexType *vf1=Vti.Vt()->V(Tetra::VofF(f0,1));
        VertexType *vf2=Vti.Vt()->V(Tetra::VofF(f0,2));
        //if there is the edge then put the three vertex in the vector
        if ((vf0==v1)||(vf1==v1)||(vf2==v1))
        {
          Faces.push_back(_triV(vf0,vf1,vf2));
          num++;
        }
    }
    sort(Faces.begin(),Faces.end());
    //now look if one face is no shared from other tetrahedron
    //2 instances of same face in vector means it is internal face
    bool isExtern=false;
    typename std::vector < _triV >::iterator TVIo;
    typename std::vector < _triV >::iterator TVIn;
    TVIo=Faces.begin();
    TVIn=Faces.begin();
    TVIn++;
    int j=0;
    while (((*TVIo)==(*TVIn))&&(j<num))
    {
        //move 2 steps each iterator to frify each pair of faces
        TVIo++;
        TVIo++;
        TVIn++;
        TVIn++;
        j++;
        j++;
    }
    if (j>=num)
      return false;
    else
      return true;
  }

}
}; // end class

template <class TetraMeshType>
class UpdateTopology: public UpdateTopologyBase<typename TetraMeshType::VertexContainer,
																										typename TetraMeshType::TetraContainer>{
public:
		static void TTTopology(TetraMeshType & tmesh){
			UpdateTopologyBase<typename TetraMeshType::VertexContainer,typename TetraMeshType::TetraContainer>::
				TTTopology(tmesh.vert,tmesh.tetra);
		}
		static void VTTopology(TetraMeshType & tmesh){	
			UpdateTopologyBase<typename TetraMeshType::VertexContainer,typename TetraMeshType::TetraContainer>::
				VTTopology(tmesh.vert,tmesh.tetra);
		}


};

/*@}*/
}	// End namespace
}	// End namespace


#endif
