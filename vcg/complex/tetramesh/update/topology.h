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

Revision 1.1  2004/16/04 14:32  pietroni
Initial commit


****************************************************************************/
#ifndef __VCG_TETRA_UPDATE_TOPOLOGY
#define __VCG_TETRA_UPDATE_TOPOLOGY
#include <algorithm>

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
		sort(vertex,vertex+3);
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
		return !(((*this)>f)&&((*this)!=f));
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
class UpdateTopology
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

private: 
  VertexContainer* _vert;
  TetraContainer* _tetra;

public:
  ///defaul constructor
  UpdateTopology(VertexContainer *v,TetraContainer *t)
	{
		_vert=v;
		_tetra=t;
	}
/***********************************************/
/** @Vertex-Tetrahedron Topology Funtions
**/
//@{

  ///create the VT topology for tetrahedrons that are into containers
void VTTopology()
	{
			vertex_iterator v;
			tetra_iterator  t;

			ClearVTTopology();

			for(t=_tetra.begin();t!=_tetra.end();++t)
			if( ! (*t).IsD())
				for(int j=0;j<4;++j)
				{
					(*t).tv[j] = (*t).V(j)->Fp();
					(*t).zv[j] = (*t).V(j)->Zp();
					(*t).V(j)->Fp() = &(*t);
					(*t).V(j)->Zp() = j;
				}
	
	}

  /// clear the Vertex-Tetra topology
 void ClearVTTopology()
	{
		vertex_iterator v;
		for(v=_vert->begin();v!=_vert->end();++v)
			{
				v->Fp() = 0;
				v->Zp() = 0;
			}

		tetra_iterator   t;
		for(t=_tetra->begin();t!=_tetra->end();++t)	
			for(int j=0;j<4;++j)
				{
					(*t).TV(j) = 0;
					(*t).ZV(j) = 0;
				}
}

///erase one tetrahedron from VTTopology of all his vertices

void DetachVTTopology(TetraType *t)
	{	
		int i;
		for(i=0;i<4;i++)
			DetachVTTopology(t->V(i),t);
	}

///erase one tetrahedron from VTTopology of one specified vertex
void DetachVTTopology(VertexType *v,TetraType *t)
	{	
		TetraType *lastt;
		int lastz;
		EdgePosT<TetraType> Et(v->Fp(),v->Zp());
		if (Et.t==t)
		{
			v->Fp()=(TetraType *)t->tv[v->Zp()];
			v->Zp()=t->zv[v->Zp()];
		}
		else
		{	
			lastz=Et.z;
			while((Et.t!=t)&&(Et.t!=NULL))
			{	
				lastz=Et.z;
				lastt=Et.t;
				Et.NextT();
			}
			//in the list of the vertex v must be present the 
			//tetrahedron that you want to detach
			assert(Et.t!=NULL);
			lastt->tv[lastz]=Et.t->tv[Et.z];
			lastt->zv[lastz]=Et.t->zv[Et.z];
		}
	}

///insert the tetrahedron t in VT topology for vertex v of index z
void InsertVTTopology(VertexType *v,int z,TetraType *t)
	{
		if( ! (*t).IsD())
				{
					(*t).tv[z] = v->Fp();
					(*t).zv[z] = v->Zp();
					 v->Fp() = &(*t);
					 v->Zp() = z;
				}
}


///insert the tetrahedron t in VT topology for all his vertices
void InsertVTTopology(TetraType *t)
	{	
		assert(!t->IsD());
		int k=0;
		for (k=0;k<4;k++)
			{
					assert(!t->V(k)->IsD());
					InsertVTTopology(t->V(k),k,t);
			}		
}
/*@}*/
/***********************************************/
/** @Tetrahedron-Tetrahedron Topology Funtions
**/
//@{
///Build the Tetrahedron-Tetrahedron Topology (by Face)
void TTTopology()
	{
		vector <Facet<VertexType,TetraType> > VF;
		VertexType* v0;
		VertexType* v1;
		VertexType* v2;
		
		for (TetraIterator ti=_tetra->begin();ti!=_tetra->end();ti++)
    {	
		 if (!(*ti).IsD())
     {
			(*ti).Z(0)=0;
			(*ti).Z(1)=1;
			(*ti).Z(2)=2;
			(*ti).Z(3)=3;
		 	(*ti).T(0)=(&(*ti));
			(*ti).T(1)=(&(*ti));
			(*ti).T(2)=(&(*ti));
			(*ti).T(3)=(&(*ti));
			
      v0=(*ti).V(Tetra4<double>::VofF(0,0));
      v1=(*ti).V(Tetra4<double>::VofF(0,1));
      v2=(*ti).V(Tetra4<double>::VofF(0,2));
		
			VF.push_back(Facet<VertexType,TetraType>(v0,v1,v2,&(*ti),0));

      v0=(*ti).V(Tetra4<double>::VofF(1,0));
      v1=(*ti).V(Tetra4<double>::VofF(1,1));
      v2=(*ti).V(Tetra4<double>::VofF(1,2));

			VF.push_back(Facet<VertexType,TetraType>(v0,v1,v2,&(*ti),1));

      v0=(*ti).V(Tetra4<double>::VofF(2,0));
      v1=(*ti).V(Tetra4<double>::VofF(2,1));
      v2=(*ti).V(Tetra4<double>::VofF(2,2));
			
			VF.push_back(Facet<VertexType,TetraType>(v0,v1,v2,&(*ti),2));

      v0=(*ti).V(Tetra4<double>::VofF(3,0));
      v1=(*ti).V(Tetra4<double>::VofF(3,1));
      v2=(*ti).V(Tetra4<double>::VofF(3,2));

			VF.push_back(Facet<VertexType,TetraType>(v0,v1,v2,&(*ti),3));
     }
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
					t0->T(faceindex0)=(t1);
					t1->T(faceindex1)=(t0);
					t0->Z(faceindex0)=(faceindex1);
					t1->Z(faceindex1)=(faceindex0);
					i++;
				}
				
		}
  }

///Test the Tetrahedron-Tetrahedron Topology (by Face)
void TestTTTopology()
	{
    int i;
		for (TetraIterator ti=_tetra->begin();ti!=_tetra->end();ti++)
    {	
			for (i=0;i<4;i++)
			{
				if ((!(*ti).IsD())&& ((*ti).T(i)!=&(*ti)))
					{	
					  assert( ((((*ti).T(i))->T((*ti).Z(i)))==&(*ti)));
            
            VertexType	*v0=(*ti).V(Tetra4<double>::VofF(i,0));
            VertexType	*v1=(*ti).V(Tetra4<double>::VofF(i,1));
            VertexType	*v2=(*ti).V(Tetra4<double>::VofF(i,2));
						
						TetraType *t1=(TetraType*)(*ti).T(i);
						int z1=(*ti).Z(i);
            
            VertexType	*vo0=(*ti).V(Tetra4<double>::VofF(z1,0));
            VertexType	*vo1=(*ti).V(Tetra4<double>::VofF(z1,1));
            VertexType	*vo2=(*ti).V(Tetra4<double>::VofF(z1,2));

						assert((v0!=v1)&&(v0!=v2)&&(v1!=v2));
						assert((vo0!=vo1)&&(vo0!=vo2)&&(vo1!=vo2));

						assert ((v0==vo0)||(v0==vo1)||(v0==vo2));
						assert ((v1==vo0)||(v1==vo1)||(v1==vo2));
						assert ((v2==vo0)||(v2==vo1)||(v2==vo2));
					}
			}
    }
 		
	}

void setExternalVertices()
  {	
    
		TetraIterator tt;
		int i;
    for (tt=_tetra.begin();tt<_tetra.end();++tt)
    {
     
			for(i=0;i<4;i++)
			{
				if ((*tt).IsBorderF(i))
				{
					(*tt).FV(i,0)->SetB();
					(*tt).FV(i,1)->SetB();
					(*tt).FV(i,2)->SetB();
				}
        else
        {
          (*tt).FV(i,0)->SetB();
					(*tt).FV(i,1)->SetB();
					(*tt).FV(i,2)->SetB();
        }
      }
       
    }
	}

/*@}*/
}; // end class


/*@}*/
}	// End namespace
}	// End namespace


#endif
