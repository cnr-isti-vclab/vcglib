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

#ifndef __VCG_TETRA_TRI_COLLAPSE
#define __VCG_TETRA_TRI_COLLAPSE


#include<vcg\simplex\face\pos.h>
#include<map>

/** \addtogroup trimesh */
/*@{*/
/// This Class is used for the edge collapse

namespace vcg{
  namespace tri{	

template <class TRI_MESH_TYPE> 
class EdgeCollapse
{
	public:
  /// The tetrahedral mesh type
  typedef	typename TRI_MESH_TYPE TriMeshType;
  /// The tetrahedron type
  typedef	typename TriMeshType::FaceType FaceType;
	/// The vertex type
	typedef	typename FaceType::VertexType VertexType;
  /// The vertex iterator type
  typedef	typename TriMeshType::VertexIterator VertexIterator;
  /// The tetra iterator type
  typedef	typename TriMeshType::FaceIterator TetraIterator;
  /// The coordinate type
	typedef	typename FaceType::VertexType::CoordType CoordType;
  /// The scalar type
  typedef	typename TriMeshType::VertexType::ScalarType ScalarType;
  ///the container of tetrahedron type
  typedef typename TriMeshType::FaceContainer FaceContainer;
  ///the container of vertex type
  typedef typename TriMeshType::VertexContainer VertexContainer;
 
  /// Default Constructor
	EdgeCollapse()
		{
		};

  ~EdgeCollapse()
		{
		};


  struct Edge{
			VertexType* v0,v1;
      EdgeMark(	const VertexType*& a,const VertexType*& b){
						assert(a!=b);
						if(a<b) 
							{v0=a;v1=b;}
						else
						{v1=a;v0=b;}
      }
			
			const bool operator <(const EdgeMark & e) const {
				return (v0==e.v0)?(v1<e.v1):(v0<e.v0);
			}

			const bool operator ==(const EdgeMark & e) const {
			return (v0==e.v0)&&(v1==e.v1);
			}

};

map<Edge,char> EdgeMark;

void orMark(Edge E,char M)
{
 map<Edge,char>::Iterator EI;
 EdgeMark.find(E);
 if (EI==EdgeMArk.end())
  EdgeMark.insert (Pair<Edge,char>(E,M));
 else
   (*EI).second()|=M;
}

bool isMarked(Edge E,char M)
{
 map<Edge,char>::Iterator EI;
 EdgeMark.find(E);
 if (EI==EdgeMArk.end())
    return false;
 else return ((*EI).second()&M);
}

///control link conditions for the collapse
bool LinkCondition(Pos<FaceType> pos) 
{ 	
    const int LINK_V0 = VertexType::NewBitFlag();
		const int LINK_V1 = VertexType::NewBitFlag();
    const int LINK_EE = VertexType::NewBitFlag();
    
		const int NOT_LINKED = ~(LINK_V0 | LINK_V1 | LINK_EE);
    
    VertexType *ve0=pos.f->V(pos.z);
    VertexType *ve1=pos.f->V((pos.z+1)%3);
    int edge =pos.z;

		VFIterator<FaceType> vf0(ve0->VFb(),ve0->VFi());
		// Clear visited and adj flag for all vertices adj to v0;
    while (!vf0.End())
    {
			vf0.f->V(0)->Flags() &= NOT_LINKED;
			vf0.f->V(1)->Flags() &= NOT_LINKED;
      vf0.f->V(2)->Flags() &= NOT_LINKED;
      vf0++;
		}

		VFIterator<FaceType> vf1(ve1->VFb(),ve1->VFi());
		// Clear visited and adj flag for all vertices adj to v0;
    while (!vf1.End())
    {
			vf1.f->V(0)->Flags() &= NOT_LINKED;
			vf1.f->V(1)->Flags() &= NOT_LINKED;
      vf1.f->V(2)->Flags() &= NOT_LINKED;
      vf1++;
		}

    vf0.f=ve0->VFb();
    vf0.f=ve0->VFi();
		// Mark vertices of v0 
		while (!vf0.End())
    {
      vf0.f->V(0)->Flags() |= (LINK_V0);
      vf0.f->V(1)->Flags() |= (LINK_V0);
      vf0.f->V(2)->Flags() |= (LINK_V0);
      orMark(Edge(vf0.f->V(0),vf0.f->V(1)),LINK_V0);
      orMark(Edge(vf0.f->V(1),vf0.f->V(2)),LINK_V0);
      orMark(Edge(vf0.f->V(2),vf0.f->V(0)),LINK_V0);
      vf0++;
    }
		
    //mark the entities on the edge
    VertexType* vt0=pos.f->V(edge);
    VertexType* vt1=pos.f->V((edge+1)%3);
    VertexType* vt2=pos.f->V((edge+2)%3);

    vt0->Flags() |= (LINK_EE);
    vt1->Flags() |= (LINK_EE);
    vt2->Flags() |= (LINK_EE);
    
    
    FaceType *opp=pos.f()->FFp(edge);
    int eopp=pos.f()->FFi(edge);
    
    VertexType* vt3=opp.f->V((eopp+2)%3);
    vt3->Flags() |= LINK_EE;

    //mark the edges
    orMark(Edge(vt0,vt1),LINK_EE);
    orMark(Edge(vt0,vt2),LINK_EE);
    orMark(Edge(vt1,vt2),LINK_EE);
    orMark(Edge(vt0,vt3),LINK_EE);
    orMark(Edge(vt1,vt3),LINK_EE);

    //and at the end I verify if the intersection is equal to the star of the edge
    vf1.f=ve1->VFb();
    vf1.f=ve1->VFi();
    bool correct=true;
    while (!vf1.End())
    {
      vt0=vf1.f->V(0);
      vt1=vf1.f->V(1);
      vt2=vf1.f->V(2);

      if ((vt0->Flags()& LINK_V0)&&(!(vt0->Flags()& LINK_EE)))
        correct=false;
      else
      if ((vt1->Flags()& LINK_V0)&&(!(vt1->Flags()& LINK_EE)))
        correct=false;
      else
      if ((vt2->Flags()& LINK_V0)&&(!(vt2->Flags()& LINK_EE)))
        correct=false;
      else
			if ((isMarked(Edge(v0,v1),LINK_V0))&&(!isMarked(Edge(v0,v1),LINK_EE)))
        correct=false;
      else
      if ((isMarked(Edge(v1,v2),LINK_V0))&&(!isMarked(Edge(v1,v2),LINK_EE)))
        correct=false;
      else
      if ((isMarked(Edge(v2,v0),LINK_V0))&&(!isMarked(Edge(v2,v0),LINK_EE)))
        correct=false;

      if (!correct)
      {
        VertexType::DeleteBitFlag(LINK_V0);
        VertexType::DeleteBitFlag(LINK_V1);
        VertexType::DeleteBitFlag(LINK_EE);
        return (false)
      }
      vf1++;
    }
    return true;
    VertexType::DeleteBitFlag(LINK_V0);
    VertexType::DeleteBitFlag(LINK_V1);
    VertexType::DeleteBitFlag(LINK_EE);
	}


};

}
}
