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
#include<vcg\simplex\face\topology.h>
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
  typedef	typename TriMeshType::FaceIterator FaceIterator;
  /// The coordinate type
	typedef	typename FaceType::VertexType::CoordType CoordType;
  /// The scalar type
  typedef	typename TriMeshType::VertexType::ScalarType ScalarType;
  ///the container of tetrahedron type
  typedef typename TriMeshType::FaceContainer FaceContainer;
  ///the container of vertex type
  typedef typename TriMeshType::VertContainer VertContainer;
  ///half edge type
	typedef typename vcg::face::Pos<FaceType> PosType;
	/// vector of pos
	typedef typename std::vector<PosType> PosVec;
	///of VFIterator
	typedef typename vcg::face::VFIterator<FaceType>  VFI;
	/// vector of VFIterator
	typedef typename std::vector<vcg::face::VFIterator<FaceType> > VFIVec;



	/// Default Constructor
	EdgeCollapse()
		{
		};

  ~EdgeCollapse()
		{
		};


  struct Edge{
			VertexType* v0,v1;
      Edge(	const VertexType*& a,const VertexType*& b){
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

std::map<Edge,char> EdgeMark;

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
//bool LinkCondition(vcg::face::Pos<FaceType> pos) 
//{ 	
//    const int LINK_V0 = VertexType::NewBitFlag();
//		const int LINK_V1 = VertexType::NewBitFlag();
//    const int LINK_EE = VertexType::NewBitFlag();
//    
//		const int NOT_LINKED = ~(LINK_V0 | LINK_V1 | LINK_EE);
//    
//    VertexType *ve0=pos.f->V(pos.z);
//    VertexType *ve1=pos.f->V((pos.z+1)%3);
//    int edge =pos.z;
//
//		VFIterator<FaceType> vf0(ve0->VFb(),ve0->VFi());
//		// Clear visited and adj flag for all vertices adj to v0;
//    while (!vf0.End())
//    {
//			vf0.f->V(0)->Flags() &= NOT_LINKED;
//			vf0.f->V(1)->Flags() &= NOT_LINKED;
//      vf0.f->V(2)->Flags() &= NOT_LINKED;
//      vf0++;
//		}
//
//		VFIterator<FaceType> vf1(ve1->VFb(),ve1->VFi());
//		// Clear visited and adj flag for all vertices adj to v0;
//    while (!vf1.End())
//    {
//			vf1.f->V(0)->Flags() &= NOT_LINKED;
//			vf1.f->V(1)->Flags() &= NOT_LINKED;
//      vf1.f->V(2)->Flags() &= NOT_LINKED;
//      vf1++;
//		}
//
//    vf0.f=ve0->VFb();
//    vf0.f=ve0->VFi();
//		// Mark vertices of v0 
//		while (!vf0.End())
//    {
//      vf0.f->V(0)->Flags() |= (LINK_V0);
//      vf0.f->V(1)->Flags() |= (LINK_V0);
//      vf0.f->V(2)->Flags() |= (LINK_V0);
//      orMark(Edge(vf0.f->V(0),vf0.f->V(1)),LINK_V0);
//      orMark(Edge(vf0.f->V(1),vf0.f->V(2)),LINK_V0);
//      orMark(Edge(vf0.f->V(2),vf0.f->V(0)),LINK_V0);
//      vf0++;
//    }
//		
//    //mark the entities on the edge
//    VertexType* vt0=pos.f->V(edge);
//    VertexType* vt1=pos.f->V((edge+1)%3);
//    VertexType* vt2=pos.f->V((edge+2)%3);
//
//    vt0->Flags() |= (LINK_EE);
//    vt1->Flags() |= (LINK_EE);
//    vt2->Flags() |= (LINK_EE);
//    
//    
//    FaceType *opp=pos.f()->FFp(edge);
//    int eopp=pos.f()->FFi(edge);
//    
//    VertexType* vt3=opp.f->V((eopp+2)%3);
//    vt3->Flags() |= LINK_EE;
//
//    //mark the edges
//    orMark(Edge(vt0,vt1),LINK_EE);
//    orMark(Edge(vt0,vt2),LINK_EE);
//    orMark(Edge(vt1,vt2),LINK_EE);
//    orMark(Edge(vt0,vt3),LINK_EE);
//    orMark(Edge(vt1,vt3),LINK_EE);
//
//    //and at the end I verify if the intersection is equal to the star of the edge
//    vf1.f=ve1->VFb();
//    vf1.f=ve1->VFi();
//    bool correct=true;
//    while (!vf1.End())
//    {
//      vt0=vf1.f->V(0);
//      vt1=vf1.f->V(1);
//      vt2=vf1.f->V(2);
//
//      if ((vt0->Flags()& LINK_V0)&&(!(vt0->Flags()& LINK_EE)))
//        correct=false;
//      else
//      if ((vt1->Flags()& LINK_V0)&&(!(vt1->Flags()& LINK_EE)))
//        correct=false;
//      else
//      if ((vt2->Flags()& LINK_V0)&&(!(vt2->Flags()& LINK_EE)))
//        correct=false;
//      else
//			if ((isMarked(Edge(v0,v1),LINK_V0))&&(!isMarked(Edge(v0,v1),LINK_EE)))
//        correct=false;
//      else
//      if ((isMarked(Edge(v1,v2),LINK_V0))&&(!isMarked(Edge(v1,v2),LINK_EE)))
//        correct=false;
//      else
//      if ((isMarked(Edge(v2,v0),LINK_V0))&&(!isMarked(Edge(v2,v0),LINK_EE)))
//        correct=false;
//
//      if (!correct)
//      {
//        VertexType::DeleteBitFlag(LINK_V0);
//        VertexType::DeleteBitFlag(LINK_V1);
//        VertexType::DeleteBitFlag(LINK_EE);
//        return (false)
//      }
//      vf1++;
//    }
//    return true;
//    VertexType::DeleteBitFlag(LINK_V0);
//    VertexType::DeleteBitFlag(LINK_V1);
//    VertexType::DeleteBitFlag(LINK_EE);
//	}

	static VFIVec & AV0(){static VFIVec av0; return av0;}
	static VFIVec & AV1(){static VFIVec av1; return av1;}
	static VFIVec & AV01(){static VFIVec av01; return av01;}


	void FindSets(PosType p)
	{
		VertexType * v0 = p.V(0);
		VertexType * v1 = p.V(1);

		AV0().clear(); // Facce incidenti in v0 
		AV1().clear(); // Facce incidenti in v1 
		AV01().clear(); // Facce incidenti in v0 e v1
		
		VFI x;

		for( x.f = v0->VFp(), x.z = v0->VFi(); x.f!=0; x++)
		{
			
			int zv1 = -1;

			for(int j=0;j<3;++j) 
				if( x.f->V(j)==&*v1 )	{
					zv1 = j;
					break;
				}
			if(zv1==-1) 	AV0().push_back( x ); // la faccia x.f non ha il vertice v1 => e' incidente solo in v0
			else    			AV01().push_back( x );
		}

		for( x.f = v1->VFp(), x.z = v1->VFi(); x.f!=0; x++ )
		{
			int zv0 = -1;

			for(int j=0;j<3;++j)
				if( x.f->V(j)==&*v0 )	{
					zv0 = j;
					break;
				}
			if(zv0==-1)	AV1().push_back( x ); // la faccia x.f non ha il vertice v1 => e' incidente solo in v0
		}
}

	bool LinkConditions(PosType  pos){
		
		const int ADJ_1 = TriMeshType::VertexType::NewBitFlag();
		const int ADJ_E = TriMeshType::VertexType::NewBitFlag();
		//enum {ADJ_1= MeshType::VertexType::USER0,
		//	    ADJ_E= MeshType::VertexType::USER0<<1} ;
		const int ALLADJ = ADJ_1|ADJ_E;
		const int NOTALLADJ = ~(ADJ_1 | ADJ_E | TriMeshType::VertexType::VISITED);
		const int NOTALLADJ1 = ~(ADJ_E | TriMeshType::VertexType::VISITED);

		//EdgePosB<MeshType::face_type::face_base> x;
		vcg::face::VFIterator<FaceType> x;
		// Clear visited and adj flag for all vertices adj to v0;
		for(x.f = pos.V(0)->VFp(), x.z = pos.V(0)->VFi(); x.f!=0; x++ )	 	{
			x.f->V1(x.z)->Flags() &= NOTALLADJ;
			x.f->V2(x.z)->Flags() &= NOTALLADJ;
		}
		// Clear visited flag for all vertices adj to v1 and set them adj1 to v1;
		for(x.f = pos.V(1)->VFp(), x.z = pos.V(1)->VFi(); x.f!=0; x++ )	 	{
			x.f->V1(x.z)->Flags() &= NOTALLADJ1;
			x.f->V2(x.z)->Flags() &= NOTALLADJ1;
		}
		// Mark vertices adj to v1 as  ADJ_1 and  adj1 to v1;
		for(x.f = pos.V(1)->VFp(), x.z = pos.V(1)->VFi(); x.f!=0; x++ )	 	{
			if(x.f->V1(x.z)==pos.V(0)) x.f->V2(x.z)->Flags() |= ADJ_E | ADJ_1;
			else x.f->V2(x.z)->Flags() |= ADJ_1;
			if(x.f->V2(x.z)==pos.V(0)) x.f->V1(x.z)->Flags() |= ADJ_E | ADJ_1;
			else x.f->V1(x.z)->Flags() |= ADJ_1;
		}
			
		// compute the number of:
		int adj01=0;  // vertices adjacents to both v0 and v1 
		int adje=0;   // vertices adjacents to an egde (usually 2)
		for(x.f = pos.V(0)->VFp(), x.z = pos.V(0)->VFi(); x.f!=0; x++ )	 	{
			if(!x.f->V1(x.z)->IsV()) {
				x.f->V1(x.z)->SetV();
				if(x.f->V1(x.z)->Flags()&ADJ_1) ++adj01;
				if(x.f->V1(x.z)->Flags()&ADJ_E) ++adje;
			}
			if(!x.f->V2(x.z)->IsV()) {
				x.f->V2(x.z)->SetV();
				if(x.f->V2(x.z)->Flags()&ADJ_1) ++adj01;
				if(x.f->V2(x.z)->Flags()&ADJ_E) ++adje;
			}
		}

		//bool val=TopoCheck2();
		//if(val != (adj01==adje)) 		printf("Wrong topo %i %i\n",adj01,adje);
		TriMeshType::VertexType::DeleteBitFlag(ADJ_E);
		TriMeshType::VertexType::DeleteBitFlag(ADJ_1);
		
    return (adj01==adje);
	}



	int DoCollapse(PosType & c, Point3<ScalarType> p)
	{
		VFIVec::iterator i;
		int n_face_del =0 ;

		for(i=AV01().begin();i!=AV01().end();++i)
		{	
			FaceType  * f = ((*i).f);
			assert(f->V((*i).z) == c.V(0));

			vcg::face::VFDetach(*f,((*i).z+1)%3);
			vcg::face::VFDetach(*((*i).f),((*i).z+2)%3);
			((*i).f)->SetD();
			n_face_del++;
		}

		for(i=AV0().begin();i!=AV0().end();++i)
		{
			(*i).f->V((*i).z) = c.V(1);									 // In tutte le facce incidenti in v0, si sostituisce v0 con v1
			(*i).f->VFp((*i).z) = (*i).f->V((*i).z)->VFp(); // e appendo la lista di facce incidenti in v1 a questa faccia
			(*i).f->VFi((*i).z) = (*i).f->V((*i).z)->VFi();
			(*i).f->V((*i).z)->VFp() = (*i).f;
			(*i).f->V((*i).z)->VFi() = (*i).z;
		}


		TriMeshType::VertexType nv;
//*		c.v[1]->Merge(*c.v[0]);
		c.V(0)->SetD();
		return n_face_del;
	}

};

}
}
#endif 