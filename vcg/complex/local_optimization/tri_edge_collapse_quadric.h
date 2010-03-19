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
Revision 1.14  2007/03/22 11:07:16  cignoni
Solved an issue related to different casting double-float between gcc 3 and gcc 4

Revision 1.13  2007/02/25 09:20:10  cignoni
Added Rad to the NormalThr Option and removed a bug in multiple exectuion of non optimal simplification (missing an isD check)

Revision 1.12  2007/01/19 09:13:14  cignoni
Added Finalize() method to the interface, corrected minor bugs on border preserving and postsimplification cleanup
Avoided double make_heap (it is done only in the local_optimization init)

Revision 1.11  2006/10/15 07:31:21  cignoni
typenames and qualifiers for gcc compliance

Revision 1.10  2006/10/09 20:12:55  cignoni
Heavyly restructured for meshlab inclusion. Now the access to the quadric elements are mediated by a static helper class.

Revision 1.9  2006/10/07 17:20:25  cignoni
Updated to the new style face->Normal() becomes Normal(face)

Revision 1.8  2005/10/02 23:19:36  cignoni
Changed the sign of the priority of a collapse. Now it is its the error as it should (and not -error)

Revision 1.7  2005/04/14 11:35:07  ponchio
*** empty log message ***

Revision 1.6  2005/01/19 10:35:28  cignoni
Better management of symmetric/asymmetric edge collapses

Revision 1.5  2004/12/10 01:07:15  cignoni
Moved param classes inside; added support for optimal placement and symmetric; added update heap also here (not only in the base class)

Revision 1.4  2004/11/23 10:34:23  cignoni
passed parameters by reference in many funcs and gcc cleaning

Revision 1.3  2004/10/25 07:07:56  ganovelli
A vcg.::Pos was used to implement the collapse type. CHanged
to vcg::Edge

Revision 1.2  2004/09/29 17:08:16  ganovelli
corrected error in -error (see localoptimization)


****************************************************************************/




#ifndef __VCG_TRIMESHCOLLAPSE_QUADRIC__
#define __VCG_TRIMESHCOLLAPSE_QUADRIC__

#include<vcg/math/quadric.h>
#include<vcg/simplex/face/pos.h>
#include<vcg/complex/trimesh/update/flag.h>
#include<vcg/complex/trimesh/update/topology.h>
#include<vcg/complex/trimesh/update/bounding.h>
#include<vcg/complex/local_optimization/tri_edge_collapse.h>
#include<vcg/complex/local_optimization.h>


namespace vcg{
namespace tri{




/** 
  This class describe Quadric based collapse operation.

	Requirements:

	Vertex 
	must have:
   incremental mark
   VF topology
  
	must have:
		members
      
      QuadricType Qd();
		 
			ScalarType W() const;
				A per-vertex Weight that can be used in simplification 
				lower weight means that error is lowered, 
				standard: return W==1.0
				
			void Merge(MESH_TYPE::vertex_type const & v);
				Merges the attributes of the current vertex with the ones of v
				(e.g. its weight with the one of the given vertex, the color ect).
				Standard: void function;

      OtherWise the class should be templated with a static helper class that helps to retrieve these functions.
      If the vertex class exposes these functions a default static helper class is provided.

*/
		//**Helper CLASSES**//
		template <class VERTEX_TYPE>
		class QInfoStandard
		{
		public:
			QInfoStandard(){};
      static void Init(){};
      static math::Quadric<double> &Qd(VERTEX_TYPE &v) {return v.Qd();}
      static math::Quadric<double> &Qd(VERTEX_TYPE *v) {return v->Qd();}
      static typename VERTEX_TYPE::ScalarType W(VERTEX_TYPE *v) {return 1.0;};
      static typename VERTEX_TYPE::ScalarType W(VERTEX_TYPE &v) {return 1.0;};
      static void Merge(VERTEX_TYPE & v_dest, VERTEX_TYPE const & v_del){};
		};


class TriEdgeCollapseQuadricParameter
{
public:
	double	QualityThr; // all 
	double	BoundaryWeight;
	double	NormalThrRad;
	double	CosineThr;
	double	QuadricEpsilon;
	double	ScaleFactor;
	bool		UseArea;
	bool		UseVertexWeight;
	bool		NormalCheck;
	bool		QualityCheck;
	bool		OptimalPlacement;
	bool		MemoryLess;
	bool		QualityWeight;
	bool		ScaleIndependent;
	//***********************
	bool    QualityQuadric; // During the initialization manage all the edges as border edges adding a set of additional quadrics that are useful mostly for keeping face aspect ratio good.
	bool		PreserveTopology; 
	bool		PreserveBoundary; 
	bool		MarkComplex;
	bool		FastPreserveBoundary; 
	bool		SafeHeapUpdate;
};


template<class TriMeshType,class MYTYPE, class HelperType = QInfoStandard<typename TriMeshType::VertexType> >
class TriEdgeCollapseQuadric: public TriEdgeCollapse< TriMeshType, MYTYPE> 
{
public:
		typedef typename vcg::tri::TriEdgeCollapse< TriMeshType, MYTYPE > TEC;
		typedef typename TEC::EdgeType EdgeType;
		typedef typename TriEdgeCollapse<TriMeshType, MYTYPE>::HeapType HeapType;
		typedef typename TriEdgeCollapse<TriMeshType, MYTYPE>::HeapElem HeapElem;
		typedef typename TriMeshType::CoordType CoordType;
		typedef typename TriMeshType::ScalarType ScalarType;
		typedef  math::Quadric< double > QuadricType;
		typedef typename TriMeshType::FaceType FaceType;
		typedef typename TriMeshType::VertexType VertexType;
    typedef TriEdgeCollapseQuadricParameter QParameter;
    typedef HelperType QH;

		static QParameter & Params(){
			static QParameter p; 
			return p;
		}
		enum Hint {
			HNHasFFTopology       = 0x0001,  // La mesh arriva con la topologia ff gia'fatta
			HNHasVFTopology       = 0x0002,  // La mesh arriva con la topologia bf gia'fatta
			HNHasBorderFlag       = 0x0004  // La mesh arriva con i flag di bordo gia' settati
		};

		static int & Hnt(){static int hnt; return hnt;}      // the current hints

		static void SetHint(Hint hn)		{	Hnt() |= hn; }
		static void ClearHint(Hint hn)	{	Hnt()&=(~hn);}
		static bool IsSetHint(Hint hn)  { return (Hnt()&hn)!=0; }

		// puntatori ai vertici che sono stati messi non-w per preservare il boundary
		static std::vector<typename TriMeshType::VertexPointer>  & WV(){
      static std::vector<typename TriMeshType::VertexPointer> _WV; return _WV;
    }; 

		inline TriEdgeCollapseQuadric(const EdgeType &p, int i)
			//:TEC(p,i){}
		{
				this->localMark = i;
				this->pos=p;
				this->_priority = ComputePriority();
		}


		inline bool IsFeasible(){
      bool res = ( !Params().PreserveTopology || LinkConditions(this->pos) );
      if(!res) ++( TriEdgeCollapse< TriMeshType,MYTYPE>::FailStat::LinkConditionEdge() );
      return res;
    }

		void Execute(TriMeshType &m)
  {	CoordType newPos;
		if(Params().OptimalPlacement) newPos= static_cast<MYTYPE*>(this)->ComputeMinimal();
    else newPos=this->pos.V(1)->P();
		//this->pos.V(1)->Qd()+=this->pos.V(0)->Qd();
    QH::Qd(this->pos.V(1))+=QH::Qd(this->pos.V(0));
		//int FaceDel=
		DoCollapse(m, this->pos, newPos); // v0 is deleted and v1 take the new position
		//m.fn-=FaceDel;
		//--m.vn;
  }

  
    
    // Final Clean up after the end of the simplification process
    static void Finalize(TriMeshType &m, HeapType& /*h_ret*/)
    {
      // if the mesh was prepared with precomputed borderflags 
      // correctly set them again.
      if(IsSetHint(HNHasBorderFlag) ) 
			  vcg::tri::UpdateFlags<TriMeshType>::FaceBorderFromVF(m);

      // If we had the boundary preservation we should clean up the writable flags
      if(Params().FastPreserveBoundary)
      {
        typename 	TriMeshType::VertexIterator  vi;
    	  for(vi=m.vert.begin();vi!=m.vert.end();++vi) 
          if(!(*vi).IsD()) (*vi).SetW();
      }
    	if(Params().PreserveBoundary)
      {
        typename 	std::vector<typename TriMeshType::VertexPointer>::iterator wvi;
        for(wvi=WV().begin();wvi!=WV().end();++wvi)
          if(!(*wvi)->IsD()) (*wvi)->SetW();
      }
    }

	static void Init(TriMeshType &m,HeapType&h_ret){

	typename 	TriMeshType::VertexIterator  vi;
	typename 	TriMeshType::FaceIterator  pf;

	EdgeType av0,av1,av01;
	Params().CosineThr=cos(Params().NormalThrRad);

	if(!IsSetHint(HNHasVFTopology) ) vcg::tri::UpdateTopology<TriMeshType>::VertexFace(m);

	if(Params().MarkComplex) 		{
		vcg::tri::UpdateTopology<TriMeshType>::FaceFace(m);
		vcg::tri::UpdateFlags<TriMeshType>::FaceBorderFromFF(m);
		vcg::tri::UpdateTopology<TriMeshType>::VertexFace(m);
	} // e' un po' piu' lenta ma marca i vertici complex
	else 
		if(!IsSetHint(HNHasBorderFlag) ) 
			vcg::tri::UpdateFlags<TriMeshType>::FaceBorderFromVF(m);

	if(Params().FastPreserveBoundary)
		{
			for(pf=m.face.begin();pf!=m.face.end();++pf)
			if( !(*pf).IsD() && (*pf).IsW() )
				for(int j=0;j<3;++j)
					if((*pf).IsB(j)) 
					{
						(*pf).V(j)->ClearW();
						(*pf).V1(j)->ClearW();
					}
			}

	if(Params().PreserveBoundary)
		{
      WV().clear();
			for(pf=m.face.begin();pf!=m.face.end();++pf)
			if( !(*pf).IsD() && (*pf).IsW() )
				for(int j=0;j<3;++j)
					if((*pf).IsB(j)) 
					{
						if((*pf).V(j)->IsW())  {(*pf).V(j)->ClearW(); WV().push_back((*pf).V(j));}
						if((*pf).V1(j)->IsW()) {(*pf).V1(j)->ClearW();WV().push_back((*pf).V1(j));}
					}
			}

		InitQuadric(m);

	// Initialize the heap with all the possible collapses 
		if(IsSymmetric()) 
    { // if the collapse is symmetric (e.g. u->v == v->u)
		  for(vi=m.vert.begin();vi!=m.vert.end();++vi) 
				if(!(*vi).IsD() && (*vi).IsRW())
						{
								vcg::face::VFIterator<FaceType> x;
								for( x.F() = (*vi).VFp(), x.I() = (*vi).VFi(); x.F()!=0; ++ x){
									x.V1()->ClearV();
									x.V2()->ClearV();
							  }
								for( x.F() = (*vi).VFp(), x.I() = (*vi).VFi(); x.F()!=0; ++x )
                {
									assert(x.F()->V(x.I())==&(*vi));
									if((x.V0()<x.V1()) && x.V1()->IsRW() && !x.V1()->IsV()){
												x.V1()->SetV();
												h_ret.push_back(HeapElem(new MYTYPE(EdgeType(x.V0(),x.V1()),TriEdgeCollapse< TriMeshType,MYTYPE>::GlobalMark() )));
												}
									if((x.V0()<x.V2()) && x.V2()->IsRW()&& !x.V2()->IsV()){
												x.V2()->SetV();
												h_ret.push_back(HeapElem(new MYTYPE(EdgeType(x.V0(),x.V2()),TriEdgeCollapse< TriMeshType,MYTYPE>::GlobalMark() )));
											}
								}
		        }	
	  }	
		else 
    { // if the collapse is A-symmetric (e.g. u->v != v->u) 
			for(vi=m.vert.begin();vi!=m.vert.end();++vi)
				if(!(*vi).IsD() && (*vi).IsRW())
					{
						vcg::face::VFIterator<FaceType> x;
						UnMarkAll(m);
						for( x.F() = (*vi).VFp(), x.I() = (*vi).VFi(); x.F()!=0; ++ x)
						{
							assert(x.F()->V(x.I())==&(*vi));
							if(x.V()->IsRW() && x.V1()->IsRW() && !IsMarked(m,x.F()->V1(x.I()))){
										h_ret.push_back( HeapElem( new MYTYPE( EdgeType (x.V(),x.V1()),TriEdgeCollapse< TriMeshType,MYTYPE>::GlobalMark())));
										}
							if(x.V()->IsRW() && x.V2()->IsRW() && !IsMarked(m,x.F()->V2(x.I()))){
										h_ret.push_back( HeapElem( new MYTYPE( EdgeType (x.V(),x.V2()),TriEdgeCollapse< TriMeshType,MYTYPE>::GlobalMark())));
									}
						}
					}	
	  }
}
	static float HeapSimplexRatio() {return IsSymmetric()?5.0f:9.0f;}
	static bool IsSymmetric() {return Params().OptimalPlacement;} 
	static bool IsVertexStable() {return !Params().OptimalPlacement;}
	static void SetDefaultParams(){
		Params().UseArea=true;
		Params().UseVertexWeight=false;
		Params().NormalCheck=false;
		Params().NormalThrRad=M_PI/2;
		Params().QualityCheck=true;
		Params().QualityThr=.1;
		Params().BoundaryWeight=.5;
		Params().QualityQuadric=false;
		Params().OptimalPlacement=true;
		Params().ScaleIndependent=true;
		Params().QualityWeight=false;
		Params().QuadricEpsilon =1e-15;
		Params().ScaleFactor=1.0;

		Params().PreserveTopology = false;
	}
	
///*
//	Funzione principale di valutazione dell'errore del collasso.
//	In pratica simula il collasso vero e proprio.
//
//	Da ottimizzare il ciclo sulle normali (deve sparire on e si deve usare per face normals)
//*/
	ScalarType ComputePriority()  {
		ScalarType error;
		typename vcg::face::VFIterator<FaceType> x;
		std::vector<CoordType> on; // original normals
		typename TriMeshType::VertexType * v[2];
		v[0] = this->pos.V(0);
		v[1] = this->pos.V(1);

		if(Params().NormalCheck){ // Compute maximal normal variation 
			// store the old normals for non-collapsed face in v0
			for(x.F() = v[0]->VFp(), x.I() = v[0]->VFi(); x.F()!=0; ++x )	 // for all faces in v0		
				if(x.F()->V(0)!=v[1] && x.F()->V(1)!=v[1] && x.F()->V(2)!=v[1] ) // skip faces with v1
					on.push_back(NormalizedNormal(*x.F()));
				// store the old normals for non-collapsed face in v1
				for(x.F() = v[1]->VFp(), x.I() = v[1]->VFi(); x.F()!=0; ++x )	 // for all faces in v1	
					if(x.F()->V(0)!=v[0] && x.F()->V(1)!=v[0] && x.F()->V(2)!=v[0] ) // skip faces with v0
						on.push_back(NormalizedNormal(*x.F()));
			}

		//// Move the two vertexe  into new position (storing the old ones)
		CoordType OldPos0=v[0]->P();
		CoordType OldPos1=v[1]->P();
		if(Params().OptimalPlacement)	 { v[0]->P() = ComputeMinimal(); v[1]->P()=v[0]->P();}
				else  v[0]->P() = v[1]->P();

		//// Rescan faces and compute quality and difference between normals
		int i;
		double ndiff,MinCos  = 1e100; // minimo coseno di variazione di una normale della faccia 
																	// (e.g. max angle) Mincos varia da 1 (normali coincidenti) a
																	// -1 (normali opposte);
		double qt,   MinQual = 1e100;
		CoordType nn;
		for(x.F() = v[0]->VFp(), x.I() = v[0]->VFi(),i=0; x.F()!=0; ++x )	// for all faces in v0		
			if(x.F()->V(0)!=v[1] && x.F()->V(1)!=v[1] && x.F()->V(2)!=v[1] )		// skip faces with v1
			{
				if(Params().NormalCheck){
					nn=NormalizedNormal(*x.F());
					ndiff=nn.dot(on[i++]);
					if(ndiff<MinCos) MinCos=ndiff;
				}
				if(Params().QualityCheck){
					qt= QualityFace(*x.F());
					if(qt<MinQual) MinQual=qt;
				}
			}
		for(x.F() = v[1]->VFp(), x.I() = v[1]->VFi(),i=0; x.F()!=0; ++x )		// for all faces in v1	
			if(x.F()->V(0)!=v[0] && x.F()->V(1)!=v[0] && x.F()->V(2)!=v[0] )			// skip faces with v0
			{
				if(Params().NormalCheck){
					nn=NormalizedNormal(*x.F());
					ndiff=nn.dot(on[i++]);
					if(ndiff<MinCos) MinCos=ndiff;
				}
				if(Params().QualityCheck){
					qt= QualityFace(*x.F());
					if(qt<MinQual) MinQual=qt;
				}
			}

    QuadricType qq=QH::Qd(v[0]);
    qq+=QH::Qd(v[1]);
		Point3d tpd=Point3d::Construct(v[1]->P());
    double QuadErr = Params().ScaleFactor*qq.Apply(tpd); 

		// All collapses involving triangles with quality larger than <QualityThr> has no penalty;
		if(MinQual>Params().QualityThr) MinQual=Params().QualityThr;  
				
		if(Params().NormalCheck){
			// All collapses where the normal vary less  than <NormalThr> (e.g. more than CosineThr)
			// have no penalty
			if(MinCos>Params().CosineThr) MinCos=Params().CosineThr;
			MinCos=(MinCos+1)/2.0; // Now it is in the range 0..1 with 0 very dangerous!
		}
		
		if(QuadErr<Params().QuadricEpsilon) QuadErr=Params().QuadricEpsilon;
		
		if( Params().UseVertexWeight ) QuadErr *= (QH::W(v[1])+QH::W(v[0]))/2;
		
		if(!Params().QualityCheck && !Params().NormalCheck) error = (ScalarType)(QuadErr);
		if( Params().QualityCheck && !Params().NormalCheck) error = (ScalarType)(QuadErr / MinQual);
		if(!Params().QualityCheck &&  Params().NormalCheck) error = (ScalarType)(QuadErr / MinCos);
		if( Params().QualityCheck &&  Params().NormalCheck) error = (ScalarType)(QuadErr / (MinQual*MinCos));
																							           
		//Rrestore old position of v0 and v1
		v[0]->P()=OldPos0;
		v[1]->P()=OldPos1;
		this->_priority = error;
		return this->_priority;
	}

//	
//static double MaxError() {return 1e100;}
//
  inline  void UpdateHeap(HeapType & h_ret)
  {
		this->GlobalMark()++;
		VertexType *v[2];
		v[0]= this->pos.V(0);
        v[1]= this->pos.V(1);	
		v[1]->IMark() = this->GlobalMark();

		// First loop around the remaining vertex to unmark visited flags
    vcg::face::VFIterator<FaceType> vfi(v[1]);	
		while (!vfi.End()){
			vfi.V1()->ClearV();
			vfi.V2()->ClearV();
			++vfi;
		}

    // Second Loop 
		vfi = face::VFIterator<FaceType>(v[1]);	
		while (!vfi.End())
    {
			assert(!vfi.F()->IsD());
      for (int j=0;j<3;j++)
			{
				if( !(vfi.V1()->IsV()) && vfi.V1()->IsRW())
				{
				  vfi.V1()->SetV();
          h_ret.push_back(HeapElem(new MYTYPE(EdgeType(vfi.V0(),vfi.V1()), this->GlobalMark())));
				  std::push_heap(h_ret.begin(),h_ret.end());
				  if(!IsSymmetric()){				
					  h_ret.push_back(HeapElem(new MYTYPE(EdgeType(vfi.V1(),vfi.V0()), this->GlobalMark())));
					  std::push_heap(h_ret.begin(),h_ret.end());
				  }
        }
				if(  !(vfi.V2()->IsV()) && vfi.V2()->IsRW())
				{
					vfi.V2()->SetV();
				  h_ret.push_back(HeapElem(new MYTYPE(EdgeType(vfi.V0(),vfi.V2()),this->GlobalMark())));
				  std::push_heap(h_ret.begin(),h_ret.end());
				  if(!IsSymmetric()){				
					  h_ret.push_back(HeapElem(new MYTYPE(EdgeType(vfi.V2(),vfi.V0()), this->GlobalMark())));
					  std::push_heap(h_ret.begin(),h_ret.end());
				  }
        }
        if(Params().SafeHeapUpdate && vfi.V1()->IsRW() && vfi.V2()->IsRW() )
        {
          h_ret.push_back(HeapElem(new MYTYPE(EdgeType(vfi.V1(),vfi.V2()),this->GlobalMark())));
				  std::push_heap(h_ret.begin(),h_ret.end());
				  if(!IsSymmetric()){				
					  h_ret.push_back(HeapElem(new MYTYPE(EdgeType(vfi.V2(),vfi.V1()), this->GlobalMark())));
					  std::push_heap(h_ret.begin(),h_ret.end());
				  }
        }
			}
      ++vfi;
    }

  }

static void InitQuadric(TriMeshType &m)
{
	typename TriMeshType::FaceIterator pf;
	typename TriMeshType::VertexIterator pv;
	int j;
  QH::Init();
	//	m.ClearFlags();
	for(pv=m.vert.begin();pv!=m.vert.end();++pv)		// Azzero le quadriche
		if( ! (*pv).IsD() && (*pv).IsW()) 	
      QH::Qd(*pv).SetZero();

		
	for(pf=m.face.begin();pf!=m.face.end();++pf)
		if( !(*pf).IsD() && (*pf).IsR() )
			if((*pf).V(0)->IsR() &&(*pf).V(1)->IsR() &&(*pf).V(2)->IsR())
					{
						QuadricType q;
						Plane3<ScalarType,false> p;
							// Calcolo piano
						p.SetDirection( ( (*pf).V(1)->cP() - (*pf).V(0)->cP() ) ^  ( (*pf).V(2)->cP() - (*pf).V(0)->cP() ));
							// Se normalizzo non dipende dall'area

						if(!Params().UseArea)	
							p.Normalize();

						p.SetOffset( p.Direction().dot((*pf).V(0)->cP()));

							// Calcolo quadrica	delle facce
						q.ByPlane(p);

						for(j=0;j<3;++j)
              if( (*pf).V(j)->IsW() )	
								{
									if(Params().QualityWeight) 
												q*=(*pf).V(j)->Q();
									QH::Qd((*pf).V(j)) += q;				// Sommo la quadrica ai vertici
								}
						
						for(j=0;j<3;++j)
							if( (*pf).IsB(j) || Params().QualityQuadric )				// Bordo!
							{
								Plane3<ScalarType,false> pb;						// Piano di bordo

								// Calcolo la normale al piano di bordo e la sua distanza
								// Nota che la lunghezza dell'edge DEVE essere Normalizzata 
								// poiche' la pesatura in funzione dell'area e'gia fatta in p.Direction() 
								// Senza la normalize il bordo e' pesato in funzione della grandezza della mesh (mesh grandi non decimano sul bordo)
								pb.SetDirection(p.Direction() ^ ( (*pf).V1(j)->cP() - (*pf).V(j)->cP() ).normalized());
								if(  (*pf).IsB(j) ) pb.SetDirection(pb.Direction()* (ScalarType)Params().BoundaryWeight);        // amplify border planes
								               else pb.SetDirection(pb.Direction()* (ScalarType)(Params().BoundaryWeight/100.0)); // and consider much less quadric for quality
								pb.SetOffset(pb.Direction().dot((*pf).V(j)->cP()));
								q.ByPlane(pb);

								if( (*pf).V (j)->IsW() )	QH::Qd((*pf).V (j)) += q;			// Sommo le quadriche
								if( (*pf).V1(j)->IsW() )	QH::Qd((*pf).V1(j)) += q;
							}
					}
   
	if(Params().ScaleIndependent)
	{
		vcg::tri::UpdateBounding<TriMeshType>::Box(m);
		//Make all quadric independent from mesh size
		Params().ScaleFactor = 1e8*pow(1.0/m.bbox.Diag(),6); // scaling factor
		//Params().ScaleFactor *=Params().ScaleFactor ;
		//Params().ScaleFactor *=Params().ScaleFactor ;
		//printf("Scale factor =%f\n",Params().ScaleFactor );
		//printf("bb (%5.2f %5.2f %5.2f)-(%5.2f %5.2f %5.2f) Diag %f\n",m.bbox.min[0],m.bbox.min[1],m.bbox.min[2],m.bbox.max[0],m.bbox.max[1],m.bbox.max[2],m.bbox.Diag());
	}				
}



//
//
//
//
//
//
//static void InitMesh(MESH_TYPE &m){
//	Params().CosineThr=cos(Params().NormalThr);
//  InitQuadric(m);
//	//m.Topology();
//	//OldInitQuadric(m,UseArea);
//	}
//
 CoordType ComputeMinimal()
{	
		typename TriMeshType::VertexType * v[2];
		v[0] = this->pos.V(0);
		v[1] = this->pos.V(1);
		QuadricType q=QH::Qd(v[0]);
		q+=QH::Qd(v[1]);
		
    Point3<QuadricType::ScalarType> x;
		
    bool rt=q.Minimum(x);
		if(!rt) { // if the computation of the minimum fails we choose between the two edge points and the middle one.
			Point3<QuadricType::ScalarType> x0=Point3d::Construct(v[0]->P());
		  Point3<QuadricType::ScalarType> x1=Point3d::Construct(v[1]->P());
      x.Import((v[0]->P()+v[1]->P())/2);
			double qvx=q.Apply(x);
			double qv0=q.Apply(x0);
			double qv1=q.Apply(x1);
      if(qv0<qvx) x=x0;
			if(qv1<qvx && qv1<qv0) x=x1;
		}
		
    return CoordType::Construct(x);
}
//
//

};
		} // namespace tri
	} // namespace vcg
#endif
