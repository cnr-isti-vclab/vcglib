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
#include<vcg/complex/trimesh/update/bounding.h>
#include<vcg/complex/local_optimization/tri_edge_collapse.h>
#include<vcg/complex/local_optimization.h>


namespace vcg{
namespace tri{

class QCollapseParameter
{
public:
	double	QualityThr; // all 
	double	BoundaryWeight;
	double	NormalThr;
	double	CosineThr;
	double	QuadricEpsilon;
	double	ScaleFactor;
	bool		UseArea;
	bool		UseVertexWeight;
	bool		NormalCheck;
	bool		QualityCheck;
	bool		OptimalPlacement;
	bool		MemoryLess;
	bool		ComplexCheck;
	bool		ScaleIndependent;
	//***********************
	bool		PreserveTopology; 
	bool		PreserveBoundary; 
	bool		MarkComplex;
	bool		FastPreserveBoundary; 
	bool		SafeHeapUpdate;
};


/** 
  This class describe Quadric based collapse operation.

	Requirements:

	Vertex 
	must have incremental mark
	must have:
		field QuadricType Q;
		member 
			ScalarType W() const;
				A per-vertex Weight that can be used in simplification 
				lower weight means that error is lowered, 
				standard: return W==1.0
				
			void Merge(MESH_TYPE::vertex_type const & v);
				Merges the attributes of the current vertex with the ones of v
				(e.g. its weight with the one of the given vertex, the color ect).
				Standard: void function;

	Faces devono avere Shared Adjacency 
		durante la init serve FF per le quadriche di bordo
		durante la semplificazione si usa VF

*/

template<class TriMeshType,class MYTYPE>
class TriEdgeCollapseQuadric: public TriEdgeCollapse< TriMeshType,MYTYPE> 
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

		static QCollapseParameter & Params(){static QCollapseParameter p; return p;}
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
		static std::vector<typename TriMeshType::VertexPointer>  & WV(){static std::vector<typename TriMeshType::VertexPointer> _WV; return _WV;}; 

		inline TriEdgeCollapseQuadric(const EdgeType &p, int i)
			//:TEC(p,i){}
		{
				localMark = i;
				pos=p;
				_priority = ComputePriority();
		}


		inline bool IsFeasible(){
      bool res = ( !Params().PreserveTopology || LinkConditions(pos) );
      if(!res) ++( TriEdgeCollapse< TriMeshType,MYTYPE>::FailStat::LinkConditionEdge() );
      return res;
    }

		void Execute(TriMeshType &m)
  {	
		CoordType newPos = ComputeMinimal();
		pos.V(1)->q+=pos.V(0)->q;
		int FaceDel=DoCollapse(pos, newPos);
		m.fn-=FaceDel;
		--m.vn;
  }

	static void Init(TriMeshType &m,HeapType&h_ret){

	typename 	TriMeshType::VertexIterator  vi;
	typename 	TriMeshType::FaceIterator  pf;

	EdgeType av0,av1,av01;
	Params().CosineThr=cos(Params().NormalThr);

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
		if(IsSymmetric()) { // if the collapse is symmetric (e.g. u->v == v->u)
		for(vi=m.vert.begin();vi!=m.vert.end();++vi) 
				if((*vi).IsRW())
						{
								vcg::face::VFIterator<FaceType> x;
								for( x.F() = (*vi).VFp(), x.I() = (*vi).VFi(); x.F()!=0; ++ x){
									x.F()->V1(x.I())->ClearV();
									x.F()->V2(x.I())->ClearV();
							}
								for( x.F() = (*vi).VFp(), x.I() = (*vi).VFi(); x.F()!=0; ++x ){
									assert(x.F()->V(x.I())==&(*vi));
									if((x.F()->V(x.I())<x.F()->V1(x.I())) && x.F()->V1(x.I())->IsRW() && !x.F()->V1(x.I())->IsV()){
												x.F()->V1(x.I())->SetV();

												h_ret.push_back(HeapElem(new MYTYPE(EdgeType(x.F()->V(x.I()),x.F()->V1(x.I())),GlobalMark())));
												}
									if((x.F()->V(x.I())<x.F()->V2(x.I())) && x.F()->V2(x.I())->IsRW()&& !x.F()->V2(x.I())->IsV()){
												x.F()->V2(x.I())->SetV();
												h_ret.push_back(HeapElem(new MYTYPE(EdgeType(x.F()->V(x.I()),x.F()->V2(x.I())),GlobalMark() )));
											}
								}
		}	
	}	
		else { // if the collapse is A-symmetric (e.g. u->v != v->u) 
			for(vi=m.vert.begin();vi!=m.vert.end();++vi)
		{
			vcg::face::VFIterator<FaceType> x;
			m.UnMarkAll();
			for( x.F() = (*vi).VFp(), x.I() = (*vi).VFi(); x.F()!=0; ++ x){
				assert(x.F()->V(x.I())==&(*vi));
				if(x.F()->V(x.I())->IsRW() && x.F()->V1(x.I())->IsRW() && !m.IsMarked(x.F()->V1(x.I()))){
							m.Mark( x.F()->V1(x.I()) );
							h_ret.push_back( HeapElem( new MYTYPE( EdgeType (x.F()->V(x.I()),x.F()->V1(x.I())), m.imark)));
							}
				if(x.F()->V(x.I())->IsRW() && x.F()->V2(x.I())->IsRW()&& !m.IsMarked(x.F()->V2(x.I()))){
							m.Mark( x.F()->V2(x.I()) );
							h_ret.push_back( HeapElem( new MYTYPE( EdgeType (x.F()->V(x.I()),x.F()->V2(x.I())), m.imark)));
						}
			}
		}	
	}
	make_heap(h_ret.begin(),h_ret.end());
}

	static bool IsSymmetric() {return Params().OptimalPlacement;} 
	static bool IsVertexStable() {return !Params().OptimalPlacement;}
	static void SetDefaultParams(){
		Params().UseArea=true;
		Params().UseVertexWeight=false;
		Params().NormalCheck=false;
		Params().NormalThr=M_PI/2;
		Params().QualityCheck=true;
		Params().QualityThr=.1;
		Params().BoundaryWeight=.5;
		Params().OptimalPlacement=true;
		Params().ScaleIndependent=true;
		Params().ComplexCheck=false;
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
		v[0] = pos.V(0);
		v[1] = pos.V(1);

		if(Params().NormalCheck){ // Compute maximal normal variation 
			// store the old normals for non-collapsed face in v0
			for(x.F() = v[0]->VFp(), x.I() = v[0]->VFi(); x.F()!=0; ++x )	 // for all faces in v0		
				if(x.F()->V(0)!=v[1] && x.F()->V(1)!=v[1] && x.F()->V(2)!=v[1] ) // skip faces with v1
					on.push_back(x.F()->NormalizedNormal());
				// store the old normals for non-collapsed face in v1
				for(x.F() = v[1]->VFp(), x.I() = v[1]->VFi(); x.F()!=0; ++x )	 // for all faces in v1	
					if(x.F()->V(0)!=v[0] && x.F()->V(1)!=v[0] && x.F()->V(2)!=v[0] ) // skip faces with v0
						on.push_back(x.F()->NormalizedNormal());
			}

		//// Move the two vertexe  into new position (storing the old ones)
		CoordType OldPos0=v[0]->P();
		CoordType OldPos1=v[1]->P();
		if(Params().OptimalPlacement)	 
			{ v[0]->P() = ComputeMinimal(); v[1]->P()=v[0]->P();}
				else  
			v[0]->P() = v[1]->P();

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
					nn=x.F()->NormalizedNormal();
					ndiff=nn*on[i++];
					if(ndiff<MinCos) MinCos=ndiff;
				}
				if(Params().QualityCheck){
					qt= x.F()->QualityFace();
					if(qt<MinQual) MinQual=qt;
				}
			}
		for(x.F() = v[1]->VFp(), x.I() = v[1]->VFi(),i=0; x.F()!=0; ++x )		// for all faces in v1	
			if(x.F()->V(0)!=v[0] && x.F()->V(1)!=v[0] && x.F()->V(2)!=v[0] )			// skip faces with v0
			{
				if(Params().NormalCheck){
					nn=x.F()->NormalizedNormal();
					ndiff=nn*on[i++];
					if(ndiff<MinCos) MinCos=ndiff;
				}
				if(Params().QualityCheck){
					qt= x.F()->QualityFace();
					if(qt<MinQual) MinQual=qt;
				}
			}

		QuadricType qq=v[0]->q;
		qq+=v[1]->q;
		double QuadErr = Params().ScaleFactor*qq.Apply(v[1]->P()); 

		// All collapses involving triangles with quality larger than <QualityThr> has no penalty;
		if(MinQual>Params().QualityThr) MinQual=Params().QualityThr;  
				
		if(Params().NormalCheck){
			// All collapses where the normal vary less  than <NormalThr> (e.g. more than CosineThr)
			// have no penalty
			if(MinCos>Params().CosineThr) MinCos=Params().CosineThr;
			MinCos=(MinCos+1)/2.0; // Now it is in the range 0..1 with 0 very dangerous!
		}
		
		if(QuadErr<Params().QuadricEpsilon) QuadErr=Params().QuadricEpsilon;
		
		if( Params().UseVertexWeight ) QuadErr *= (v[1]->W()+v[0]->W())/2;
		
		if(!Params().QualityCheck && !Params().NormalCheck) error = QuadErr;
		if( Params().QualityCheck && !Params().NormalCheck) error = QuadErr / MinQual;
		if(!Params().QualityCheck &&  Params().NormalCheck) error = QuadErr / MinCos;
		if( Params().QualityCheck &&  Params().NormalCheck) error = QuadErr / (MinQual*MinCos);
																							           
		//Rrestore old position of v0 and v1
		v[0]->P()=OldPos0;
		v[1]->P()=OldPos1;
		_priority = -error;
		return _priority;
	}

//	
//static double MaxError() {return 1e100;}
//

static void InitQuadric(TriMeshType &m)
{
	typename TriMeshType::FaceIterator pf;
	typename TriMeshType::VertexIterator pv;
	int j;
	//	m.ClearFlags();
	for(pv=m.vert.begin();pv!=m.vert.end();++pv)		// Azzero le quadriche
		if( ! (*pv).IsD() && (*pv).IsW()) 	
			(*pv).q.Zero();

		
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

						p.SetOffset( p.Direction() * (*pf).V(0)->cP());

							// Calcolo quadrica	delle facce
						q.ByPlane(p);

						for(j=0;j<3;++j)
							if( (*pf).V(j)->IsW() )	(*pf).V(j)->q += q;				// Sommo la quadrica ai vertici
						
						for(j=0;j<3;++j)
							if( (*pf).IsB(j))				// Bordo!
							{
								Plane3<ScalarType,false> pb;						// Piano di bordo

								// Calcolo la normale al piano di bordo e la sua distanza
								// Nota che la lunghezza dell'edge DEVE essere Normalizzata 
								// poiche' la pesatura in funzione dell'area e'gia fatta in p.Direction() 
								// Senza la normalize il bordo e' pesato in funzione della grandezza della mesh (mesh grandi non decimano sul bordo)
								pb.SetDirection(p.Direction() ^ ( (*pf).V1(j)->cP() - (*pf).V(j)->cP() ).Normalize());
								pb.SetDirection(pb.Direction()* Params().BoundaryWeight);  // amplify border planes
								pb.SetOffset(pb.Direction() * (*pf).V(j)->cP());
								q.ByPlane(pb);

								if( (*pf).V (j)->IsW() )	(*pf).V (j)->q += q;			// Sommo le quadriche
								if( (*pf).V1(j)->IsW() )	(*pf).V1(j)->q += q;
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
		

	if(Params().ComplexCheck)
	{
		// secondo loop per diminuire quadriche complex (se non c'erano i complex si poteva fare in un giro solo)
		//for(pf=m.face.begin();pf!=m.face.end();++pf)
		//if( !(*pf).IsD() && (*pf).IsR() )
		//	if((*pf).V(0)->IsR() &&(*pf).V(1)->IsR() &&(*pf).V(2)->IsR())
		//			{
		//				for(j=0;j<3;++j)
		//				if((*pf).IsCF(j))				// Complex!
		//				{
		//					if( (*pf).V (j)->IsW() )	(*pf).V (j)->q *= 0.01;			// Scalo le quadriche
		//					if( (*pf).V1(j)->IsW() )	(*pf).V1(j)->q *= 0.01;
		//				}
		//		}
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
		v[0] = pos.V(0);
		v[1] = pos.V(1);
		QuadricType q=v[0]->q;
		q+=v[1]->q;
		
		CoordType x;
		bool rt=q.Minimum(x);
		if(!rt) {
			x=(v[0]->P()+v[1]->P())/2;
			double qvx=q.Apply(x);
			double qv0=q.Apply(v[0]->P());
			double qv1=q.Apply(v[1]->P());
			if(qv0<qvx) x=v[0]->P();
			if(qv1<qvx && qv1<qv0) x=v[1]->P();
		}
		
		return x;
}
//
//

};
		} // namespace tri
	} // namespace vcg
#endif
