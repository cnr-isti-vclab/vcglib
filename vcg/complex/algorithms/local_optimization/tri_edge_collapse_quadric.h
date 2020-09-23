/****************************************************************************
* VCGLib                                                            o o     *
* Visual and Computer Graphics Library                            o     o   *
*                                                                _   O  _   *
* Copyright(C) 2004-2016                                           \/)\/    *
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

#ifndef __VCG_TRIMESHCOLLAPSE_QUADRIC__
#define __VCG_TRIMESHCOLLAPSE_QUADRIC__

#include<vcg/math/quadric.h>
#include<vcg/complex/algorithms/update/bounding.h>
#include<vcg/complex/algorithms/local_optimization/tri_edge_collapse.h>
#include<vcg/complex/algorithms/local_optimization.h>
#include<vcg/complex/algorithms/stat.h>


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
      QInfoStandard(){}
      static void Init(){}
      static math::Quadric<double> &Qd(VERTEX_TYPE &v) {return v.Qd();}
      static math::Quadric<double> &Qd(VERTEX_TYPE *v) {return v->Qd();}
      static typename VERTEX_TYPE::ScalarType W(VERTEX_TYPE * /*v*/) {return 1.0;}
      static typename VERTEX_TYPE::ScalarType W(VERTEX_TYPE &/*v*/) {return 1.0;}
      static void Merge(VERTEX_TYPE & /*v_dest*/, VERTEX_TYPE const & /*v_del*/){}
        };


class TriEdgeCollapseQuadricParameter : public BaseParameterClass
{
public:
  double    BoundaryQuadricWeight = 0.5;
  bool      FastPreserveBoundary  = false;
  bool      AreaCheck           = false;
  bool      HardQualityCheck = false;
  double    HardQualityThr = 0.1;
  bool      HardNormalCheck =  false;
  bool      NormalCheck           = false;
  double    NormalThrRad          = M_PI/2.0;
  double    CosineThr             = 0 ; // ~ cos(pi/2) 
  bool      OptimalPlacement =true;
  bool      SVDPlacement = false;
  bool      PreserveTopology =false;
  bool      PreserveBoundary = false;
  double    QuadricEpsilon = 1e-15;
  bool      QualityCheck =true;
  double    QualityThr =.3;     // Collapsed that generate faces with quality LOWER than this value are penalized. So higher the value -> better the quality of the accepted triangles
  bool      QualityQuadric =false; // During the initialization manage all the edges as border edges adding a set of additional quadrics that are useful mostly for keeping face aspect ratio good.
  double    QualityQuadricWeight = 0.001f; // During the initialization manage all the edges as border edges adding a set of additional quadrics that are useful mostly for keeping face aspect ratio good.
  bool      QualityWeight=false;
  double    QualityWeightFactor=100.0;
  double    ScaleFactor=1.0;
  bool      ScaleIndependent=true;
  bool      UseArea =true;
  bool      UseVertexWeight=false;  

  TriEdgeCollapseQuadricParameter() {}
};


template<class TriMeshType, class VertexPair, class MYTYPE, class HelperType = QInfoStandard<typename TriMeshType::VertexType> >
class TriEdgeCollapseQuadric: public TriEdgeCollapse< TriMeshType, VertexPair, MYTYPE>
{
public:
  typedef typename vcg::tri::TriEdgeCollapse< TriMeshType, VertexPair, MYTYPE > TEC;
  typedef typename TriEdgeCollapse<TriMeshType, VertexPair, MYTYPE>::HeapType HeapType;
  typedef typename TriEdgeCollapse<TriMeshType, VertexPair, MYTYPE>::HeapElem HeapElem;
  typedef typename TriMeshType::CoordType CoordType;
  typedef typename TriMeshType::ScalarType ScalarType;
  typedef typename TriMeshType::FaceType FaceType;
  typedef typename TriMeshType::VertexType VertexType;
  typedef typename TriMeshType::VertexIterator VertexIterator;
  typedef typename TriMeshType::FaceIterator FaceIterator;
  typedef typename vcg::face::VFIterator<FaceType> VFIterator;
  typedef  math::Quadric< double > QuadricType;
  typedef TriEdgeCollapseQuadricParameter QParameter;
  typedef HelperType QH;
  
  CoordType optimalPos;  // Local storage of the once computed optimal position of the collapse.
  
  // Pointer to the vector that store the Write flags. Used to preserve them if you ask to preserve for the boundaries.
  static std::vector<typename TriMeshType::VertexPointer>  & WV(){
    static std::vector<typename TriMeshType::VertexPointer> _WV; return _WV;
  }
  
  inline TriEdgeCollapseQuadric(){}
  
  inline TriEdgeCollapseQuadric(const VertexPair &p, int i, BaseParameterClass *pp)
  {
    this->localMark = i;
    this->pos=p;
    this->_priority = ComputePriority(pp);
  }
  

  inline bool IsFeasible(BaseParameterClass *_pp){
    QParameter *pp=(QParameter *)_pp;
    if(!pp->PreserveTopology) return true;
    
    bool res = ( EdgeCollapser<TriMeshType, VertexPair>::LinkConditions(this->pos) );
    if(!res) ++( TEC::FailStat::LinkConditionEdge() );
    return res;
  }
  
  void ComputePosition(BaseParameterClass *_pp)
  {
    QParameter *pp=(QParameter *)_pp;
    CoordType newPos = (this->pos.V(0)->P()+this->pos.V(1)->P())/2.0;
    if(pp->OptimalPlacement==false)
      newPos=this->pos.V(1)->P();      
    else 
    {
      if((QH::Qd(this->pos.V(0)).Apply(newPos) + QH::Qd(this->pos.V(1)).Apply(newPos)) > 2.0*pp->QuadricEpsilon)              
      {
        QuadricType q=QH::Qd(this->pos.V(0));
        q+=QH::Qd(this->pos.V(1));
        
        Point3<QuadricType::ScalarType> x;
        if(pp->SVDPlacement)
          q.MinimumClosestToPoint(x,Point3d::Construct(newPos));
        else 
          q.Minimum(x);
        newPos = CoordType::Construct(x);  
      }      
    }
    this->optimalPos = newPos;
  }
  
  void Execute(TriMeshType &m, BaseParameterClass * /*_pp*/)
  {
    CoordType newPos = this->optimalPos;
    QH::Qd(this->pos.V(1))+=QH::Qd(this->pos.V(0)); // v0 is deleted and v1 take the new position
    EdgeCollapser<TriMeshType,VertexPair>::Do(m, this->pos, newPos); 
  }
  
  // Final Clean up after the end of the simplification process
  static void Finalize(TriMeshType &m, HeapType& /*h_ret*/, BaseParameterClass *_pp)
  {
    QParameter *pp=(QParameter *)_pp;
    
    // If we had the boundary preservation we should clean up the writable flags
    if(pp->FastPreserveBoundary)
    {
      typename 	TriMeshType::VertexIterator  vi;
      for(vi=m.vert.begin();vi!=m.vert.end();++vi)
        if(!(*vi).IsD()) (*vi).SetW();
    }
    if(pp->PreserveBoundary)
    {
      typename 	std::vector<typename TriMeshType::VertexPointer>::iterator wvi;
      for(wvi=WV().begin();wvi!=WV().end();++wvi)
        if(!(*wvi)->IsD()) (*wvi)->SetW();
    }
  }
  
  static void Init(TriMeshType &m, HeapType &h_ret, BaseParameterClass *_pp)
  {
    QParameter *pp=(QParameter *)_pp;    
    pp->CosineThr=cos(pp->NormalThrRad);
    h_ret.clear();
    vcg::tri::UpdateTopology<TriMeshType>::VertexFace(m);
    vcg::tri::UpdateFlags<TriMeshType>::FaceBorderFromVF(m);
    
    if(pp->FastPreserveBoundary)
    {
      for(auto pf=m.face.begin();pf!=m.face.end();++pf)
        if( !(*pf).IsD() && (*pf).IsW() )
          for(int j=0;j<3;++j)
            if((*pf).IsB(j))
            {
              (*pf).V(j)->ClearW();
              (*pf).V1(j)->ClearW();
            }
    }
    
    if(pp->PreserveBoundary)
    {
      WV().clear();
      for(auto pf=m.face.begin();pf!=m.face.end();++pf)
        if( !(*pf).IsD() && (*pf).IsW() )
          for(int j=0;j<3;++j)
            if((*pf).IsB(j))
            {
              if((*pf).V(j)->IsW())  {(*pf).V(j)->ClearW(); WV().push_back((*pf).V(j));}
              if((*pf).V1(j)->IsW()) {(*pf).V1(j)->ClearW();WV().push_back((*pf).V1(j));}
            }
    }
    
    InitQuadric(m,pp);
    
    // Initialize the heap with all the possible collapses
    if(IsSymmetric(pp))
    { // if the collapse is symmetric (e.g. u->v == v->u)
      for(auto vi=m.vert.begin();vi!=m.vert.end();++vi)
        if(!(*vi).IsD() && (*vi).IsRW())
        {
          vcg::face::VFIterator<FaceType> x;
          for( x.F() = (*vi).VFp(), x.I() = (*vi).VFi(); x.F()!=0; ++ x){
            x.V1()->ClearV();
            x.V2()->ClearV();
          }
          for( x.F() = (*vi).VFp(), x.I() = (*vi).VFi(); x.F()!=0; ++x )
          {
            if((x.V0()<x.V1()) && x.V1()->IsRW() && !x.V1()->IsV()){
              x.V1()->SetV();
              h_ret.push_back(HeapElem(new MYTYPE(VertexPair(x.V0(),x.V1()),TriEdgeCollapse< TriMeshType,VertexPair,MYTYPE>::GlobalMark(),_pp )));
            }
            if((x.V0()<x.V2()) && x.V2()->IsRW()&& !x.V2()->IsV()){
              x.V2()->SetV();
              h_ret.push_back(HeapElem(new MYTYPE(VertexPair(x.V0(),x.V2()),TriEdgeCollapse< TriMeshType,VertexPair,MYTYPE>::GlobalMark(),_pp )));
            }
          }
        }
    }
    else
    { // if the collapse is A-symmetric (e.g. u->v != v->u)
      for(auto vi=m.vert.begin();vi!=m.vert.end();++vi)
        if(!(*vi).IsD() && (*vi).IsRW())
        {
          vcg::face::VFIterator<FaceType> x;
          UnMarkAll(m);
          for( x.F() = (*vi).VFp(), x.I() = (*vi).VFi(); x.F()!=0; ++ x)
          {
            if(x.V()->IsRW() && x.V1()->IsRW() && !IsMarked(m,x.F()->V1(x.I()))){
              h_ret.push_back( HeapElem( new MYTYPE( VertexPair (x.V(),x.V1()),TriEdgeCollapse< TriMeshType,VertexPair,MYTYPE>::GlobalMark(),_pp)));
            }
            if(x.V()->IsRW() && x.V2()->IsRW() && !IsMarked(m,x.F()->V2(x.I()))){
              h_ret.push_back( HeapElem( new MYTYPE( VertexPair (x.V(),x.V2()),TriEdgeCollapse< TriMeshType,VertexPair,MYTYPE>::GlobalMark(),_pp)));
            }
          }
        }
    }
  }
//  static float HeapSimplexRatio(BaseParameterClass *_pp) {return IsSymmetric(_pp)?5.0f:9.0f;}
  static float HeapSimplexRatio(BaseParameterClass *_pp) {return IsSymmetric(_pp)?4.0f:8.0f;}
  static bool IsSymmetric(BaseParameterClass *_pp) {return ((QParameter *)_pp)->OptimalPlacement;}
  static bool IsVertexStable(BaseParameterClass *_pp) {return !((QParameter *)_pp)->OptimalPlacement;}


/** Evaluate the priority (error) for an edge collapse
  *
  * It simulate the collapse and compute the quadric error 
  * generated by this collapse. This error is weighted with 
  * - aspect ratio of involved triangles
  * - normal variation
  */
  ScalarType ComputePriority(BaseParameterClass *_pp)
  {
    QParameter *pp=(QParameter *)_pp;
    
    VertexType * v[2];
    v[0] = this->pos.V(0);
    v[1] = this->pos.V(1);
    
    std::vector<CoordType> origNormalVec; // vector with incident faces original normals 
    if(pp->NormalCheck){ // Collect Original Normals
      for(VFIterator x(v[0]); !x.End(); ++x )	 // for all faces in v0
        if( x.V1()!=v[1] && x.V2()!=v[1] )     // skip faces with v1
          origNormalVec.push_back(NormalizedTriangleNormal(*x.F()));
      for(VFIterator x(v[1]); !x.End(); ++x )	 // for all faces in v1
        if( x.V1()!=v[0] && x.V2()!=v[0] )     // skip faces with v0
          origNormalVec.push_back(NormalizedTriangleNormal(*x.F()));
    }
    
    ScalarType origArea=0;
    if(pp->AreaCheck){ // Collect Original Area
      for(VFIterator x(v[0]); !x.End(); ++x )	 // for all faces in v0
        origArea += DoubleArea(*x.F());
      for(VFIterator x(v[1]); !x.End(); ++x )	 // for all faces in v1
        if( x.V1()!=v[0] && x.V2()!=v[0] )     // skip faces with v0
          origArea += DoubleArea(*x.F());
    }
    
    ScalarType origQual= std::numeric_limits<double>::max(); 
    if(pp->HardQualityCheck){ // Collect original quality
      for(VFIterator x(v[0]); !x.End(); ++x )	 // for all faces in v0
          origQual=std::min(origQual, QualityFace(*x.F()));
      for(VFIterator x(v[1]); !x.End(); ++x )	 // for all faces in v1
        if( x.V1()!=v[0] && x.V2()!=v[0] )     // skip faces with v0
          origQual=std::min(origQual, QualityFace(*x.F()));
    }
 
    
    //// Move the two vertexes into new position (storing the old ones)
    CoordType OldPos0=v[0]->P();
    CoordType OldPos1=v[1]->P();
    ComputePosition(_pp);      
    // Now Simulate the collapse 
    v[0]->P() = v[1]->P() =  this->optimalPos;    
     
    //// Rescan faces and compute the worst quality and normals that happens after collapse
    ScalarType MinCos  = std::numeric_limits<double>::max();  // Cosine of the angle variation: -1 ~ very bad to 1~perfect
    if(pp->NormalCheck){    
      int i=0;
      for(VFIterator x(v[0]); !x.End(); ++x )  // for all faces in v0
        if( x.V1()!=v[1] && x.V2()!=v[1] )     // skipping faces with v1
        {
          CoordType nn=NormalizedTriangleNormal(*x.F());
          MinCos=std::min(MinCos,nn.dot(origNormalVec[i++]));
        }
      for(VFIterator x(v[1]); !x.End(); ++x )	 // for all faces in v1
        if( x.V1()!=v[0] && x.V2()!=v[0] ) // skip faces with v0
        {
          CoordType nn=NormalizedTriangleNormal(*x.F());
          MinCos=std::min(MinCos,nn.dot(origNormalVec[i++]));
        }
    }      
    
    ScalarType newQual = std::numeric_limits<ScalarType>::max();  // 
    if(pp->QualityCheck){ 
      for(VFIterator x(v[0]); !x.End(); ++x )  // for all faces in v0
        if( x.V1()!=v[1] && x.V2()!=v[1] )   
          newQual=std::min(newQual,QualityFace(*x.F()));
      for(VFIterator x(v[1]); !x.End(); ++x )	 // for all faces in v1
        if( x.V1()!=v[0] && x.V2()!=v[0] ) // skip faces with v0
          newQual=std::min(newQual,QualityFace(*x.F()));
    }
            
    ScalarType newArea=0;
    if(pp->AreaCheck){ // Collect Area
      for(VFIterator x(v[0]); !x.End(); ++x )	 // for all faces in v0
          newArea += DoubleArea(*x.F());
      for(VFIterator x(v[1]); !x.End(); ++x )	 // for all faces in v1
        if( x.V1()!=v[0] && x.V2()!=v[0] )     // skip faces with v0
          newArea += DoubleArea(*x.F());
    }         
    
    QuadricType qq=QH::Qd(v[0]);
    qq+=QH::Qd(v[1]);

    double QuadErr = pp->ScaleFactor*qq.Apply(Point3d::Construct(v[1]->P()));
    
    assert(!math::IsNAN(QuadErr));
    // All collapses involving triangles with quality larger than <QualityThr> have no penalty;
    if(newQual>pp->QualityThr) newQual=pp->QualityThr;
    
    if(pp->NormalCheck){     
      // All collapses where the normal vary less than <NormalThr> (e.g. more than CosineThr)
      // have no particular penalty
      if(MinCos>pp->CosineThr) MinCos=pp->CosineThr;
      MinCos=fabs((MinCos+1.0)/2.0); // Now it is in the range 0..1 with 0 very bad!
      assert(MinCos >=0 && MinCos<1.1 );
    }

    
    QuadErr= std::max(QuadErr,pp->QuadricEpsilon);
    if(QuadErr <= pp->QuadricEpsilon) 
    {
      QuadErr *= Distance(OldPos0,OldPos1);  
    }

    if( pp->UseVertexWeight ) QuadErr *= (QH::W(v[1])+QH::W(v[0]))/2;
    
    ScalarType error;
    if(!pp->QualityCheck && !pp->NormalCheck) error = (ScalarType)(QuadErr);
    if( pp->QualityCheck && !pp->NormalCheck) error = (ScalarType)(QuadErr / newQual);
    if(!pp->QualityCheck &&  pp->NormalCheck) error = (ScalarType)(QuadErr / MinCos);
    if( pp->QualityCheck &&  pp->NormalCheck) error = (ScalarType)(QuadErr / (newQual*MinCos));

    if(pp->AreaCheck && ((fabs(origArea-newArea)/(origArea+newArea))>0.01) )
        error = std::numeric_limits<ScalarType>::max();
    
    if(pp->HardQualityCheck && 
       (newQual < pp->HardQualityThr && newQual < origQual*0.9) )
      error = std::numeric_limits<ScalarType>::max();
    
    if(pp->HardNormalCheck)
      if(CheckForFlip())  error = std::numeric_limits<ScalarType>::max();
    
    // Restore old position of v0 and v1
    v[0]->P()=OldPos0;
    v[1]->P()=OldPos1;
    
    this->_priority = error;
    return this->_priority;
  }
  
  
  bool CheckForFlippedFaceOverVertex(VertexType *vp, ScalarType angleThrRad =  math::ToRad(150.))
  {
    std::map<VertexType *, CoordType>  edgeNormMap; 
    ScalarType maxAngle=0;
    
    for(VFIterator x(vp); !x.End(); ++x )	 // for all faces in v1
    {
      if(QualityFace(*x.F()) <0.01 ) return true; 
        for(int i=0;i<2;++i)
        { 
          VertexType *vv= i==0?x.V1():x.V2();
          assert(vv!=vp);
          auto ni = edgeNormMap.find(vv);
          if(ni==edgeNormMap.end()) edgeNormMap[vv] = NormalizedTriangleNormal(*x.F());
          else maxAngle = std::max(maxAngle,AngleN(NormalizedTriangleNormal(*x.F()),ni->second));
        }
    }
    
    return (maxAngle > angleThrRad);          
  }  
  
  // This function return true if, after an edge collapse, 
  // among the surviving faces, there are two adjacent faces forming a 
  // diedral angle larger than the given threshold
  // It assumes that the two vertexes of the collapsing edge 
  // have been already moved to the new position but the topolgy has not yet been changed (e.g. there are two zero-area faces)
  
  bool CheckForFlip(ScalarType angleThrRad =  math::ToRad(150.))
  {
    std::map<VertexType *, CoordType>  edgeNormMap; 
    VertexType * v[2];
    v[0] = this->pos.V(0);
    v[1] = this->pos.V(1);
    ScalarType maxAngle=0;
    assert (v[0]->P()==v[1]->P());
    
    for(VFIterator x(v[0]); !x.End(); ++x )	 // for all faces in v0
      if( x.V1()!=v[1] && x.V2()!=v[1] )     // skip faces with v1
      {
        if(QualityFace(*x.F()) <0.01 ) return true;         
        for(int i=0;i<2;++i)
        { 
          VertexType *vv= (i==0)?x.V1():x.V2();
          assert(vv!=v[0]);
          auto ni = edgeNormMap.find(vv);
          if(ni==edgeNormMap.end()) edgeNormMap[vv] = NormalizedTriangleNormal(*x.F());
          else maxAngle = std::max(maxAngle,AngleN(NormalizedTriangleNormal(*x.F()),ni->second));
        }
      }
    for(VFIterator x(v[1]); !x.End(); ++x )	 // for all faces in v1
      if( x.V1()!=v[0] && x.V2()!=v[0] )     // skip faces with v0
      {
        if(QualityFace(*x.F()) <0.01 ) return true;         
        for(int i=0;i<2;++i)
        { 
          VertexType *vv= i==0?x.V1():x.V2();
          assert(vv!=v[1]);
          auto ni = edgeNormMap.find(vv);
          if(ni==edgeNormMap.end()) edgeNormMap[vv] = NormalizedTriangleNormal(*x.F());
          else maxAngle = std::max(maxAngle,AngleN(NormalizedTriangleNormal(*x.F()),ni->second));
        }
    }
    return (maxAngle > angleThrRad);      
  }
  

  
  
  inline void AddCollapseToHeap(HeapType & h_ret, VertexType *v0, VertexType *v1, BaseParameterClass *_pp)
  {
    QParameter *pp=(QParameter *)_pp;    
    ScalarType maxAdmitErr = std::numeric_limits<ScalarType>::max();
    h_ret.push_back(HeapElem(new MYTYPE(VertexPair(v0,v1), this->GlobalMark(),_pp)));
    if(h_ret.back().pri > maxAdmitErr) {
      delete h_ret.back().locModPtr;
      h_ret.pop_back(); 
    }
    else
      std::push_heap(h_ret.begin(),h_ret.end());
    
    if(!IsSymmetric(pp)){
      h_ret.push_back(HeapElem(new MYTYPE(VertexPair(v1,v0), this->GlobalMark(),_pp)));
      if(h_ret.back().pri > maxAdmitErr) {
        delete h_ret.back().locModPtr;
        h_ret.pop_back(); 
      }
      else      
        std::push_heap(h_ret.begin(),h_ret.end());
    }
  }
  
  inline  void UpdateHeap(HeapType & h_ret, BaseParameterClass *_pp)
  {
    this->GlobalMark()++;
    VertexType *v[2];
    v[0]= this->pos.V(0);
    v[1]= this->pos.V(1);
    v[1]->IMark() = this->GlobalMark();

    // First loop around the surviving vertex to unmark the Visit flags
    for(VFIterator vfi(v[1]); !vfi.End(); ++vfi ) {
      vfi.V1()->ClearV();
      vfi.V2()->ClearV();
      vfi.V1()->IMark() = this->GlobalMark();
      vfi.V2()->IMark() = this->GlobalMark();      
    }

    // Second Loop
    for(VFIterator vfi(v[1]); !vfi.End(); ++vfi ) {
      if( !(vfi.V1()->IsV()) && vfi.V1()->IsRW())
      {
        vfi.V1()->SetV();
        AddCollapseToHeap(h_ret,vfi.V0(),vfi.V1(),_pp);
      }
      if(  !(vfi.V2()->IsV()) && vfi.V2()->IsRW())
      {
        vfi.V2()->SetV();
        AddCollapseToHeap(h_ret,vfi.V2(),vfi.V0(),_pp);
      }
      if(vfi.V1()->IsRW() && vfi.V2()->IsRW() )
        AddCollapseToHeap(h_ret,vfi.V1(),vfi.V2(),_pp);
    } // end second loop around surviving vertex.
  }

  static void InitQuadric(TriMeshType &m,BaseParameterClass *_pp)
  {
    QParameter *pp=(QParameter *)_pp;
    QH::Init();
    for(VertexIterator pv=m.vert.begin();pv!=m.vert.end();++pv)
      if( ! (*pv).IsD() && (*pv).IsW())
        QH::Qd(*pv).SetZero();    
    
    for(FaceIterator fi=m.face.begin();fi!=m.face.end();++fi)
      if( !(*fi).IsD() && (*fi).IsR() )
        if((*fi).V(0)->IsR() &&(*fi).V(1)->IsR() &&(*fi).V(2)->IsR())
        {
          Plane3<ScalarType,false> facePlane;
          facePlane.SetDirection( ( (*fi).V(1)->cP() - (*fi).V(0)->cP() ) ^  ( (*fi).V(2)->cP() - (*fi).V(0)->cP() ));
          if(!pp->UseArea)
            facePlane.Normalize();
          facePlane.SetOffset( facePlane.Direction().dot((*fi).V(0)->cP()));                   

          QuadricType q;
          q.ByPlane(facePlane);          
          
          // The basic < add face quadric to each vertex > loop
          for(int j=0;j<3;++j)
            if( (*fi).V(j)->IsW() )
              QH::Qd((*fi).V(j)) += q;
          
          for(int j=0;j<3;++j)
            if( (*fi).IsB(j) || pp->QualityQuadric )
            {
              Plane3<ScalarType,false> borderPlane; 
              QuadricType bq;
              // Border quadric record the squared distance from the plane orthogonal to the face and passing 
              // through the edge. 
              borderPlane.SetDirection(facePlane.Direction() ^ (( (*fi).V1(j)->cP() - (*fi).V(j)->cP() ).normalized()));
              if(  (*fi).IsB(j) ) borderPlane.SetDirection(borderPlane.Direction()* (ScalarType)(pp->BoundaryQuadricWeight ));        // amplify border planes
              else                borderPlane.SetDirection(borderPlane.Direction()* (ScalarType)(pp->QualityQuadricWeight ));   // and consider much less quadric for quality
              borderPlane.SetOffset(borderPlane.Direction().dot((*fi).V(j)->cP()));
              bq.ByPlane(borderPlane);
              
              if( (*fi).V (j)->IsW() )	QH::Qd((*fi).V (j)) += bq;
              if( (*fi).V1(j)->IsW() )	QH::Qd((*fi).V1(j)) += bq;
            }
        }
    
    if(pp->ScaleIndependent)
    {
      vcg::tri::UpdateBounding<TriMeshType>::Box(m);
      //Make all quadric independent from mesh size
      pp->ScaleFactor = 1e8*pow(1.0/m.bbox.Diag(),6); // scaling factor
    }

    if(pp->QualityWeight) // we map quality range into a squared 01 and than this into the 1..QualityWeightFactor range
    {
      ScalarType minQ, maxQ;
      tri::Stat<TriMeshType>::ComputePerVertexQualityMinMax(m,minQ,maxQ);      
      for(VertexIterator vi=m.vert.begin();vi!=m.vert.end();++vi)
        if( ! (*vi).IsD() && (*vi).IsW())
        {
          const double quality01squared = pow((double)((vi->Q()-minQ)/(maxQ-minQ)),2.0);
          QH::Qd(*vi) *= 1.0 + quality01squared * (pp->QualityWeightFactor-1.0); 
        }
    }
  }
  
  CoordType ComputeMinimal()
  {
  }
  
CoordType ComputeMinimalOld()
{
   VertexType* &v0 = this->pos.V(0);
   VertexType* &v1 = this->pos.V(1);
   QuadricType q=QH::Qd(v0);
   q+=QH::Qd(v1);
   
   Point3<QuadricType::ScalarType> x;
   
   bool rt=q.Minimum(x);
   if(!rt) { // if the computation of the minimum fails we choose between the two edge points and the middle one.
     Point3<QuadricType::ScalarType> x0=Point3d::Construct(v0->P());
     Point3<QuadricType::ScalarType> x1=Point3d::Construct(v1->P());
     x.Import((v1->P()+v1->P())/2);
     double qvx=q.Apply(x);
     double qv0=q.Apply(x0);
     double qv1=q.Apply(x1);
     if(qv0<qvx) x=x0;
     if(qv1<qvx && qv1<qv0) x=x1;
   }
   
   return CoordType::Construct(x);
 }
};

} // namespace tri
} // namespace vcg
#endif
