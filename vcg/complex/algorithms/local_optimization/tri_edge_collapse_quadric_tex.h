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




#ifndef __VCG_TRIMESHCOLLAPSE_QUADRIC_TEX_
#define __VCG_TRIMESHCOLLAPSE_QUADRIC_TEX_

#include <vcg/complex/algorithms/local_optimization.h>
#include <vcg/complex/algorithms/local_optimization/tri_edge_collapse_quadric.h>
#include <vcg/container/simple_temporary_data.h>
#include <vcg/math/quadric5.h>
namespace vcg
{
namespace tri
{


class TriEdgeCollapseQuadricTexParameter : public BaseParameterClass
{
public:
  double  BoundaryWeight;
  double  CosineThr;
  float   ExtraTCoordWeight;
  bool    NormalCheck;
  double	  NormalThrRad;
  bool    OptimalPlacement;
  bool		  PreserveBoundary;
  bool		  PreserveTopology;
  double	  QuadricEpsilon;
  double	  QualityThr;
  bool	    QualityQuadric;
  bool    SafeHeapUpdate;
  double	  ScaleFactor;
  bool	    ScaleIndependent;
  bool	    UseArea;
  bool	    UseVertexWeight;

  TriEdgeCollapseQuadricTexParameter()
  {
    SetDefaultParams();
  }

  void SetDefaultParams()
  {
    this->BoundaryWeight=.5;
    this->CosineThr = cos(M_PI/2);
    this->ExtraTCoordWeight=0.0;
    this->NormalCheck=false;
    this->NormalThrRad=M_PI/2;
    this->OptimalPlacement=true;
    this->PreserveBoundary = false;
    this->PreserveTopology = false;
    this->QuadricEpsilon =1e-15;
    this->QualityThr=.1;
    this->QualityQuadric = false;
    this->SafeHeapUpdate=false;
    this->ScaleFactor=1.0;
    this->ScaleIndependent=true;
    this->UseArea=true;
    this->UseVertexWeight=false;
  }
};


// This is a static class that contains the references for the simple temporary data for the current mesh.
// for each vertex we keep a classic Quadric3D  and a vector of pairs texcoord+Quadric5D

template <class MeshType>
class QuadricTexHelper
    {
    public:
  typedef typename MeshType::VertexType VertexType;

  typedef	SimpleTempData<typename MeshType::VertContainer, std::vector<std::pair<vcg::TexCoord2f ,Quadric5<double> > > > Quadric5Temp;
  typedef	SimpleTempData<typename MeshType::VertContainer, math::Quadric<double> > QuadricTemp;

      QuadricTexHelper(){}


        static void Init(){}

    // it allocs the std::pair for the vertex relativly to the texture coord parameter
    static void Alloc(VertexType *v,vcg::TexCoord2f &coord)
    {
       std::vector<std::pair<vcg::TexCoord2f ,Quadric5<double> > > &qv = Vd(v);
       Quadric5<double> newq5;
       newq5.Zero();
       vcg::TexCoord2f newcoord;
       newcoord.u() = coord.u();
       newcoord.v() = coord.v();

       newq5.Sum3(Qd3(v),coord.u(),coord.v());

       qv.push_back(std::pair<vcg::TexCoord2f ,Quadric5<double> >(newcoord,newq5));
    }

    static void SumAll(VertexType *v,vcg::TexCoord2f &coord, Quadric5<double>& q)
    {
       std::vector<std::pair<vcg::TexCoord2f ,Quadric5<double> > > &qv = Vd(v);

       for(size_t i = 0; i < qv.size(); i++)
       {
         vcg::TexCoord2f &f = qv[i].first;
         if((f.u() == coord.u()) && (f.v() == coord.v()))
           qv[i].second += q;
         else
           qv[i].second.Sum3(Qd3(v),f.u(),f.v());
       }
    }

    static bool Contains(VertexType *v,vcg::TexCoord2f &coord)
    {
       std::vector<std::pair<vcg::TexCoord2f ,Quadric5<double> > > &qv = Vd(v);

       for(size_t i = 0; i < qv.size(); i++)
       {
         vcg::TexCoord2f &f = qv[i].first;
         if((f.u() == coord.u()) && (f.v() == coord.v()))
           return true;
       }

       return false;
    }

    static Quadric5<double> &Qd(VertexType *v,const vcg::TexCoord2f &coord)
    {
       std::vector<std::pair<vcg::TexCoord2f ,Quadric5<double> > > &qv = Vd(v);
       //assert(coord.N()>=0);
       for(size_t i = 0; i < qv.size(); i++)
       {
         vcg::TexCoord2f &f = qv[i].first;
         if((f.u() == coord.u()) && (f.v() == coord.v()))
           return qv[i].second;
       }

       assert(0);       
       return qv[0].second; 
    }
      static math::Quadric<double> &Qd3(VertexType *v) {return TD3()[*v];}
    static math::Quadric<double> &Qd3(VertexType &v) {return TD3()[v];}

    static std::vector<std::pair<vcg::TexCoord2f ,Quadric5<double> > > &Vd(VertexType *v){return (TD()[*v]);}
      static typename VertexType::ScalarType W(VertexType * /*v*/) {return 1.0;}
      static typename VertexType::ScalarType W(VertexType & /*v*/) {return 1.0;}
      static void Merge(VertexType & /*v_dest*/, VertexType const & /*v_del*/){}
      static  Quadric5Temp* &TDp() {static  Quadric5Temp *td; return td;}
      static  Quadric5Temp &TD() {return *TDp();}
      static  QuadricTemp* &TDp3() {static  QuadricTemp *td3; return td3;}
      static  QuadricTemp &TD3() {return *TDp3();}
    };







template<class TriMeshType, class VertexPair, class MYTYPE, class HelperType = tri::QInfoStandard<typename TriMeshType::VertexType>  >
class TriEdgeCollapseQuadricTex: public vcg::tri::TriEdgeCollapse< TriMeshType, VertexPair, MYTYPE>
{
  typedef HelperType QH;
  typedef typename tri::TriEdgeCollapse<TriMeshType, VertexPair, MYTYPE>::HeapType HeapType;
  typedef typename tri::TriEdgeCollapse<TriMeshType, VertexPair, MYTYPE>::HeapElem HeapElem;
  typedef typename TriMeshType::FaceType FaceType;
  typedef typename TriMeshType::VertexType VertexType;
  typedef typename TriMeshType::CoordType CoordType;
  typedef typename TriMeshType::CoordType::ScalarType ScalarType;
  typedef typename TriMeshType::VertexPointer VertexPointer;


  public:

  inline TriEdgeCollapseQuadricTex(const VertexPair &p, int mark, BaseParameterClass *_pp)
  {
      TriEdgeCollapseQuadricTexParameter *pp = (TriEdgeCollapseQuadricTexParameter *)_pp;
      this->localMark = mark;
      this->pos=p;
      this->_priority = ComputePriority(pp);
  }

// puntatori ai vertici che sono stati messi non-w per preservare il boundary
  static std::vector<VertexPointer>  & WV(){
      static std::vector<VertexPointer> _WV; return _WV;
    };



  static TriEdgeCollapseQuadricTexParameter & Params(){static TriEdgeCollapseQuadricTexParameter p; return p;}

      // Final Clean up after the end of the simplification process
    static void Finalize(TriMeshType &m,HeapType & /*h_ret*/, BaseParameterClass *_pp)
    {
      TriEdgeCollapseQuadricTexParameter *pp = (TriEdgeCollapseQuadricTexParameter *)_pp;
      vcg::tri::UpdateFlags<TriMeshType>::FaceBorderFromVF(m);

      // If we had the boundary preservation we should clean up the writable flags
      if(pp->PreserveBoundary)
      {
        typename 	std::vector<VertexPointer>::iterator wvi;
        for(wvi=WV().begin();wvi!=WV().end();++wvi)
          if(!(*wvi)->IsD()) (*wvi)->SetW();
      }
    }




  inline static int matchVertexID(FaceType *f,VertexType *v)
  {
    if(f->V(0)==v) return 0;
    if(f->V(1)==v) return 1;
    if(f->V(2)==v) return 2;

    assert(0); return -1;
  }

  inline int GetTexCoords(vcg::TexCoord2f &tcoord0_1, vcg::TexCoord2f &tcoord1_1,vcg::TexCoord2f &tcoord0_2,vcg::TexCoord2f &tcoord1_2)
  {
    int ncoords = 0;
    tcoord0_1.P()=Point2f(0.5f,0.5f);
    tcoord1_1.P()=Point2f(0.5f,0.5f);
    tcoord0_2.P()=Point2f(0.5f,0.5f);
    tcoord1_2.P()=Point2f(0.5f,0.5f);

    vcg::face::VFIterator<FaceType> vfi(this->pos.V(0));

    for(vfi.F() = this->pos.V(0)->VFp(), vfi.I() = this->pos.V(0)->VFi(); vfi.F()!=0; ++vfi )	// for all faces in v0
      if(vfi.F()->V(0)==this->pos.V(1) || vfi.F()->V(1)==this->pos.V(1) || vfi.F()->V(2)==this->pos.V(1) ) // and v1
      {
        if(ncoords == 0)
        {
          tcoord0_1 = vfi.F()->WT(matchVertexID(vfi.F(),this->pos.V(0)));
          tcoord1_1 = vfi.F()->WT(matchVertexID(vfi.F(),this->pos.V(1)));
        }
        else
        {
          tcoord0_2 = vfi.F()->WT(matchVertexID(vfi.F(),this->pos.V(0)));
          tcoord1_2 = vfi.F()->WT(matchVertexID(vfi.F(),this->pos.V(1)));

          if(
            (tcoord0_1.u() == tcoord0_2.u()) &&
            (tcoord0_1.v() == tcoord0_2.v()) &&
            (tcoord1_1.u() == tcoord1_2.u()) &&
            (tcoord1_1.v() == tcoord1_2.v())
            )
            return 1;
          else
            return 2;
        }
        ncoords++;
      }

    return ncoords;
  }


    ScalarType ComputePriority(BaseParameterClass *_pp)
    {
      TriEdgeCollapseQuadricTexParameter *pp = (TriEdgeCollapseQuadricTexParameter *)_pp;
      Quadric5<double> qsum1;
      Quadric5<double> qsum2;
      double min1[5];
      double min2[5];
      vcg::TexCoord2f tcoord0_1;
      vcg::TexCoord2f tcoord1_1;
      vcg::TexCoord2f tcoord0_2;
      vcg::TexCoord2f tcoord1_2;
      int ncoords;

      ncoords = GetTexCoords(tcoord0_1,tcoord1_1,tcoord0_2,tcoord1_2);

      return (ScalarType)ComputeMinimalsAndPriority(min1,min2,qsum1,qsum2,tcoord0_1,tcoord1_1,tcoord0_2,tcoord1_2,ncoords,pp);
    }


    /*
    * the very important function that says how good is a collapse.
    Originally is should be just the quadric error (e.g. you should choose the collapse that make the minimal quadric error)
    but important correcting factors has to be applyed
    - quality of the involved triangles
    - normal checking
    */
    ScalarType ComputeTexPriority(const double vv[5],Quadric5<double> &qsum, BaseParameterClass *_pp)
    {
      TriEdgeCollapseQuadricTexParameter *pp = (TriEdgeCollapseQuadricTexParameter *)_pp;
      VertexType * v[2];
      v[0] = this->pos.V(0);
      v[1] = this->pos.V(1);

      assert(!math::IsNAN(vv[0]));
      assert(!math::IsNAN(vv[1]));
      assert(!math::IsNAN(vv[2]));
      assert(!math::IsNAN(vv[3]));
      assert(!math::IsNAN(vv[4]));

      //// Move the two vertexe  into new position (storing the old ones)
      CoordType OldPos0=v[0]->P();
      CoordType OldPos1=v[1]->P();

      v[0]->P() = CoordType(vv[0],vv[1],vv[2]);
      //v[0]->P() = (v[0]->P()+v[1]->P())/2.0;
      v[1]->P()=v[0]->P();

      double QuadErr = qsum.Apply(vv);

      //// Rescan faces and compute quality and difference between normals
      double qt,   MinQual = 1e100;
      vcg::face::VFIterator<FaceType> x(this->pos.V(0));

      double ndiff,MinCos  = 1e100; // minimo coseno di variazione di una normale della faccia
                                  // (e.g. max angle) Mincos varia da 1 (normali coincidenti) a
                                  // -1 (normali opposte);

      for(x.F() = v[0]->VFp(), x.I() = v[0]->VFi(); x.F()!=0; ++x )	// for all faces in v0
        if(x.F()->V(0)!=v[1] && x.F()->V(1)!=v[1] && x.F()->V(2)!=v[1] )		// skip faces with v1
        {
          qt= QualityFace(*x.F());
          if(qt<MinQual) MinQual=qt;
          if(pp->NormalCheck){
              CoordType nn=TriangleNormal(*x.F()).Normalize();
              ndiff=nn.dot(x.F()->N()) / x.F()->N().Norm();
              if(ndiff<MinCos) MinCos=ndiff;
              assert(!math::IsNAN(ndiff));
              }
        }
      for(x.F() = v[1]->VFp(), x.I() = v[1]->VFi(); x.F()!=0; ++x )		// for all faces in v1
        if(x.F()->V(0)!=v[0] && x.F()->V(1)!=v[0] && x.F()->V(2)!=v[0] )			// skip faces with v0
        {
          qt= QualityFace(*x.F());
          if(qt<MinQual) MinQual=qt;
          if(pp->NormalCheck){
              CoordType nn=TriangleNormal(*x.F()).Normalize();
              ndiff=nn.dot(x.F()->N() / x.F()->N().Norm());
              if(ndiff<MinCos) MinCos=ndiff;
              assert(!math::IsNAN(ndiff));
          }

        }


      // All collapses involving triangles with quality larger than <QualityThr> have no penalty;
      if(MinQual>pp->QualityThr) MinQual=pp->QualityThr;
      if(QuadErr<1e-15) QuadErr=1e-15; // Do not use zero quality penalties

      this->_priority = (ScalarType)(QuadErr / MinQual);  // the priority of collapses that create low quality triangles has a penalty (it is increased)


      if(pp->NormalCheck){
          if(MinCos<pp->CosineThr) this->_priority *=1000;  // gross penalty to collapses that move too much the original normals.
      }



      //Restore old position of v0 and v1
      v[0]->P()=OldPos0;
      v[1]->P()=OldPos1;
      return this->_priority;
    }

    inline ScalarType ComputeMinimalsAndPriority(double dest_1[5],
                  double dest_2[5],
                  Quadric5<double> &qsum_1,
                  Quadric5<double> &qsum_2,
                  const vcg::TexCoord2f &tcoord0_1,
                  const vcg::TexCoord2f &tcoord1_1,
                  const vcg::TexCoord2f &tcoord0_2,
                  const vcg::TexCoord2f &tcoord1_2,
                   int ncoords,
                   BaseParameterClass *_pp)
    {
      TriEdgeCollapseQuadricTexParameter *pp = (TriEdgeCollapseQuadricTexParameter *)_pp;
      double tmp1[5];
      double tmp2[5];
      ScalarType priority1;
      ScalarType priority2;

      assert(!math::IsNAN(tcoord0_1.u()));
      assert(!math::IsNAN(tcoord0_1.v()));
      assert(!math::IsNAN(tcoord1_1.u()));
      assert(!math::IsNAN(tcoord1_1.v()));
      assert(!math::IsNAN(tcoord0_2.u()));
      assert(!math::IsNAN(tcoord0_2.v()));
      assert(!math::IsNAN(tcoord1_2.u()));
      assert(!math::IsNAN(tcoord1_2.v()));


      tmp1[0] = this->pos.V(0)->P().X();
      tmp1[1] = this->pos.V(0)->P().Y();
      tmp1[2] = this->pos.V(0)->P().Z();
      tmp1[3] = tcoord0_1.u();
      tmp1[4] = tcoord0_1.v();

      tmp2[0] = this->pos.V(1)->P().X();
      tmp2[1] = this->pos.V(1)->P().Y();
      tmp2[2] = this->pos.V(1)->P().Z();
      tmp2[3] = tcoord1_1.u();
      tmp2[4] = tcoord1_1.v();

      assert(QH::Qd(this->pos.V(0),tcoord0_1).IsValid());
      assert(QH::Qd(this->pos.V(1),tcoord1_1).IsValid());

      qsum_1 = QH::Qd(this->pos.V(0),tcoord0_1);
      qsum_1 += QH::Qd(this->pos.V(1),tcoord1_1);

      ComputeMinimal(dest_1,tmp1,tmp2,qsum_1,pp);
      priority1 = ComputeTexPriority(dest_1,qsum_1,pp);

      if(ncoords < 2)
        return priority1*(1 + (pp->ExtraTCoordWeight)*(QH::Vd(this->pos.V(0)).size()+ QH::Vd(this->pos.V(1)).size() - 2));


      tmp1[3] = tcoord0_2.u();
      tmp1[4] = tcoord0_2.v();

      tmp2[3] = tcoord1_2.u();
      tmp2[4] = tcoord1_2.v();

      assert(QH::Qd(this->pos.V(0),tcoord0_2).IsValid());
      assert(QH::Qd(this->pos.V(1),tcoord1_2).IsValid());

      qsum_2 = QH::Qd(this->pos.V(0),tcoord0_2);
      qsum_2 += QH::Qd(this->pos.V(1),tcoord1_2);

      ComputeMinimal(dest_2,tmp1,tmp2,qsum_2,pp);
      priority2 = ComputeTexPriority(dest_2,qsum_2,pp);

      if(priority1 > priority2)
      {
        ComputeMinimalWithGeoContraints(dest_2,tmp1,tmp2,qsum_2,dest_1,pp);
        priority2 = ComputeTexPriority(dest_2,qsum_2,pp);
      }
      else
      {
        ComputeMinimalWithGeoContraints(dest_1,tmp1,tmp2,qsum_1,dest_2,pp);
        priority1 = ComputeTexPriority(dest_1,qsum_1,pp);
      }


      this->_priority = std::max(priority1, priority2)*(1 + (pp->ExtraTCoordWeight)*(QH::Vd(this->pos.V(0)).size()+QH::Vd(this->pos.V(1)).size() - 2));

      return this->_priority;
    }

    inline void ComputeMinimal(double vv[5],const double v0[5],const double v1[5], const Quadric5<double> qsum,BaseParameterClass *_pp)
    {
      tri::TriEdgeCollapseQuadricTexParameter *pp =(tri::TriEdgeCollapseQuadricTexParameter *)_pp;
      bool rt=qsum.Minimum(vv);
      // if the computation of the minimum fails we choose between the two edge points and the middle one.
      // Switch to this branch also in the case of not using the optimal placement.
        if(!rt || !pp->OptimalPlacement )
        {
          vv[0] = (v0[0] + v1[0])/2;
          vv[1] = (v0[1] + v1[1])/2;
          vv[2] = (v0[2] + v1[2])/2;
          vv[3] = (v0[3] + v1[3])/2;
          vv[4] = (v0[4] + v1[4])/2;

          // In the case of not using the optimal placement we have to be sure that the middle value is discarded.
          double qvx= std::numeric_limits<float>::max();
          if(pp->OptimalPlacement)
            qvx = qsum.Apply(vv);


          double qv0=qsum.Apply(v0);
          double qv1=qsum.Apply(v1);


          if(qv0<qvx)
          {
            vv[0] = v0[0];
            vv[1] = v0[1];
            vv[2] = v0[2];
            vv[3] = v0[3];
            vv[4] = v0[4];
          }

          if(qv1<qvx && qv1<qv0)
          {
            vv[0] = v1[0];
            vv[1] = v1[1];
            vv[2] = v1[2];
            vv[3] = v1[3];
            vv[4] = v1[4];
          }
        }

        assert(!math::IsNAN(vv[0]));
        assert(!math::IsNAN(vv[1]));
        assert(!math::IsNAN(vv[2]));
        assert(!math::IsNAN(vv[3]));
        assert(!math::IsNAN(vv[4]));
    }

    inline void ComputeMinimalWithGeoContraints(double vv[5],const double v0[5],const double v1[5], const Quadric5<double> qsum, const double geo[5],BaseParameterClass *_pp)
    {
    tri::TriEdgeCollapseQuadricTexParameter *pp =(tri::TriEdgeCollapseQuadricTexParameter *)_pp;
      bool rt=qsum.MinimumWithGeoContraints(vv,geo);
      // if the computation of the minimum fails we choose between the two edge points and the middle one.
      // Switch to this branch also in the case of not using the optimal placement.
      if(!rt || !pp->OptimalPlacement) {
          double qvx = std::numeric_limits<float>::max();
          vv[0] = geo[0];
          vv[1] = geo[1];
          vv[2] = geo[2];
          if(pp->OptimalPlacement)
            {
              vv[3] = (v0[3] + v1[3])/2;
              vv[4] = (v0[4] + v1[4])/2;
              qvx=qsum.Apply(vv);
            }
          vv[3] = v0[3];
          vv[4] = v0[4];

          double qv0=qsum.Apply(vv);

          vv[3] = v1[3];
          vv[4] = v1[4];

          double qv1=qsum.Apply(v1);

          vv[3] = (v0[3] + v1[3])/2;
          vv[4] = (v0[4] + v1[4])/2;

          if(qv0<qvx)
          {
            vv[3] = v0[3];
            vv[4] = v0[4];
          }

          if(qv1<qvx && qv1<qv0)
          {
            vv[3] = v1[3];
            vv[4] = v1[4];
          }
        }

    }

  static void InitQuadric(TriMeshType &m,BaseParameterClass *_pp)
  {
  tri::TriEdgeCollapseQuadricTexParameter *pp =(tri::TriEdgeCollapseQuadricTexParameter *)_pp;
    typename TriMeshType::FaceIterator pf;
    HelperType::Init();

    for(pf=m.face.begin();pf!=m.face.end();++pf)
      if( !(*pf).IsD() && (*pf).IsR() )
        if((*pf).V(0)->IsR() &&(*pf).V(1)->IsR() &&(*pf).V(2)->IsR())
            {
              Quadric5<double> q;
              q.byFace(*pf, QH::Qd3((*pf).V(0)), QH::Qd3((*pf).V(1)), QH::Qd3((*pf).V(2)),pp->QualityQuadric,pp->BoundaryWeight);

              for(int j=0;j<3;++j)
                if( (*pf).V(j)->IsW())
                {
                  if(!HelperType::Contains((*pf).V(j),(*pf).WT(j)))
                  {
                  HelperType::Alloc((*pf).V(j),(*pf).WT(j));
                  }
                  assert(!math::IsNAN((*pf).WT(j).u()));
                  assert(!math::IsNAN((*pf).WT(j).v()));
                  HelperType::SumAll((*pf).V(j),(*pf).WT(j),q);
                }

            }

  }

    static void Init(TriMeshType &m,HeapType&h_ret,BaseParameterClass *_pp)
    {
    tri::TriEdgeCollapseQuadricTexParameter *pp =(tri::TriEdgeCollapseQuadricTexParameter *)_pp;
    typename 	TriMeshType::VertexIterator  vi;
    typename 	TriMeshType::FaceIterator  pf;

      vcg::tri::UpdateTopology<TriMeshType>::VertexFace(m);
      vcg::tri::UpdateFlags<TriMeshType>::FaceBorderFromVF(m);

    if(pp->PreserveBoundary )
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

      InitQuadric(m,pp);

    // Initialize the heap with all the possible collapses
        for(vi=m.vert.begin();vi!=m.vert.end();++vi)
          if((*vi).IsRW())
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
                          h_ret.push_back(HeapElem(new MYTYPE(VertexPair(x.V0(),x.V1()),TriEdgeCollapse< TriMeshType,VertexPair,MYTYPE>::GlobalMark(),pp )));
                          }
                    if((x.V0()<x.V2()) && x.V2()->IsRW()&& !x.V2()->IsV()){
                          x.V2()->SetV();
                          h_ret.push_back(HeapElem(new MYTYPE(VertexPair(x.V0(),x.V2()),TriEdgeCollapse< TriMeshType,VertexPair,MYTYPE>::GlobalMark(),pp )));
                        }
                  }
          }
    make_heap(h_ret.begin(),h_ret.end());
  }

  inline  void UpdateHeap(HeapType & h_ret,BaseParameterClass *_pp)
  {
    tri::TriEdgeCollapseQuadricTexParameter *pp =(tri::TriEdgeCollapseQuadricTexParameter *)_pp;
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

            h_ret.push_back(HeapElem(new MYTYPE(VertexPair(vfi.V0(),vfi.V1()), this->GlobalMark(),pp)));
            std::push_heap(h_ret.begin(),h_ret.end());
          }

          if(  !(vfi.V2()->IsV()) && vfi.V2()->IsRW())
          {
            vfi.V2()->SetV();

            h_ret.push_back(HeapElem(new MYTYPE(VertexPair(vfi.V0(),vfi.V2()),this->GlobalMark(),pp)));
            std::push_heap(h_ret.begin(),h_ret.end());
          }
        }
        ++vfi;
    }
  }

  void Execute(TriMeshType &m, BaseParameterClass *_pp)
  {
  tri::TriEdgeCollapseQuadricTexParameter *pp =(tri::TriEdgeCollapseQuadricTexParameter *)_pp;
  Quadric5<double> qsum1;
  Quadric5<double> qsum2;
  double min1[5];
  double min2[5];
  vcg::TexCoord2f tcoord0_1;
  vcg::TexCoord2f tcoord1_1;
  vcg::TexCoord2f tcoord0_2;
  vcg::TexCoord2f tcoord1_2;
  vcg::TexCoord2f newtcoord1;
  vcg::TexCoord2f newtcoord2;
  std::vector<std::pair<vcg::TexCoord2f ,Quadric5<double> > > qv;
  int ncoords;
  VertexType * v[2];
  v[0] = this->pos.V(0);
  v[1] = this->pos.V(1);

  vcg::math::Quadric<double> qsum3 = QH::Qd3(v[0]);
  qsum3 += QH::Qd3(v[1]);

  ncoords = GetTexCoords(tcoord0_1,tcoord1_1,tcoord0_2,tcoord1_2);

  ComputeMinimalsAndPriority(min1,min2,qsum1,qsum2,tcoord0_1,tcoord1_1,tcoord0_2,tcoord1_2,ncoords,pp);

   CoordType newPos(min1[0],min1[1],min1[2]); /* it's the same as min2[0],min2[1],min2[2] since the geometrical
                        constraint has been imposed during the re-computation of the other minimal */


  EdgeCollapser<TriMeshType,VertexPair>::Do(m, this->pos, newPos);
  //DoCollapse(m, this->pos, newPos ); // v0 is deleted and v1 take the new position

  vcg::TexCoord2f newtcoord;
  Quadric5<double> newq;



  newtcoord.u() = (float)min1[3];
  newtcoord.v() = (float)min1[4];
  assert(!math::IsNAN(newtcoord.u()));
  assert(!math::IsNAN(newtcoord.v()));
  newtcoord1 = newtcoord;
  newq = qsum1;

  qv.push_back(std::pair<vcg::TexCoord2f ,Quadric5<double> >(newtcoord,newq));

  if(ncoords > 1)
  {
    newtcoord.u() = min2[3];
    newtcoord.v() = min2[4];
    newtcoord2 = newtcoord;
    newq = qsum2;

    qv.push_back(std::pair<vcg::TexCoord2f ,Quadric5<double> >(newtcoord2,newq));
  }


  vcg::face::VFIterator<FaceType> vfi(v[1]);
  while (!vfi.End())
  {
    vcg::TexCoord2f & tcoords = vfi.F()->WT(matchVertexID(vfi.F(),v[1]));

    if(
      ((tcoords.u() == tcoord0_1.u()) && (tcoords.v() == tcoord0_1.v())) ||
      ((tcoords.u() == tcoord1_1.u()) && (tcoords.v() == tcoord1_1.v()))
      )
    {
      tcoords.u() = newtcoord1.u();
      tcoords.v() = newtcoord1.v();
    }
    else if(
      (ncoords > 1) &&
      (
      ((tcoords.u() == tcoord0_2.u()) && (tcoords.v() == tcoord0_2.v())) ||
      ((tcoords.u() == tcoord1_2.u()) && (tcoords.v() == tcoord1_2.v()))
      )
      )
    {
      tcoords.u()= newtcoord2.u();
      tcoords.v()= newtcoord2.v();
    }
    else
    {
      newtcoord = tcoords;

      if(QH::Contains(v[0],tcoords))
      {
        newq = QH::Qd(v[0],tcoords);
        newq.Sum3(QH::Qd3(v[1]),tcoords.u(),tcoords.v());
      }
      else if(QH::Contains(v[1],tcoords))
      {
        newq = QH::Qd(v[1],tcoords);
        newq.Sum3(QH::Qd3(v[0]),tcoords.u(),tcoords.v());
      }
      else
        assert(0);

      qv.push_back(std::pair<vcg::TexCoord2f ,Quadric5<double> >(newtcoord,newq));
    }

    ++vfi;
  }
  QH::Qd3(v[1]) = qsum3;
  QH::Vd(v[1]) = qv;

  }

};



  } // namespace tri
    } // namespace vcg
#endif
