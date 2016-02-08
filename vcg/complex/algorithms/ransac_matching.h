/****************************************************************************
* VCGLib                                                            o o     *
* Visual and Computer Graphics Library                            o     o   *
*                                                                _   O  _   *
* Copyright(C) 2004-2012                                           \/)\/    *
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

#ifndef RANSAC_MATCHING_H
#define RANSAC_MATCHING_H
#include<vcg/complex/algorithms/point_sampling.h>
#include<vcg/complex/algorithms/ransac_matching.h>
#include<vcg/space/index/kdtree/kdtree.h>
#include<vcg/space/point_matching.h>
namespace vcg
{


template <class MeshType> 
class BaseFeature
{
public:
  BaseFeature():_v(0) {}
  typename MeshType::VertexType *_v;  
  typename MeshType::CoordType P() {return _v->cP();}    
};
template <class MeshType>
class BaseFeatureSet 
{
public: 
  typedef   BaseFeature<MeshType> FeatureType;  
  typedef typename MeshType::VertexType      VertexType;
 
  std::vector<FeatureType> fixFeatureVec;
  std::vector<FeatureType> movFeatureVec;
  
  FeatureType &ff(int i) { return fixFeatureVec[i]; }
  FeatureType &mf(int i) { return movFeatureVec[i]; }
  int ffNum() const { return fixFeatureVec.size(); }
  
  void Init(MeshType &fix, MeshType &mov,  
            std::vector<VertexType *> &fixSampleVec, std::vector<VertexType *> &movSampleVec)
  {
    this->fixFeatureVec.resize(fixSampleVec.size()/20);
    for(int i=0;i<fixSampleVec.size()/20;++i) 
      this->fixFeatureVec[i]._v = fixSampleVec[i];
      
    this->movFeatureVec.resize(movSampleVec.size()/20);
    for(int i=0;i<movSampleVec.size()/20;++i) 
      this->movFeatureVec[i]._v = movSampleVec[i];
    
    printf("Generated %i Features on Fix\n",this->fixFeatureVec.size());
    printf("Generated %i Features on Mov\n",this->movFeatureVec.size());
  }
 
 // Returns the indexes of all the fix features matching a given one (from mov usually) 
 void getMatchingFeatureVec(FeatureType &q, vector<int> &mfiVec)
 {
  mfiVec.resize(movFeatureVec.size());
  
  for(int i=0;i<movFeatureVec.size();++i)
    mfiVec[i]=i;
 } 
};


template <class MeshType> 
class NDFeature 
{
public:
  typedef typename MeshType::ScalarType      ScalarType;
  
  typename MeshType::VertexType *_v;  
  typename MeshType::CoordType nd; //  
  
  typename MeshType::CoordType P() {return _v->cP();}    
  
  static void EvalNormalVariation(MeshType &m, ScalarType dist)
  {
    tri::UpdateNormal<MeshType>::PerVertexNormalized(m);
    
    VertexConstDataWrapper<MeshType > ww(m);
    KdTree<ScalarType> tree(ww); 
      
    for(int i=0;i<m.vn;++i)
    {
      std::vector<unsigned int> ptIndVec;
      std::vector<ScalarType> sqDistVec;    
      tree.doQueryDist(m.vert[i].P(),dist, ptIndVec, sqDistVec);
      ScalarType varSum=0;
      for(int j=0;j<sqDistVec.size();++j)
      {
        varSum += Distance(m.vert[i].N(),m.vert[ptIndVec[j]].N()); 
      }
     m.vert[i].Q()=varSum/ScalarType(ptIndVec.size());
    }  
    tri::UpdateColor<MeshType>::PerVertexQualityGray(m,0,0);
  }
  
};








template <class MeshType, class FeatureSetType>
class RansacFramework
{
  typedef typename FeatureSetType::FeatureType       FeatureType;
  typedef typename MeshType::CoordType       CoordType;
  typedef typename MeshType::BoxType         BoxType;
  typedef typename MeshType::ScalarType      ScalarType;
  typedef typename MeshType::VertexType      VertexType;
  typedef typename MeshType::VertexPointer   VertexPointer;
  typedef typename MeshType::VertexIterator  VertexIterator;
  typedef typename MeshType::EdgeType        EdgeType;
  typedef typename MeshType::EdgeIterator    EdgeIterator;
  typedef typename MeshType::FaceType        FaceType;
  typedef typename MeshType::FacePointer     FacePointer;
  typedef typename MeshType::FaceIterator    FaceIterator;
  typedef typename MeshType::FaceContainer   FaceContainer;
  typedef Matrix44<ScalarType>               Matrix44Type;
  
public:
  class Param
  {
  public:
    Param()
    {
      iterMax=100;
      samplingRadiusPerc=0.005;
      samplingRadiusAbs=0;
      evalSize=1000;
      inlierRatioThr=0.3;
      inlierDistanceThrPerc = 1.5; // the distance between a transformed mov sample and the corresponding on fix should be 1.5 * sampling dist.
      congruenceThrPerc = 2.0; // the distance between two matching features must be  within 2.0 * sampling distance 
      minFeatureDistancePerc = 4.0; // the distance between two chosen features must be  at least 4.0 * sampling distance 
    }
   
    ScalarType inlierRatioThr;
    ScalarType inlierDistanceThrPerc;
    ScalarType congruenceThrPerc;
    ScalarType minFeatureDistancePerc;
    ScalarType samplingRadiusPerc;
    ScalarType samplingRadiusAbs;
    int iterMax;
    int evalSize;
    
    ScalarType inlierSquareThr() const { return pow(samplingRadiusAbs* inlierDistanceThrPerc,2); }
  };
  
  class Candidate
  {
  public:
    int fixInd[3];
    int movInd[3];  
    int inlierNum;
    int evalSize;
    Matrix44Type Tr;
    ScalarType err() const {return float(inlierNum)/float(evalSize);}
    bool operator <(const Candidate &cc) const
    {
      return this->err() > cc.err();
    }
    
  };

  FeatureSetType FS;
  std::vector<Point3f> fixConsensusVec, movConsensusVec;
  KdTree<ScalarType> *consensusTree;
  
  
  // Given three pairs of sufficiently different distances (e.g. the edges of a scalene triangle)
  // it finds the permutation that brings the vertexes so that the distances match.
  // The meaning of the permutation vector nm0,nm1,nm2 is that the (N)ew index of (M)ov vertx i is the value of nmi 
  
  bool FindPermutation(int d01, int d02, int d12, int m01, int m02, int m12, int nm[], Param &pp)
  {
    ScalarType eps = pp.samplingRadiusAbs*2.0;
        
    if(fabs(d01-m01)<eps) {
      if(fabs(d02-m02)<eps) {
        if(fabs(d12-m12)<eps){ nm[0]=0;nm[1]=1;nm[2]=2; return true; }
            else return false;
      }
      if(fabs(d02-m12)<eps) {
        if(fabs(d12-m02)<eps){ nm[0]=1;nm[1]=0;nm[2]=2; return true; }
            else return false;        
      }        
    }
    
    if(fabs(d01-m02)<eps) {
      if(fabs(d02-m01)<eps) {
        if(fabs(d12-m12)<eps){ nm[0]=0;nm[1]=2;nm[2]=1; return true; }
            else return false;
      }
      if(fabs(d02-m12)<eps) {
        if(fabs(d12-m01)<eps){ nm[0]=2;nm[1]=0;nm[2]=1; return true; }
            else return false;        
      }        
    }

    if(fabs(d01-m12)<eps) {
      if(fabs(d02-m01)<eps) {
        if(fabs(d12-m02)<eps){ nm[0]=1;nm[1]=2;nm[2]=0; return true; }
            else return false;
      }
      if(fabs(d02-m02)<eps) {
        if(fabs(d12-m01)<eps){ nm[0]=2;nm[1]=1;nm[2]=0; return true; }
            else return false;        
      }        
    }
    return false;
  }
    
  // The main loop. 
  // Choose three points on fix that make a scalene triangle 
  // and search on mov three other points with matchng distances 
  
  void Process_SearchEvaluateTriple (vector<Candidate> &cVec, Param &pp)
  {
    math::MarsenneTwisterRNG rnd;
//    ScalarType congruenceEps = pow(pp.samplingRadiusAbs * pp.congruenceThrPerc,2.0f);
    ScalarType congruenceEps = pp.samplingRadiusAbs * pp.congruenceThrPerc;
    ScalarType minFeatureDistEps = pp.samplingRadiusAbs * pp.minFeatureDistancePerc;
    printf("Starting search congruenceEps = samplingRadiusAbs * 3.0 = %6.2f \n",congruenceEps);
    int iterCnt=0;
    
    while ( (iterCnt < pp.iterMax) && (cVec.size()<100) )
    {
      Candidate c;
      // Choose a random pair of features from fix 
      c.fixInd[0] = rnd.generate(FS.ffNum());
      c.fixInd[1] = rnd.generate(FS.ffNum());
      ScalarType d01 = Distance(FS.ff(c.fixInd[0]).P(),FS.ff(c.fixInd[1]).P());
      if( d01 > minFeatureDistEps )
      {
        c.fixInd[2] = rnd.generate(FS.ffNum());
        ScalarType d02=Distance(FS.ff(c.fixInd[0]).P(),FS.ff(c.fixInd[2]).P());
        ScalarType d12=Distance(FS.ff(c.fixInd[1]).P(),FS.ff(c.fixInd[2]).P());
          
        if( ( d02 > minFeatureDistEps ) &&  // Sample are sufficiently distant
            ( d12 > minFeatureDistEps ) && 
            ( fabs(d01-d02) > congruenceEps ) && // and they make a scalene triangle
            ( fabs(d01-d12) > congruenceEps ) && 
            ( fabs(d12-d02) > congruenceEps ) )
        {
          // Find a congruent triple on mov 
          printf("Starting search of a [%i] congruent triple for %4i %4i %4i - %6.2f %6.2f %6.2f\n",iterCnt,c.fixInd[0],c.fixInd[1],c.fixInd[2],d01,d02,d12);
          // As a first Step we ask for three vectors of matching features;
          
          std::vector<int> movFeatureVec0; 
          FS.getMatchingFeatureVec(FS.ff(c.fixInd[0]), movFeatureVec0);
          std::vector<int> movFeatureVec1; 
          FS.getMatchingFeatureVec(FS.ff(c.fixInd[1]), movFeatureVec1);
          std::vector<int> movFeatureVec2; 
          FS.getMatchingFeatureVec(FS.ff(c.fixInd[2]), movFeatureVec2);
          
          int congrNum=0;
          int congrGoodNum=0;
          for(int i=0;i<movFeatureVec0.size();++i)
          { 
            if(cVec.size()>100) break;
            c.movInd[0]=movFeatureVec0[i];
            for(int j=0;j<movFeatureVec1.size();++j)
            {               
              if(cVec.size()>100) break;             
              c.movInd[1]=movFeatureVec1[j];              
              ScalarType m01 = Distance(FS.mf(c.movInd[0]).P(),FS.mf(c.movInd[1]).P());
              if( (fabs(m01-d01)<congruenceEps) )
              {
//                printf("- Found a congruent pair %i %i %6.2f\n", c.movInd[0],c.movInd[1], m01);                
                ++congrNum;
                for(int k=0;k<movFeatureVec2.size();++k)
                { 
                  if(cVec.size()>100) break;                  
                  c.movInd[2]=movFeatureVec2[k];
                  ScalarType m02=Distance(FS.mf(c.movInd[0]).P(),FS.mf(c.movInd[2]).P());
                  ScalarType m12=Distance(FS.mf(c.movInd[1]).P(),FS.mf(c.movInd[2]).P());
                  if( (fabs(m02-d02)<congruenceEps)  && (fabs(m12-d12)<congruenceEps ) )
                  {
                    c.Tr = GenerateMatchingMatrix(c,pp);
                    
                    EvaluateMatrix(c,pp);
                    if(c.err() > pp.inlierRatioThr ){
                      printf("- - Found  %i th good congruent triple %i %i %i -- %f / %i \n", cVec.size(), c.movInd[0],c.movInd[1],c.movInd[2],c.err(),pp.evalSize);
                      ++congrGoodNum;                      
                      cVec.push_back(c);
                    }
                  }
                }
              }                
            }
          }
          printf("Completed Search of congruent triple (found %i / %i good/congruent)\n",congrGoodNum,congrNum);
        }               
      }
      ++iterCnt;
    } // end While

    printf("Found %i candidates \n",cVec.size());
    sort(cVec.begin(),cVec.end());
    printf("best candidate %f \n",cVec[0].err());
    
    pp.evalSize = pp.evalSize*10;
    
    for(int i=0;i<cVec.size();++i)
      EvaluateMatrix(cVec[i],pp);
    
    sort(cVec.begin(),cVec.end());
    
    printf("After re-evaluation best is %f",cVec[0].err());
        
      
    
  } // end Process
  
  int EvaluateMatrix(Candidate &c, Param &pp)
  {
    c.inlierNum=0;
    c.evalSize=pp.evalSize;    
    
    ScalarType sqThr = pp.inlierSquareThr();
    Distribution<ScalarType> H;
    for(int i=0;i<pp.evalSize;++i)
    {
      Point3f qp = c.Tr*movConsensusVec[i];
      uint ind;
      ScalarType squareDist;
      consensusTree->doQueryClosest(qp,ind,squareDist);
      if(squareDist < sqThr)
        ++c.inlierNum;
    }
  }
  
  int DumpInlier(MeshType &m, Candidate &c, Param &pp)
  {
    ScalarType sqThr = pp.inlierSquareThr();
    for(int i=0;i<pp.evalSize;++i)
    {
      Point3f qp = c.Tr*movConsensusVec[i];
      uint ind;
      ScalarType squareDist;
      consensusTree->doQueryClosest(qp,ind,squareDist);
      if(squareDist < sqThr)
        tri::Allocator<MeshType>::AddVertex(m,qp);
    }
    
  }
  

// Find the transformation that matches the mov onto the fix
// eg M * piMov = piFix 

Matrix44f GenerateMatchingMatrix(Candidate &c, Param pp)
{
  std::vector<Point3f> pFix(3);
  pFix[0]= FS.ff(c.fixInd[0]).P();
  pFix[1]= FS.ff(c.fixInd[1]).P();
  pFix[2]= FS.ff(c.fixInd[2]).P();
  
  std::vector<Point3f> pMov(3);
  pMov[0]= FS.mf(c.movInd[0]).P();
  pMov[1]= FS.mf(c.movInd[1]).P();
  pMov[2]= FS.mf(c.movInd[2]).P();

  Point3f upFix = vcg::Normal(pFix[0],pFix[1],pFix[2]);
  Point3f upMov = vcg::Normal(pMov[0],pMov[1],pMov[2]);  
  
  upFix.Normalize(); 
  upMov.Normalize();
  
  upFix *= Distance(pFix[0],pFix[1]);
  upMov *= Distance(pMov[0],pMov[1]);
  
  for(int i=0;i<3;++i) pFix.push_back(pFix[i]+upFix);
  for(int i=0;i<3;++i) pMov.push_back(pMov[i]+upMov);
  
  Matrix44f res;
  ComputeRigidMatchMatrix(pFix,pMov,res);
  return res;  
}


void Init(MeshType &fixM, MeshType &movM, Param &pp)
{
  // First a bit of Sampling
  typedef tri::TrivialPointerSampler<MeshType> BaseSampler;
  typename tri::SurfaceSampling<MeshType, BaseSampler>::PoissonDiskParam pdp;
  pdp.randomSeed = 0;
  pdp.bestSampleChoiceFlag = true;
  pdp.bestSamplePoolSize = 20;
  int t0=clock();
  pp.samplingRadiusAbs = pp.samplingRadiusPerc *fixM.bbox.Diag();
  BaseSampler pdSampler;
  std::vector<VertexType *> fixSampleVec;
  tri::SurfaceSampling<MeshType,BaseSampler>::PoissonDiskPruning(pdSampler, fixM, pp.samplingRadiusAbs,pdp);
  std::swap(pdSampler.sampleVec,fixSampleVec);
  std::vector<VertexType *> movSampleVec;  
  tri::SurfaceSampling<MeshType,BaseSampler>::PoissonDiskPruning(pdSampler, movM, pp.samplingRadiusAbs,pdp);
  std::swap(pdSampler.sampleVec,movSampleVec);
  int t1=clock();
  printf("Poisson Sampling of surfaces %5.2f ( %iv and %iv) \n",float(t1-t0)/CLOCKS_PER_SEC,fixSampleVec.size(),movSampleVec.size());
  printf("Sampling Radius %f \n",pp.samplingRadiusAbs);
  
  for(int i=0;i<fixSampleVec.size();++i) 
    this->fixConsensusVec.push_back(fixSampleVec[i]->P());    

  for(int i=0;i<movSampleVec.size();++i) 
    this->movConsensusVec.push_back(movSampleVec[i]->P());
  
  FS.Init(fixM, movM, fixSampleVec, movSampleVec);
    
  std::random_shuffle(movConsensusVec.begin(),movConsensusVec.end());
  
  VectorConstDataWrapper<std::vector<CoordType> > ww(fixConsensusVec);
  consensusTree = new  KdTree<ScalarType>(ww); 
}


};

} //end namespace vcg


#endif // RANSAC_MATCHING_H
