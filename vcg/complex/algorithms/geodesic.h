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
#include <assert.h>
#include <vcg/math/base.h>
#include <vcg/container/simple_temporary_data.h>
#include <vcg/simplex/face/pos.h>
#include <vcg/simplex/face/topology.h>
#include <vcg/complex/algorithms/update/quality.h>
#include <deque>
#include <vector>
#include <list>
#include <functional>

/*
class for computing approximated geodesic distances on a mesh.

basic example: farthest vertex from a specified one
MyMesh m;
MyMesh::VertexPointer seed,far;
MyMesh::ScalarType dist;

vcg::Geo<MyMesh> g;
g.FarthestVertex(m,seed,far,d);

*/
#ifndef __VCGLIB_GEODESIC
#define __VCGLIB_GEODESIC

namespace vcg{
namespace tri{

template <class MeshType>
struct EuclideanDistance{
  typedef typename MeshType::VertexType VertexType;
  typedef typename MeshType::ScalarType  ScalarType;

  EuclideanDistance(){}
  ScalarType operator()(const VertexType * v0, const VertexType * v1) const
  {return vcg::Distance(v0->cP(),v1->cP());}
};

template <class MeshType, class DistanceFunctor = EuclideanDistance<MeshType> >
class Geo{

public:

  typedef typename MeshType::VertexType VertexType;
  typedef typename MeshType::VertexIterator VertexIterator;
  typedef typename MeshType::VertexPointer VertexPointer;
  typedef typename MeshType::FaceType  FaceType;
  typedef typename MeshType::CoordType  CoordType;
  typedef typename MeshType::ScalarType  ScalarType;



/* Auxiliary class for keeping the heap of vertices to visit and their estimated distance */
  struct VertDist{
    VertDist(){}
    VertDist(VertexPointer _v, ScalarType _d):v(_v),d(_d){}

    VertexPointer v;
    ScalarType d;
  };


  /* Temporary data to associate to all the vertices: estimated distance and boolean flag */
  struct TempData{
    TempData(){}
    TempData(const ScalarType & _d):d(_d),source(0),parent(0){}

    ScalarType d;
    VertexPointer source;//closest source
    VertexPointer parent;
  };

  typedef SimpleTempData<std::vector<VertexType>, TempData >  TempDataType;


  struct pred: public std::binary_function<VertDist,VertDist,bool>{
    pred(){}
    bool operator()(const VertDist& v0, const VertDist& v1) const
    {return (v0.d > v1.d);}
  };

  struct pred_addr: public std::binary_function<VertDist,VertDist,bool>{
    pred_addr(){}
    bool operator()(const VertDist& v0, const VertDist& v1) const
    {return (v0.v > v1.v);}
  };

  //************** calcolo della distanza di pw in base alle distanze note di pw1 e curr
  //************** sapendo che (curr,pw,pw1) e'una faccia della mesh
  //************** (vedi figura in file distance.gif)
  static ScalarType Distance(const VertexPointer &pw,
                             const VertexPointer &pw1,
                             const VertexPointer &curr,
                             const ScalarType &d_pw1,
                             const ScalarType &d_curr)
  {
    ScalarType curr_d=0;

    ScalarType ew_c  = DistanceFunctor()(pw,curr);
    ScalarType ew_w1 = DistanceFunctor()(pw,pw1);
    ScalarType ec_w1 = DistanceFunctor()(pw1,curr);
    CoordType w_c =  (pw->cP()-curr->cP()).Normalize() * ew_c;
    CoordType w_w1 = (pw->cP() - pw1->cP()).Normalize() * ew_w1;
    CoordType w1_c =  (pw1->cP() - curr->cP()).Normalize() * ec_w1;

    ScalarType	alpha,alpha_, beta,beta_,theta,h,delta,s,a,b;

    alpha = acos((w_c.dot(w1_c))/(ew_c*ec_w1));
    s = (d_curr + d_pw1+ec_w1)/2;
    a = s/ec_w1;
    b = a*s;
    alpha_ = 2*acos ( std::min<ScalarType>(1.0,sqrt(  (b- a* d_pw1)/d_curr)));

    if ( alpha+alpha_ > M_PI){
      curr_d = d_curr + ew_c;
    }else
    {
      beta_ = 2*acos ( std::min<ScalarType>(1.0,sqrt(  (b- a* d_curr)/d_pw1)));
      beta  = acos((w_w1).dot(-w1_c)/(ew_w1*ec_w1));

      if ( beta+beta_ > M_PI)
        curr_d = d_pw1  + ew_w1;
      else
      {
        theta	= ScalarType(M_PI)-alpha-alpha_;
        delta	= cos(theta)* ew_c;
        h		= sin(theta)* ew_c;
        curr_d = sqrt( pow(h,2)+ pow(d_curr + delta,2));
      }
    }
    return (curr_d);
  }

/*
This is the low level version of the geodesic computation framework.
Starting from the seeds, it assign a distance value to each vertex. The distance of a vertex is its
approximated geodesic distance to the closest seeds.
This is function is not meant to be called (although is not prevented). Instead, it is invoked by
wrapping function.
*/
  static  VertexPointer Visit(
      MeshType & m,
      std::vector<VertDist> & seedVec, // the set of seed to start from
      bool farthestOnBorder = false,
      ScalarType distance_threshold  = std::numeric_limits<ScalarType>::max(),                    // cut off distance (do no compute anything farther than this value)
      typename MeshType::template PerVertexAttributeHandle<VertexPointer> * vertSource = NULL,    // if present we put in this attribute the closest source for each vertex
      typename MeshType::template PerVertexAttributeHandle<VertexPointer> * vertParent = NULL,    // if present we put in this attribute the parent in the path that goes from the vertex to the closest source
      std::vector<VertexPointer> *InInterval=NULL)
  {
    std::vector<VertDist> frontier;
    VertexPointer farthest=0,pw,pw1;

    //Requirements
    assert(HasPerVertexVFAdjacency(m) && HasPerFaceVFAdjacency(m));
    assert(!seedVec.empty());

    TempDataType TD(m.vert, std::numeric_limits<ScalarType>::max());

    typename std::vector <VertDist >::iterator ifr;
    for(ifr = seedVec.begin(); ifr != seedVec.end(); ++ifr){
      (*ifr).d = 0.0;
      TD[(*ifr).v].d = 0.0;
      TD[(*ifr).v].source  = (*ifr).v;
      TD[(*ifr).v].parent  = (*ifr).v;
      frontier.push_back(VertDist((*ifr).v,0.0));
    }
    // initialize Heap
    make_heap(frontier.begin(),frontier.end(),pred());

    ScalarType curr_d,d_curr = 0.0,d_heap;
    ScalarType max_distance=0.0;

    while(!frontier.empty() && max_distance < distance_threshold)
    {
      pop_heap(frontier.begin(),frontier.end(),pred());
      VertexPointer curr = (frontier.back()).v;
      if (InInterval!=NULL) InInterval->push_back(curr);

      if(vertSource!=NULL)  (*vertSource)[curr] = TD[curr].source;
      if(vertParent!=NULL)  (*vertParent)[curr] = TD[curr].parent;

      d_heap = (frontier.back()).d;
      frontier.pop_back();

      assert(TD[curr].d <= d_heap);
      if(TD[curr].d < d_heap )// a vertex whose distance has been improved after it was inserted in the queue
        continue;
      assert(TD[curr].d == d_heap);

      d_curr =  TD[curr].d;

      bool isLeaf = (!farthestOnBorder || curr->IsB());

      face::VFIterator<FaceType> x;int k;

      for( x.f = curr->VFp(), x.z = curr->VFi(); x.f!=0; ++x )
        for(k=0;k<2;++k)
        {
          if(k==0) {
            pw = x.f->V1(x.z);
            pw1=x.f->V2(x.z);
          }
          else {
            pw = x.f->V2(x.z);
            pw1=x.f->V1(x.z);
          }

          const ScalarType & d_pw1  =  TD[pw1].d;
          {
            const ScalarType inter  = DistanceFunctor()(curr,pw1);//(curr->P() - pw1->P()).Norm();
            const ScalarType tol = (inter + d_curr + d_pw1)*.0001f;

            if (	(TD[pw1].source != TD[curr].source)||// not the same source
                    (inter + d_curr < d_pw1  +tol   ) ||
                    (inter + d_pw1  < d_curr +tol  ) ||
                    (d_curr + d_pw1  < inter +tol  )   // triangular inequality
                    )
              curr_d = d_curr + DistanceFunctor()(pw,curr);//(pw->P()-curr->P()).Norm();
            else
              curr_d = Distance(pw,pw1,curr,d_pw1,d_curr);
          }

          if(TD[pw].d > curr_d){
            TD[pw].d = curr_d;
            TD[pw].source = TD[curr].source;
            TD[pw].parent = curr;
            frontier.push_back(VertDist(pw,curr_d));
            push_heap(frontier.begin(),frontier.end(),pred());
          }
          if(isLeaf){
            if(d_curr > max_distance){
              max_distance = d_curr;
              farthest = curr;
            }
          }
        }
    }// end while

    // Copy found distance onto the Quality (\todo parametric!)
    if (InInterval==NULL)
    {
      for(VertexIterator vi = m.vert.begin(); vi != m.vert.end(); ++vi) if(!(*vi).IsD())
        (*vi).Q() =  TD[&(*vi)].d;
    }
    else
    {
      assert(InInterval->size()>0);
      for(size_t i=0;i<InInterval->size();i++)
        (*InInterval)[i]->Q() =  TD[(*InInterval)[i]].d;
    }

    return farthest;
  }


public:
  /*
      Given a mesh and a vector of pointers to seed vertices, this function assigns the approximated geodesic
      distance from the closest source to all the mesh vertices  within the
      specified interval and returns the found vertices writing on their Quality field the distance.
      Optionally for each vertex it can store, in a passed attribute, its corresponding seed vertex.
      To allocate such an attribute:

      typename MeshType::template PerVertexAttributeHandle<VertexPointer> sources;
      sources =  tri::Allocator<CMeshO>::AddPerVertexAttribute<VertexPointer> (m,"sources");

            */
  static bool FarthestVertex( MeshType & m,
                              std::vector<VertexPointer> & seedVec,
                              VertexPointer & farthest_vert,
                              ScalarType distance_thr  = std::numeric_limits<ScalarType>::max(),
                              typename MeshType::template PerVertexAttributeHandle<VertexPointer> * sourceSeed = NULL,
                              typename MeshType::template PerVertexAttributeHandle<VertexPointer> * parentSeed = NULL,
                              std::vector<VertexPointer> *InInterval=NULL)
  {
    typename std::vector<VertexPointer>::iterator fi;
    std::vector<VertDist> vdSeedVec;
    if(seedVec.empty())	return false;
    for( fi  = seedVec.begin(); fi != seedVec.end() ; ++fi)
        vdSeedVec.push_back(VertDist(*fi,0.0));
    farthest_vert = Visit(m, vdSeedVec, false, distance_thr, sourceSeed, parentSeed, InInterval);
    return true;
  }
  /*
  Given a mesh and  a  pointers to a vertex-source (source), assigns the approximated geodesic
  distance from the vertex-source to all the mesh vertices and returns the pointer to the farthest
  Note: it updates the field Q() of the vertices
  */
  static bool FarthestVertex( MeshType & m, VertexPointer seed, ScalarType distance_thr  = std::numeric_limits<ScalarType>::max())
  {
    std::vector<VertexPointer>  seedVec(1,seed);
    VertexPointer v0;
    return FarthestVertex(m,seedVec,v0,distance_thr);
  }


  /*
  Same as FarthestPoint but the returned pointer is to a border vertex
  Note: update the field Q() of the vertices
  */
  static void FarthestBVertex(MeshType & m,
                              std::vector<VertexPointer> & seedVec,
                              VertexPointer & farthest,
                              ScalarType & distance,
                              typename MeshType::template PerVertexAttributeHandle<VertexPointer> * sources = NULL
      )
  {
    std::vector<VertDist>fr;
    for(typename std::vector<VertexPointer>::iterator fi  = seedVec.begin(); fi != seedVec.end() ; ++fi)
      fr.push_back(VertDist(*fi,0));
    farthest =  Visit(m,fr,distance,true,sources);
  }
  /*
            Same as FarthestPoint but the returned pointer is to a border vertex
            Note: update the field Q() of the vertices
            */
  static void FarthestBVertex(	MeshType & m,
                                VertexPointer seed,
                                VertexPointer & farthest,
                                ScalarType & distance,
                                typename MeshType::template PerVertexAttributeHandle<VertexPointer> * sources = NULL)
  {
    std::vector<VertexPointer>  fro(1,seed);
    VertexPointer v0;
    FarthestBVertex(m,fro,v0,distance,sources);
    farthest = v0;
  }

  /*
            Assigns to each vertex of the mesh its distance to the closest vertex on the border
            Note: update the field Q() of the vertices
            Note: it needs the border bit set.
            */
  static bool DistanceFromBorder(	MeshType & m, typename MeshType::template PerVertexAttributeHandle<VertexPointer> * sources = NULL
      ){
    std::vector<VertexPointer> fro;
    VertexIterator vi;
    VertexPointer farthest;
    for(vi = m.vert.begin(); vi != m.vert.end(); ++vi)
      if( (*vi).IsB())
        fro.push_back(&(*vi));
    if(fro.empty()) return false;

    tri::UpdateQuality<MeshType>::VertexConstant(m,0);

    return FarthestVertex(m,fro,farthest,std::numeric_limits<ScalarType>::max(),sources);
  }

};// end class
}// end namespace tri
}// end namespace vcg
#endif
