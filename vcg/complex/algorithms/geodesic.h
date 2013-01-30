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
#include <vcg/simplex/face/pos.h>
#include <vcg/simplex/face/topology.h>
#include <vcg/complex/algorithms/update/quality.h>
#include <deque>
#include <functional>
#ifndef __VCGLIB_GEODESIC
#define __VCGLIB_GEODESIC

namespace vcg{
namespace tri{

template <class MeshType>
struct EuclideanDistance{
  typedef typename MeshType::VertexType VertexType;
  typedef typename MeshType::ScalarType  ScalarType;
  typedef typename MeshType::FacePointer FacePointer;

  EuclideanDistance(){}

  ScalarType operator()(const VertexType * v0, const VertexType * v1) const
  {return vcg::Distance(v0->cP(),v1->cP());}

  ScalarType operator()(const FacePointer f0, const FacePointer f1) const
  {return vcg::Distance(Barycenter(*f0),Barycenter(*f1));}
};

/*! \brief class for computing approximate geodesic distances on a mesh

  require VF Adjacency relation
\sa trimesh_geodesic.cpp
*/

template <class MeshType, class DistanceFunctor = EuclideanDistance<MeshType> >
class Geodesic{

public:

  typedef typename MeshType::VertexType VertexType;
  typedef typename MeshType::VertexIterator VertexIterator;
  typedef typename MeshType::VertexPointer VertexPointer;
  typedef typename MeshType::FacePointer FacePointer;
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


  struct DIJKDist{
    DIJKDist(VertexPointer _v):v(_v){}
    VertexPointer v;

    bool operator < (const DIJKDist &o) const
    {
      if( v->Q() != o.v->Q())
        return v->Q() > o.v->Q();
      return v<o.v;
    }
   };

  /* Auxiliary class for keeping the heap of vertices to visit and their estimated distance */
    struct FaceDist{
      FaceDist(FacePointer _f):f(_f){}
      FacePointer f;
      bool operator < (const FaceDist &o) const
      {
        if( f->Q() != o.f->Q())
          return f->Q() > o.f->Q();
        return f<o.f;
      }
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
    if(!HasVFAdjacency(m)) throw vcg::MissingComponentException("VFAdjacency");
    if(!HasPerVertexQuality(m)) throw vcg::MissingComponentException("VertexQuality");
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
      if(TD[curr].d < d_heap ) // a vertex whose distance has been improved after it was inserted in the queue
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
  /*! \brief Given a set of source vertices compute the approximate geodesic distance to all the other vertices

\param m the mesh
\param seedVec a vector of Vertex pointers with the \em sources of the flood fill
\param maxDistanceThr max distance that we travel on the mesh starting from the sources
\param withinDistanceVec a pointer to a vector for storing the vertexes reached within the passed maxDistanceThr
\param sourceSeed pointer to the handle to keep for each vertex its seed
\param parentSeed pointer to the handle to keep for each vertex its parent in the closest tree

Given a mesh and a vector of pointers to seed vertices, this function compute the approximated geodesic
distance from the given sources to all the mesh vertices within the given maximum distance threshold.
The computed distance is stored in the vertex::Quality component.
Optionally for each vertex it can store, in a passed attribute, the corresponding seed vertex
(e.g. the vertex of the source set closest to him) and the 'parent' in a tree forest that connects each vertex to the closest source.

To allocate the attributes:
\code
      typename MeshType::template PerVertexAttributeHandle<VertexPointer> sourcesHandle;
      sourcesHandle =  tri::Allocator<CMeshO>::AddPerVertexAttribute<MeshType::VertexPointer> (m,"sources");

      typename MeshType::template PerVertexAttributeHandle<VertexPointer> parentHandle;
      parentHandle =  tri::Allocator<CMeshO>::AddPerVertexAttribute<MeshType::VertexPointer> (m,"parent");
\endcode

It requires VF adjacency relation (e.g. vertex::VFAdj and face::VFAdj components)
It requires per vertex Quality (e.g. vertex::Quality component)

\warning that this function has ALWAYS at least a linear cost (it use additional attributes that have a linear initialization)
\todo make it O(output) by using incremental mark and persistent attributes.
            */
  static bool Compute( MeshType & m,
                       const std::vector<VertexPointer> & seedVec,
                       ScalarType maxDistanceThr  = std::numeric_limits<ScalarType>::max(),
                       std::vector<VertexPointer> *withinDistanceVec=NULL,
                       typename MeshType::template PerVertexAttributeHandle<VertexPointer> * sourceSeed = NULL,
                       typename MeshType::template PerVertexAttributeHandle<VertexPointer> * parentSeed = NULL
                       )
  {
    if(seedVec.empty())	return false;
    std::vector<VertDist> vdSeedVec;
    typename std::vector<VertexPointer>::const_iterator fi;
    for( fi  = seedVec.begin(); fi != seedVec.end() ; ++fi)
        vdSeedVec.push_back(VertDist(*fi,0.0));

    Visit(m, vdSeedVec, false, maxDistanceThr, sourceSeed, parentSeed, withinDistanceVec);
    return true;
  }

  /* \brief Assigns to each vertex of the mesh its distance to the closest vertex on the boundary

It is just a simple wrapper of the basic Compute()

            Note: update the field Q() of the vertices
            Note: it needs the border bit set.
            */
  static bool DistanceFromBorder(	MeshType & m, typename MeshType::template PerVertexAttributeHandle<VertexPointer> * sources = NULL)
  {
    std::vector<VertexPointer> fro;
    for(VertexIterator vi = m.vert.begin(); vi != m.vert.end(); ++vi)
      if( (*vi).IsB())
        fro.push_back(&(*vi));
    if(fro.empty()) return false;

    tri::UpdateQuality<MeshType>::VertexConstant(m,0);
    return Compute(m,fro,std::numeric_limits<ScalarType>::max(),0,sources);
  }


  static bool ConvertPerVertexSeedToPerFaceSeed(MeshType &m, const std::vector<VertexPointer> &vertexSeedVec,
                                                 std::vector<FacePointer> &faceSeedVec)
  {
    tri::RequireVFAdjacency(m);
    tri::RequirePerFaceMark(m);

    faceSeedVec.clear();
    tri::UnMarkAll(m);
    for(size_t i=0;i<vertexSeedVec.size();++i)
    {
      for(face::VFIterator<FaceType> vfi(vertexSeedVec[i]);!vfi.End();++vfi)
      {
        if(tri::IsMarked(m,vfi.F())) return false;
        faceSeedVec.push_back(vfi.F());
        tri::Mark(m,vfi.F());
      }
    }
    return true;
  }


  static void PerFaceDijsktraCompute(MeshType &m, const std::vector<FacePointer> &seedVec,
                                     ScalarType maxDistanceThr  = std::numeric_limits<ScalarType>::max(),
                                     std::vector<FacePointer> *InInterval=NULL,
                                     FacePointer FaceTarget=NULL,
                                     bool avoid_selected=false)
  {
    tri::RequireFFAdjacency(m);
    tri::RequirePerFaceMark(m);
    tri::RequirePerFaceQuality(m);

    typename MeshType::template PerFaceAttributeHandle<FacePointer> sourceHandle
        = tri::Allocator<MeshType>::template GetPerFaceAttribute<FacePointer> (m,"sources");

    typename MeshType::template PerFaceAttributeHandle<FacePointer> parentHandle
        = tri::Allocator<MeshType>::template GetPerFaceAttribute<FacePointer> (m,"parent");

    std::vector<FaceDist> Heap;
    tri::UnMarkAll(m);
    for(size_t i=0;i<seedVec.size();++i)
    {
      tri::Mark(m,seedVec[i]);
      seedVec[i]->Q()=0;
      sourceHandle[seedVec[i]]=seedVec[i];
      parentHandle[seedVec[i]]=seedVec[i];
      Heap.push_back(FaceDist(seedVec[i]));
      if (InInterval!=NULL) InInterval->push_back(seedVec[i]);
    }

    std::make_heap(Heap.begin(),Heap.end());
    while(!Heap.empty())
    {
      pop_heap(Heap.begin(),Heap.end());
      FacePointer curr = (Heap.back()).f;
      if ((FaceTarget!=NULL)&&(curr==FaceTarget))return;
      Heap.pop_back();

      for(int i=0;i<3;++i)
      {
        if(!face::IsBorder(*curr,i) )
        {
          FacePointer nextF = curr->FFp(i);
          ScalarType nextDist = curr->Q() + DistanceFunctor()(curr,nextF);
          if( (nextDist < maxDistanceThr) &&
              (!tri::IsMarked(m,nextF) ||  nextDist < nextF->Q()) )
          {
            nextF->Q() = nextDist;
            if ((avoid_selected)&&(nextF->IsS()))continue;
            tri::Mark(m,nextF);
            Heap.push_back(FaceDist(nextF));
            push_heap(Heap.begin(),Heap.end());
            if (InInterval!=NULL) InInterval->push_back(nextF);
            sourceHandle[nextF] = sourceHandle[curr];
            parentHandle[nextF] = curr;
//            printf("Heapsize %i nextDist = %f curr face %i next face %i \n",Heap.size(), nextDist, tri::Index(m,curr), tri::Index(m,nextF));
          }
        }
      }
    }
  }




  static void PerVertexDijsktraCompute(MeshType &m, const std::vector<VertexPointer> &seedVec,
                                     ScalarType maxDistanceThr  = std::numeric_limits<ScalarType>::max(),
                                     std::vector<VertexPointer> *InInterval=NULL,bool avoid_selected=false,
                                     VertexPointer target=NULL)
  {
    tri::RequireVFAdjacency(m);
    tri::RequirePerVertexMark(m);
    tri::RequirePerVertexQuality(m);

    typename MeshType::template PerVertexAttributeHandle<VertexPointer> sourceHandle
        = tri::Allocator<MeshType>::template GetPerVertexAttribute<VertexPointer> (m,"sources");

    typename MeshType::template PerVertexAttributeHandle<VertexPointer> parentHandle
        = tri::Allocator<MeshType>::template GetPerVertexAttribute<VertexPointer> (m,"parent");

    std::vector<DIJKDist> Heap;
    tri::UnMarkAll(m);

    for(size_t i=0;i<seedVec.size();++i)
    {
      assert(!tri::IsMarked(m,seedVec[i]));
      tri::Mark(m,seedVec[i]);
      seedVec[i]->Q()=0;
      sourceHandle[seedVec[i]]=seedVec[i];
      parentHandle[seedVec[i]]=seedVec[i];
      Heap.push_back(DIJKDist(seedVec[i]));
      if (InInterval!=NULL) InInterval->push_back(seedVec[i]);
    }

    std::make_heap(Heap.begin(),Heap.end());
    while(!Heap.empty())
    {
      pop_heap(Heap.begin(),Heap.end());
      VertexPointer curr = (Heap.back()).v;
      if ((target!=NULL)&&(target==curr))return;
      Heap.pop_back();
      std::vector<VertexPointer> vertVec;
      face::VVStarVF<FaceType>(curr,vertVec);
      for(size_t i=0;i<vertVec.size();++i)
      {
        VertexPointer nextV = vertVec[i];
        if ((avoid_selected)&&(nextV->IsS()))continue;
        ScalarType nextDist = curr->Q() + DistanceFunctor()(curr,nextV);
        if( (nextDist < maxDistanceThr) &&
            (!tri::IsMarked(m,nextV) ||  nextDist < nextV->Q()) )
        {
          nextV->Q() = nextDist;
          tri::Mark(m,nextV);
          Heap.push_back(DIJKDist(nextV));
          push_heap(Heap.begin(),Heap.end());
          if (InInterval!=NULL) InInterval->push_back(nextV);
          sourceHandle[nextV] = sourceHandle[curr];
          parentHandle[nextV] = curr;
//          printf("Heapsize %i nextDist = %f curr vert %i next vert %i \n",Heap.size(), nextDist, tri::Index(m,curr), tri::Index(m,nextV));
        }
      }
    }
  }


};// end class
}// end namespace tri
}// end namespace vcg
#endif
