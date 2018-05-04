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
#ifndef CUT_TREE_H
#define CUT_TREE_H

#include<vcg/complex/complex.h>
#include <vcg/space/index/kdtree/kdtree.h>
#include<vcg/complex/algorithms/update/quality.h>
#include<vcg/complex/algorithms/update/color.h>

namespace vcg {
namespace tri {

template <class MeshType>
class CutTree
{
public:
  typedef typename MeshType::ScalarType     ScalarType;
  typedef typename MeshType::CoordType     CoordType;
  typedef typename MeshType::VertexType     VertexType;
  typedef typename MeshType::VertexPointer  VertexPointer;
  typedef typename MeshType::VertexIterator VertexIterator;
  typedef typename MeshType::EdgeIterator   EdgeIterator;
  typedef typename MeshType::EdgeType       EdgeType;
  typedef typename MeshType::FaceType       FaceType;
  typedef typename MeshType::FacePointer    FacePointer;
  typedef typename MeshType::FaceIterator   FaceIterator;
  typedef Box3<ScalarType>               Box3Type;
  typedef typename face::Pos<FaceType> PosType;
  typedef typename tri::UpdateTopology<MeshType>::PEdge PEdge;
  
  MeshType &base; 
  
  CutTree(MeshType &_m) :base(_m){}
  
  
// Perform a simple optimization of the three applying simple shortcuts:
// if the endpoints of two consecutive edges are connected by an edge existing on base mesh just use that edges
  
void OptimizeTree(KdTree<ScalarType> &kdtree, MeshType &t)
{
  tri::Allocator<MeshType>::CompactEveryVector(t);  
  int lastEn=t.en;
  do
  {
    lastEn=t.en;
    tri::UpdateTopology<MeshType>::VertexEdge(t);
    
    // First simple loop that search for 2->1 moves. 
    for(VertexIterator vi=t.vert.begin();vi!=t.vert.end();++vi)
    {
      std::vector<VertexType *> starVec;
      edge::VVStarVE(&*vi,starVec);
      if(starVec.size()==2)  // middle vertex has to be 1-manifold
      {
        PosType pos;
        if(ExistEdge(kdtree,starVec[0]->P(),starVec[1]->P(),pos))
          edge::VEEdgeCollapse(t,&*vi);
      }
    }
    tri::Allocator<MeshType>::CompactEveryVector(t);    
  }
  while(t.en<lastEn);
}

// Given two points return true if on the base mesh there exist an edge with that two coords
// if return true the pos indicate the found edge. 
bool ExistEdge(KdTree<ScalarType> &kdtree, CoordType &p0, CoordType &p1, PosType &fpos)
{
  ScalarType locEps = SquaredDistance(p0,p1)/100000.0;
  
  VertexType *v0=0,*v1=0;
  unsigned int veInd;
  ScalarType sqdist;
  kdtree.doQueryClosest(p0,veInd,sqdist);
  if(sqdist<locEps) 
    v0 = &base.vert[veInd];
  kdtree.doQueryClosest(p1,veInd,sqdist);
  if(sqdist<locEps) 
    v1 = &base.vert[veInd];
  if(v0 && v1)
  {
    fpos =PosType(v0->VFp(),v0);
    assert(fpos.V()==v0);
    PosType startPos=fpos;
    do
    {
      fpos.FlipE(); fpos.FlipF();
      if(fpos.VFlip()== v1) return true;
    } while(startPos!=fpos);    
  }
  return false;
}


int findNonVisitedEdgesDuringRetract(VertexType * vp, EdgeType * &ep)
{
  std::vector<EdgeType *> starVec;
  edge::VEStarVE(&*vp,starVec);
  int cnt =0;
  for(size_t i=0;i<starVec.size();++i) {
    if(!starVec[i]->IsV()) {
      cnt++;
      ep = starVec[i];
    }
  }  
  return cnt;
}

bool IsBoundaryVertexOnBase(KdTree<ScalarType> &kdtree, const CoordType &p)
{
  VertexType *v0=0;
  unsigned int veInd;
  ScalarType sqdist;
  kdtree.doQueryClosest(p,veInd,sqdist);
  if(sqdist>0) { assert(0); } 
  v0 = &base.vert[veInd];
  return v0->IsB();
}

/**
 * @brief Retract
 * @param t the edgemesh containing the visit tree. 
 *  
 */
void Retract(KdTree<ScalarType> &kdtree, MeshType &t)
{
  printf("Retracting a tree of %i edges and %i vertices\n",t.en,t.vn);
  tri::UpdateTopology<MeshType>::VertexEdge(t);
  tri::Allocator<MeshType>::CompactEveryVector(t);
  std::stack<VertexType *> vertStack;

  // Put on the stack all the vertex with just a single incident edge. 
  ForEachVertex(t, [&](VertexType &v){
    if(edge::VEDegree<EdgeType>(&v) ==1)  
      vertStack.push(&v);
  });
  
  tri::UpdateFlags<MeshType>::EdgeClearV(t);
  tri::UpdateFlags<MeshType>::VertexClearV(t);
  
  int unvisitedEdgeNum = t.en;
  while((!vertStack.empty()) && (unvisitedEdgeNum > 2) )
  {
    VertexType *vp = vertStack.top();
    vertStack.pop();
    vp->C()=Color4b::Blue;
    EdgeType *ep=0;
    int eCnt =  findNonVisitedEdgesDuringRetract(vp,ep);
    if(eCnt==1) // We have only one non visited edge over vp
    {
      assert(!ep->IsV());
      ep->SetV();
      --unvisitedEdgeNum;
      VertexType *otherVertP;
      if(ep->V(0)==vp) otherVertP = ep->V(1);
      else otherVertP = ep->V(0);
      vertStack.push(otherVertP);
    }
  }
  assert(unvisitedEdgeNum >0);
  for(size_t i =0; i<t.edge.size();++i){
    PosType fpos;
    if( ExistEdge(kdtree, t.edge[i].P(0), t.edge[i].P(1), fpos)){
      if(fpos.IsBorder()) {
        t.edge[i].SetV();
      }
    }
    else assert(0);
  }
  
  // All the boundary edges are in the initial tree so the clean boundary loops chains remains as irreducible loops
  // We delete them (leaving dangling edges with a vertex on the boundary)
  for(size_t i =0; i<t.edge.size();++i){    
    if (t.edge[i].IsV()) 
      tri::Allocator<MeshType>::DeleteEdge(t,t.edge[i]) ;
  }
  assert(t.en >0);
  tri::Clean<MeshType>::RemoveUnreferencedVertex(t);
  tri::Allocator<MeshType>::CompactEveryVector(t);
}

/** \brief Main function
 * 
 * It builds a cut tree that open the mesh into a topological disk
 * 
 * 
 */
void Build(MeshType &dualMesh, int startingFaceInd=0)
{
  tri::UpdateTopology<MeshType>::FaceFace(base);
  tri::UpdateTopology<MeshType>::VertexFace(base); 
  
  BuildVisitTree(dualMesh,startingFaceInd);
//  BuildDijkstraVisitTree(dualMesh,startingFaceInd);

  VertexConstDataWrapper<MeshType > vdw(base);
  KdTree<ScalarType> kdtree(vdw);  
  Retract(kdtree,dualMesh);  
  OptimizeTree(kdtree, dualMesh);
  tri::UpdateBounding<MeshType>::Box(dualMesh);      
}

/* Auxiliary class for keeping the heap of vertices to visit and their estimated distance */
  struct FaceDist{
    FaceDist(FacePointer _f):f(_f),dist(_f->Q()){}
    FacePointer f;
    ScalarType dist; 
    bool operator < (const FaceDist &o) const
    {
      if( dist != o.dist)
        return dist > o.dist;
      return f<o.f;
    }
  };


void BuildDijkstraVisitTree(MeshType &dualMesh, int startingFaceInd=0, ScalarType maxDistanceThr=std::numeric_limits<ScalarType>::max())
{
  tri::RequireFFAdjacency(base);
  tri::RequirePerFaceMark(base);
  tri::RequirePerFaceQuality(base);
  typename MeshType::template PerFaceAttributeHandle<FacePointer> parentHandle
      = tri::Allocator<MeshType>::template GetPerFaceAttribute<FacePointer>(base, "parent");

  std::vector<FacePointer> seedVec;
  seedVec.push_back(&base.face[startingFaceInd]);
   
  std::vector<FaceDist> Heap;
  tri::UnMarkAll(base);
  tri::UpdateQuality<MeshType>::FaceConstant(base,0);
  ForEachVertex(base, [&](VertexType &v){
    tri::Allocator<MeshType>::AddVertex(dualMesh,v.cP());
  });
  
  // Initialize the face heap; 
  // All faces in the heap are already marked; Q() store the distance from the source faces; 
  for(size_t i=0;i<seedVec.size();++i)
  {
    seedVec[i]->Q()=0;
    Heap.push_back(FaceDist(seedVec[i]));
  }
  // Main Loop
  int boundary=0;
  std::make_heap(Heap.begin(),Heap.end());
  
  int vCnt=0;
  int eCnt=0;
  int fCnt=0;
  
  // The main idea is that in the heap we maintain all the faces to be visited. 
  int nonDiskCnt=0;
  while(!Heap.empty() && nonDiskCnt<10)
  {
    int eulerChi= vCnt-eCnt+fCnt;
    if(eulerChi==1) nonDiskCnt=0;
    else ++nonDiskCnt;
//    printf("HeapSize %i: %i - %i + %i = %i\n",Heap.size(), vCnt,eCnt,fCnt,eulerChi);
    pop_heap(Heap.begin(),Heap.end());
    FacePointer currFp = (Heap.back()).f;
    if(tri::IsMarked(base,currFp))
    {
//      printf("Found an already visited face %f %f \n",Heap.back().dist, Heap.back().f->Q());
      //assert(Heap.back().dist != currFp->Q());
            
      Heap.pop_back(); 
      continue;
    }
    Heap.pop_back();
    ++fCnt;
    eCnt+=3;
    tri::Mark(base,currFp);
    
//    printf("pop face %i \n", tri::Index(base,currFp));
    for(int i=0;i<3;++i)
    {
      if(!currFp->V(i)->IsV()) {++vCnt; currFp->V(i)->SetV();}
      
      FacePointer nextFp = currFp->FFp(i);
      if( tri::IsMarked(base,nextFp) )
      {
        eCnt-=1;
        printf("is marked\n");
        if(nextFp != parentHandle[currFp] )
        {
          if(currFp>nextFp){
            tri::Allocator<MeshType>::AddEdge(dualMesh,tri::Index(base,currFp->V0(i)), tri::Index(base,currFp->V1(i)));
          }
        }
      }
      else // add it to the heap;
      {
//        printf("is NOT marked\n");
        parentHandle[nextFp] = currFp;
        ScalarType nextDist = currFp->Q() + Distance(Barycenter(*currFp),Barycenter(*nextFp));
        int adjMarkedNum=0; 
        for(int k=0;k<3;++k) if(tri::IsMarked(base,nextFp->FFp(k))) ++adjMarkedNum;
        if(nextDist < maxDistanceThr || adjMarkedNum>1)        
        {
          nextFp->Q() = nextDist;
          Heap.push_back(FaceDist(nextFp));
          push_heap(Heap.begin(),Heap.end());
        }
        else {
//          printf("boundary %i\n",++boundary);
          tri::Allocator<MeshType>::AddEdge(dualMesh,tri::Index(base,currFp->V0(i)), tri::Index(base,currFp->V1(i)));
        }
      }
    }
  } // End while
  printf("fulltree %i vn %i en \n",dualMesh.vn, dualMesh.en);
  int dupVert=tri::Clean<MeshType>::RemoveDuplicateVertex(dualMesh,false);   printf("Removed %i dup vert\n",dupVert);
  int dupEdge=tri::Clean<MeshType>::RemoveDuplicateEdge(dualMesh);   printf("Removed %i dup edges %i\n",dupEdge,dualMesh.EN());
  tri::Clean<MeshType>::RemoveUnreferencedVertex(dualMesh);   
  
  tri::io::ExporterPLY<MeshType>::Save(dualMesh,"fulltree.ply",tri::io::Mask::IOM_EDGEINDEX);   
  tri::UpdateColor<MeshType>::PerFaceQualityRamp(base);
  tri::io::ExporterPLY<MeshType>::Save(base,"colored_Bydistance.ply",tri::io::Mask::IOM_FACECOLOR);    
}

// \brief This function build a cut tree. 
//
// First we make a bread first FF face visit. 
// Each time that we encounter a visited face we add to the tree the edge 
// that brings to the already visited face.
// this structure build a dense graph and we retract this graph retracting each 
// leaf until we remains with just the loops that cuts the object. 

void BuildVisitTree(MeshType &dualMesh, int startingFaceInd=0)
{
  tri::UpdateFlags<MeshType>::FaceClearV(base);
  tri::UpdateFlags<MeshType>::VertexBorderFromFaceAdj(base);
  std::vector<face::Pos<FaceType> > visitStack; // the stack contain the pos on the 'starting' face. 
  
  base.face[startingFaceInd].SetV();
  for(int i=0;i<3;++i)
    visitStack.push_back(PosType(&(base.face[startingFaceInd]),i,base.face[startingFaceInd].V(i)));

  int cnt=1;
  
  while(!visitStack.empty())
  {
    std::swap(visitStack.back(),visitStack[rand()%visitStack.size()]);
    PosType c = visitStack.back();
    visitStack.pop_back();
    assert(c.F()->IsV());
    c.F()->C() = Color4b::ColorRamp(0,base.fn,cnt);
    c.FlipF();
    if(!c.F()->IsV())
    {
      ++cnt;
      c.F()->SetV();
      c.FlipE();c.FlipV();
      visitStack.push_back(c);
      c.FlipE();c.FlipV();
      visitStack.push_back(c);
    }
    else
    {
      tri::Allocator<MeshType>::AddEdge(dualMesh,c.V()->P(),c.VFlip()->P());
    }
  }
  assert(cnt==base.fn);
 
  tri::Clean<MeshType>::RemoveDuplicateVertex(dualMesh);    
  tri::io::ExporterPLY<MeshType>::Save(dualMesh,"fulltree.ply",tri::io::Mask::IOM_EDGEINDEX);    
} 

};
} // end namespace tri
} // end namespace vcg

#endif // CUT_TREE_H
