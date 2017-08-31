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
  typedef Box3  <ScalarType>               Box3Type;
  typedef typename vcg::GridStaticPtr<FaceType, ScalarType> MeshGrid;  
  typedef typename vcg::GridStaticPtr<EdgeType, ScalarType> EdgeGrid;
  typedef typename face::Pos<FaceType> PosType;
  typedef typename tri::UpdateTopology<MeshType>::PEdge PEdge;
  
  MeshType &base; 
//  MeshGrid uniformGrid;
  
//  Param par; 
  CutTree(MeshType &_m) :base(_m){}
  
  
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
      if(starVec.size()==2)
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


// \brief This function build a cut tree. 
//
// First we make a bread first FF face visit. 
// Each time that we encounter a visited face we add to the tree the edge 
// that brings to the already visited face.
// this structure build a dense graph and we retract this graph retracting each 
// leaf until we remains with just the loops that cuts the object. 

void BuildVisitTree(MeshType &dualMesh, int startingFaceInd=0)
{
  tri::UpdateTopology<MeshType>::FaceFace(base);
  tri::UpdateTopology<MeshType>::VertexFace(base); 
  
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
  
  VertexConstDataWrapper<MeshType > vdw(base);
  KdTree<ScalarType> kdtree(vdw);
  
  tri::Clean<MeshType>::RemoveDuplicateVertex(dualMesh);    
//  tri::io::ExporterPLY<MeshType>::Save(dualMesh,"fulltree.ply",tri::io::Mask::IOM_EDGEINDEX);  
  
  Retract(kdtree,dualMesh);
  OptimizeTree(kdtree, dualMesh);
  tri::UpdateBounding<MeshType>::Box(dualMesh);    
} 

};
} // end namespace tri
} // end namespace vcg

#endif // CUT_TREE_H
