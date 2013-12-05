/****************************************************************************
* MeshLab                                                           o o     *
* A versatile mesh processing toolbox                             o     o   *
*                                                                _   O  _   *
* Copyright(C) 2005                                                \/)\/    *
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

#ifndef VORONOI_PROCESSING_H
#define VORONOI_PROCESSING_H

#include <vcg/complex/algorithms/geodesic.h>
#include <vcg/complex/algorithms/update/color.h>
namespace vcg
{
namespace tri
{

template <class MeshType>
class ClusteringSampler
{
public:
  typedef typename MeshType::VertexType	VertexType;

  ClusteringSampler(std::vector<VertexType *> &_vec): sampleVec(_vec)
  {
    sampleVec = _vec;
  }

  std::vector<VertexType *> &sampleVec;

  void AddVert(const VertexType &p)
  {
    sampleVec.push_back((VertexType *)(&p));
  }
}; // end class ClusteringSampler


struct VoronoiProcessingParameter
{
  enum {
    None=0,
    DistanceFromSeed=1,
    DistanceFromBorder=2,
    RegionArea=3
  };

  VoronoiProcessingParameter()
  {
    colorStrategy = DistanceFromSeed;
    areaThresholdPerc=0;
    deleteUnreachedRegionFlag=false;
    fixSelectedSeed=false;
    collapseShortEdge=false;
    collapseShortEdgePerc = 0.01f;
    triangulateRegion=false;
    unbiasedSeedFlag = true;
    geodesicRelaxFlag = true;
  }
  int colorStrategy;
  float areaThresholdPerc;
  bool deleteUnreachedRegionFlag;
  bool unbiasedSeedFlag;
  bool fixSelectedSeed;   /// the vertexes that are selected are used as fixed seeds:
                          /// They will not move during relaxing
                          /// and they will be always included in the final triangulation (if on a border).
  bool triangulateRegion;
  bool collapseShortEdge;
  float collapseShortEdgePerc;
  bool geodesicRelaxFlag;
};

template <class MeshType, class DistanceFunctor = EuclideanDistance<MeshType> >
class VoronoiProcessing
{
  typedef typename MeshType::CoordType				CoordType;
  typedef typename MeshType::ScalarType				ScalarType;
  typedef typename MeshType::VertexType				VertexType;
  typedef typename MeshType::VertexPointer		VertexPointer;
  typedef typename MeshType::VertexIterator		VertexIterator;
  typedef typename MeshType::FacePointer			FacePointer;
  typedef typename MeshType::FaceIterator			FaceIterator;
  typedef typename MeshType::FaceType					FaceType;
  typedef typename MeshType::FaceContainer		FaceContainer;
  typedef typename tri::Geodesic<MeshType>::VertDist VertDist;

public:


// Given a vector of point3f it finds the closest vertices on the mesh.
static void SeedToVertexConversion(MeshType &m,std::vector<CoordType> &seedPVec,std::vector<VertexType *> &seedVVec)
{
    typedef typename vcg::SpatialHashTable<VertexType, ScalarType> HashVertexGrid;
    seedVVec.clear();

    HashVertexGrid HG;
    HG.Set(m.vert.begin(),m.vert.end());

    const float dist_upper_bound=m.bbox.Diag()/10.0;

    typename std::vector<CoordType>::iterator pi;
    for(pi=seedPVec.begin();pi!=seedPVec.end();++pi)
        {
            float dist;
            VertexPointer vp;
            vp=tri::GetClosestVertex<MeshType,HashVertexGrid>(m, HG, *pi, dist_upper_bound, dist);
            if(vp)
                {
                    seedVVec.push_back(vp);
                }
        }
}

typedef typename MeshType::template PerVertexAttributeHandle<VertexPointer> PerVertexPointerHandle;
typedef typename MeshType::template PerFaceAttributeHandle<VertexPointer> PerFacePointerHandle;

static void ComputePerVertexSources(MeshType &m, std::vector<VertexType *> &seedVec, DistanceFunctor &df)
{
  tri::Allocator<MeshType>::DeletePerVertexAttribute(m,"sources"); // delete any conflicting handle regardless of the type...
  PerVertexPointerHandle vertexSources =  tri::Allocator<MeshType>:: template AddPerVertexAttribute<VertexPointer> (m,"sources");

  tri::Allocator<MeshType>::DeletePerFaceAttribute(m,"sources"); // delete any conflicting handle regardless of the type...
  PerFacePointerHandle faceSources =  tri::Allocator<MeshType>:: template AddPerFaceAttribute<VertexPointer> (m,"sources");

  assert(tri::Allocator<MeshType>::IsValidHandle(m,vertexSources));

  tri::Geodesic<MeshType>::Compute(m,seedVec,df,std::numeric_limits<ScalarType>::max(),0,&vertexSources);
}

static void VoronoiColoring(MeshType &m, bool frontierFlag=true)
{
  PerVertexPointerHandle sources =  tri::Allocator<MeshType>:: template GetPerVertexAttribute<VertexPointer> (m,"sources");
  assert(tri::Allocator<MeshType>::IsValidHandle(m,sources));

  if(frontierFlag)
  {
    //static_cast<VertexPointer>(NULL) has been introduced just to avoid an error in the MSVS2010's compiler confusing pointer with int. You could use nullptr to avoid it, but it's not supported by all compilers.
    //The error should have been removed from MSVS2012
    std::pair<float,VertexPointer> zz(0.0f,static_cast<VertexPointer>(NULL));
    std::vector< std::pair<float,VertexPointer> > regionArea(m.vert.size(),zz);
    std::vector<VertexPointer> frontierVec;
    GetAreaAndFrontier(m, sources,  regionArea, frontierVec);
    tri::Geodesic<MeshType>::Compute(m,frontierVec);
  }
  tri::UpdateColor<MeshType>::PerVertexQualityRamp(m);
}

static void VoronoiAreaColoring(MeshType &m,std::vector<VertexType *> &seedVec,
                                std::vector< std::pair<float,VertexPointer> > &regionArea)
{
  PerVertexPointerHandle vertexSources =  tri::Allocator<MeshType>:: template GetPerVertexAttribute<VertexPointer> (m,"sources");
  float meshArea = tri::Stat<MeshType>::ComputeMeshArea(m);
  float expectedArea = meshArea/float(seedVec.size());
  for(size_t i=0;i<m.vert.size();++i)
      m.vert[i].C()=Color4b::ColorRamp(expectedArea *0.75f ,expectedArea*1.25f, regionArea[tri::Index(m,vertexSources[i])].first);
}
// It associates the faces with a given vertex according to the vertex associations
//
// It READS  the PerVertex attribute 'sources'
// It WRITES the PerFace attribute 'sources'

static void FaceAssociateRegion(MeshType &m)
{
  PerFacePointerHandle   faceSources =  tri::Allocator<MeshType>:: template GetPerFaceAttribute<VertexPointer> (m,"sources");
  PerVertexPointerHandle vertexSources =  tri::Allocator<MeshType>:: template GetPerVertexAttribute<VertexPointer> (m,"sources");
  for(FaceIterator fi=m.face.begin();fi!=m.face.end();++fi)
  {
    faceSources[fi]=0;
    std::vector<VertexPointer> vp(3);
    for(int i=0;i<3;++i) vp[i]=vertexSources[fi->V(i)];

    for(int i=0;i<3;++i) // First try to associate to the most reached vertex
    {
      if(vp[0]==vp[1] && vp[0]==vp[2]) faceSources[fi] = vp[0];
      else
      {
        if(vp[0]==vp[1] && vp[0]->Q()< vp[2]->Q()) faceSources[fi] = vp[0];
        if(vp[0]==vp[2] && vp[0]->Q()< vp[1]->Q()) faceSources[fi] = vp[0];
        if(vp[1]==vp[2] && vp[1]->Q()< vp[0]->Q()) faceSources[fi] = vp[1];
      }
    }
  }
  tri::UpdateTopology<MeshType>::FaceFace(m);
  int unassCnt=0;
  do
  {
    unassCnt=0;
    for(FaceIterator fi=m.face.begin();fi!=m.face.end();++fi)
    {
      if(faceSources[fi]==0)
      {
        std::vector<VertexPointer> vp(3);
        for(int i=0;i<3;++i)
          vp[i]=faceSources[fi->FFp(i)];

        if(vp[0]!=0 && (vp[0]==vp[1] || vp[0]==vp[2]))
          faceSources[fi] = vp[0];
        else if(vp[1]!=0 && (vp[1]==vp[2]))
          faceSources[fi] = vp[1];
        else
          faceSources[fi] = std::max(vp[0],std::max(vp[1],vp[2]));
        if(faceSources[fi]==0) unassCnt++;
      }
    }
  }
  while(unassCnt>0);
}

// Select all the faces with a given source vertex <vp>
// It reads the PerFace attribute 'sources'

static int FaceSelectAssociateRegion(MeshType &m, VertexPointer vp)
{
  PerFacePointerHandle sources =  tri::Allocator<MeshType>:: template FindPerFaceAttribute<VertexPointer> (m,"sources");
  assert(tri::Allocator<MeshType>::IsValidHandle(m,sources));
  tri::UpdateSelection<MeshType>::Clear(m);
  int selCnt=0;
  for(FaceIterator fi=m.face.begin();fi!=m.face.end();++fi)
  {
    if(sources[fi]==vp)
    {
      fi->SetS();
      ++selCnt;
    }
  }
  return selCnt;
}

// Given a seed <vp>, it selects all the faces that have the minimal distance vertex sourced by the given <vp>.
// <vp> can be null (it search for unreached faces...)
// returns the number of selected faces;
//
// It reads the PerVertex attribute 'sources'
static int FaceSelectRegion(MeshType &m, VertexPointer vp)
{
  PerVertexPointerHandle sources =  tri::Allocator<MeshType>:: template GetPerVertexAttribute<VertexPointer> (m,"sources");
  assert(tri::Allocator<MeshType>::IsValidHandle(m,sources));
  tri::UpdateSelection<MeshType>::Clear(m);
  int selCnt=0;
  for(FaceIterator fi=m.face.begin();fi!=m.face.end();++fi)
  {
    int minInd = 0; float minVal=std::numeric_limits<float>::max();
    for(int i=0;i<3;++i)
    {
      if((*fi).V(i)->Q()<minVal)
      {
        minInd=i;
        minVal=(*fi).V(i)->Q();
      }
    }

    if(	sources[(*fi).V(minInd)] == vp)
    {
      fi->SetS();
      selCnt++;
    }
  }
  return selCnt;
}

/// Given a mesh with for each vertex the link to the closest seed
/// (e.g. for all vertexes we know what is the corresponding voronoi region)
/// we compute:
///  area of all the voronoi regions
///  the vector of the frontier vertexes (e.g. vert of faces shared by at least two regions)
///
///  Area is computed only for triangles that fully belong to a given source.

static void GetAreaAndFrontier(MeshType &m, PerVertexPointerHandle &sources,
                               std::vector< std::pair<float, VertexPointer> > &regionArea, // for each seed we store area
                               std::vector<VertexPointer> &frontierVec)
{
  tri::UpdateFlags<MeshType>::VertexClearV(m);
  frontierVec.clear();
  for(FaceIterator fi=m.face.begin();fi!=m.face.end();++fi)
  {
    VertexPointer s0 = sources[(*fi).V(0)];
    VertexPointer s1 = sources[(*fi).V(1)];
    VertexPointer s2 = sources[(*fi).V(2)];
    assert(s0 && s1 && s2);
    if((s0 != s1) || (s0 != s2) )
    {
      for(int i=0;i<3;++i)
        if(!fi->V(i)->IsV())
        {
          frontierVec.push_back(fi->V(i));
          fi->V(i)->SetV();
        }
    }
    else // the face belongs to a single region; accumulate area;
    {
      if(s0 != 0)
      {
        int seedIndex = tri::Index(m,s0);
        regionArea[seedIndex].first+=DoubleArea(*fi)*0.5f;
        regionArea[seedIndex].second=s0;
      }
    }
  }
}


/// Given a mesh with for each vertex the link to the closest seed
/// we compute:
///  the vector of the corner faces (ie the faces shared exactly by three regions)
///  the vector of the frontier faces that are on the boundary.

static void GetFaceCornerVec(MeshType &m, PerVertexPointerHandle &sources,
                               std::vector<FacePointer> &cornerVec,
                               std::vector<FacePointer> &borderCornerVec)
{
  tri::UpdateFlags<MeshType>::VertexClearV(m);
  cornerVec.clear();
  for(FaceIterator fi=m.face.begin();fi!=m.face.end();++fi)
  {
    VertexPointer s0 = sources[(*fi).V(0)];
    VertexPointer s1 = sources[(*fi).V(1)];
    VertexPointer s2 = sources[(*fi).V(2)];
    assert(s0 && s1 && s2);
    if(s1!=s2 && s0!=s1 && s0!=s2) {
        cornerVec.push_back(&*fi);
      }
    else
    {
      if(isBorderCorner(&*fi,sources))
          borderCornerVec.push_back(&*fi);
    }
  }
}

static bool isBorderCorner(FaceType *f, typename MeshType::template PerVertexAttributeHandle<VertexPointer> &sources)
{
  for(int i=0;i<3;++i)
  {
    if(sources[(*f).V0(i)] != sources[(*f).V1(i)] && f->IsB(i))
      return true;
  }
  return false;
}

static VertexPointer CommonSourceBetweenBorderCorner(FacePointer f0, FacePointer f1,  typename MeshType::template PerVertexAttributeHandle<VertexPointer> &sources)
{
  assert(isBorderCorner(f0,sources));
  assert(isBorderCorner(f1,sources));
  int b0 =-1,b1=-1;
  for(int i=0;i<3;++i)
  {
    if(face::IsBorder(*f0,i)) b0=i;
    if(face::IsBorder(*f1,i)) b1=i;
  }
  assert(b0!=-1 && b1!=-1);

  if( (sources[f0->V0(b0)] == sources[f1->V0(b1)]) || (sources[f0->V0(b0)] == sources[f1->V1(b1)]) )
    return sources[f0->V0(b0)];

  if( (sources[f0->V1(b0)] == sources[f1->V0(b1)]) || (sources[f0->V1(b0)] == sources[f1->V1(b1)]) )
    return sources[f0->V1(b0)];

  assert(0);
  return 0;
}

static void ConvertVoronoiDiagramToMesh(MeshType &m,
                                        MeshType &outMesh, MeshType &outPoly,
                                        std::vector<VertexType *> &seedVec,
                                        DistanceFunctor &df, VoronoiProcessingParameter &vpp )
{
  typename MeshType::template PerVertexAttributeHandle<VertexPointer> sources;
  sources = tri::Allocator<MeshType>:: template GetPerVertexAttribute<VertexPointer> (m,"sources");

  tri::Geodesic<MeshType>::Compute(m,seedVec, df,std::numeric_limits<ScalarType>::max(),0,&sources);
  outMesh.Clear();
  outPoly.Clear();
  tri::UpdateTopology<MeshType>::FaceFace(m);
  tri::UpdateFlags<MeshType>::FaceBorderFromFF(m);

  std::map<VertexPointer, int> seedMap;  // It says if a given vertex of m is a seed (and what)
  for(size_t i=0;i<m.vert.size();++i)
    seedMap[&(m.vert[i])]=-1;
  for(size_t i=0;i<seedVec.size();++i)
    seedMap[seedVec[i]]=i;

  std::vector<FacePointer> innerCornerVec,   // Faces adjacent to three different regions
                           borderCornerVec;  // Faces that are on the border and adjacent to at least two regions.
  GetFaceCornerVec(m, sources, innerCornerVec, borderCornerVec);

  std::map<FacePointer,int> vertexIndCornerMap; // Given a cornerFace (border or inner) what is the corresponding vertex?
  for(size_t i=0;i<m.face.size();++i)
    vertexIndCornerMap[&(m.face[i])]=-1;

  // First add all the needed vertices: seeds and corners
  for(size_t i=0;i<seedVec.size();++i)
    tri::Allocator<MeshType>::AddVertex(outMesh, seedVec[i]->P(),Color4b::White);

  for(size_t i=0;i<innerCornerVec.size();++i){
    tri::Allocator<MeshType>::AddVertex(outMesh, vcg::Barycenter(*(innerCornerVec[i])),Color4b::Gray);
    vertexIndCornerMap[innerCornerVec[i]] = outMesh.vn-1;
  }
  for(size_t i=0;i<borderCornerVec.size();++i){
    Point3f edgeCenter;
    for(int j=0;j<3;++j) if(face::IsBorder(*(borderCornerVec[i]),j))
      edgeCenter=(borderCornerVec[i]->P0(j)+borderCornerVec[i]->P1(j))/2.0f;
    tri::Allocator<MeshType>::AddVertex(outMesh, edgeCenter,Color4b::Gray);
    vertexIndCornerMap[borderCornerVec[i]] = outMesh.vn-1;
  }
  tri::Append<MeshType,MeshType>::MeshCopy(outPoly,outMesh);

  // There is a voronoi edge if there are two corner face that share two sources.
  // In such a case we add a pair of triangles with an edge connecting these two corner faces
  // and with the two involved sources
  // For each pair of adjacent seed we store the first of the two corner that we encounter.
  std::map<std::pair<VertexPointer,VertexPointer>, FacePointer > VoronoiEdge;

  // 1) Build internal triangles
  // Loop build all the triangles connecting seeds with internal corners
  // we loop over the all the voronoi corner (triangles with three different sources)
  // we build
  for(size_t i=0;i<innerCornerVec.size();++i)
  {
    for(int j=0;j<3;++j)
    {
      VertexPointer v0 = sources[innerCornerVec[i]->V0(j)];
      VertexPointer v1 = sources[innerCornerVec[i]->V1(j)];
      if(v1<v0) std::swap(v0,v1); assert(v1!=v0);

      if(VoronoiEdge[std::make_pair(v0,v1)] == 0)
        VoronoiEdge[std::make_pair(v0,v1)] = innerCornerVec[i];
      else
      {
        FacePointer otherCorner = VoronoiEdge[std::make_pair(v0,v1)];
        VertexPointer corner0 = &(outMesh.vert[vertexIndCornerMap[innerCornerVec[i]]]);
        VertexPointer corner1 = &(outMesh.vert[vertexIndCornerMap[otherCorner]]);
        tri::Allocator<MeshType>::AddFace(outMesh,&(outMesh.vert[seedMap[v0]]), corner0, corner1);
        tri::Allocator<MeshType>::AddFace(outMesh,&(outMesh.vert[seedMap[v1]]), corner1, corner0);
      }
    }
  }

  // 2) build the boundary facets:
  // We loop over border corners and build triangles with seed vertex
  // we do **only** triangles with a  bordercorner and a internal 'corner'
  for(size_t i=0;i<borderCornerVec.size();++i)
  {
    VertexPointer s0 = sources[borderCornerVec[i]->V(0)]; // All bordercorner faces have only two different regions
    VertexPointer s1 = sources[borderCornerVec[i]->V(1)];
    if(s1==s0)    s1 = sources[borderCornerVec[i]->V(2)];
    if(s1<s0) std::swap(s0,s1); assert(s1!=s0);

    FacePointer innerCorner = VoronoiEdge[std::make_pair(s0,s1)] ;
    if(innerCorner)
    {
      VertexPointer corner0 = &(outMesh.vert[vertexIndCornerMap[innerCorner]]);
      VertexPointer corner1 = &(outMesh.vert[vertexIndCornerMap[borderCornerVec[i]]]);
      tri::Allocator<MeshType>::AddFace(outMesh,&(outMesh.vert[seedMap[s0]]), corner0, corner1);
      tri::Allocator<MeshType>::AddFace(outMesh,&(outMesh.vert[seedMap[s1]]), corner0, corner1);
    }
  }

  // Final pass
  tri::UpdateFlags<MeshType>::FaceClearV(m);
  bool AllFaceVisited = false;
  while(!AllFaceVisited)
  {
    // search for a unvisited boundary face
    face::Pos<FaceType> pos,startPos;
    AllFaceVisited=true;
    for(size_t i=0; (AllFaceVisited) && (i<borderCornerVec.size()); ++i)
      if(!borderCornerVec[i]->IsV())
      {
        for(int j=0;j<3;++j)
          if(face::IsBorder(*(borderCornerVec[i]),j))
          {
            pos.Set(borderCornerVec[i],j,borderCornerVec[i]->V(j));
            AllFaceVisited =false;
            printf("SearchForBorder\n");
          }
      }
    if(AllFaceVisited) break;
    assert(pos.IsBorder());
    startPos=pos;
    bool foundBorderSeed=false;
    FacePointer curBorderCorner = pos.F();
    do
    {
      pos.F()->SetV();
      pos.NextB();
      if(sources[pos.V()]==pos.V())
        foundBorderSeed=true;
      assert(isBorderCorner(curBorderCorner,sources));
      if(isBorderCorner(pos.F(),sources))
        if(pos.F() != curBorderCorner)
        {
          VertexPointer curReg = CommonSourceBetweenBorderCorner(curBorderCorner, pos.F(),sources);
          VertexPointer curSeed = &(outMesh.vert[seedMap[curReg]]);
          int otherCorner0 = vertexIndCornerMap[pos.F() ];
          int otherCorner1 = vertexIndCornerMap[curBorderCorner];
          VertexPointer corner0 = &(outMesh.vert[otherCorner0]);
          VertexPointer corner1 = &(outMesh.vert[otherCorner1]);
          if(!foundBorderSeed)
            tri::Allocator<MeshType>::AddFace(outMesh,curSeed,corner0,corner1);
          foundBorderSeed=false;
          curBorderCorner=pos.F();
        }
    }
    while(pos!=startPos);
  }

  //**************** CLEANING ***************
  // 1) reorient
  bool oriented,orientable;
  tri::UpdateTopology<MeshType>::FaceFace(outMesh);
  tri::Clean<MeshType>::OrientCoherentlyMesh(outMesh,oriented,orientable);
  assert(orientable);
  // check that the normal of the input mesh are consistent with the result
  tri::UpdateNormal<MeshType>::PerVertexNormalizedPerFaceNormalized(outMesh);
  tri::UpdateNormal<MeshType>::PerVertexNormalizedPerFaceNormalized(m);
  if(seedVec[0]->N() * outMesh.vert[0].N() < 0 )
    tri::Clean<MeshType>::FlipMesh(outMesh);

  tri::UpdateTopology<MeshType>::FaceFace(outMesh);
  tri::UpdateFlags<MeshType>::FaceBorderFromFF(outMesh);

  // 2) Remove Flips
  tri::UpdateNormal<MeshType>::PerFaceNormalized(outMesh);
  tri::UpdateFlags<MeshType>::FaceClearV(outMesh);
  for(FaceIterator fi=outMesh.face.begin();fi!=outMesh.face.end();++fi)
  {
    int badDiedralCnt=0;
    for(int i=0;i<3;++i)
      if(fi->N() * fi->FFp(i)->N() <0 ) badDiedralCnt++;

    if(badDiedralCnt == 2) fi->SetV();
  }
  for(FaceIterator fi=outMesh.face.begin();fi!=outMesh.face.end();++fi)
    if(fi->IsV())  Allocator<MeshType>::DeleteFace(outMesh,*fi);
  tri::Allocator<MeshType>::CompactEveryVector(outMesh);
  tri::UpdateTopology<MeshType>::FaceFace(outMesh);
  tri::UpdateFlags<MeshType>::FaceBorderFromFF(outMesh);
  tri::UpdateFlags<MeshType>::VertexBorderFromFace(outMesh);

  // 3) set up faux bits
  for(FaceIterator fi=outMesh.face.begin();fi!=outMesh.face.end();++fi)
    for(int i=0;i<3;++i)
    {
      size_t v0 = tri::Index(outMesh,fi->V0(i) );
      size_t v1 = tri::Index(outMesh,fi->V1(i) );
      if (v0 < seedVec.size() && !(seedVec[v0]->IsB() && fi->IsB(i))) fi->SetF(i);
      if (v1 < seedVec.size() && !(seedVec[v1]->IsB() && fi->IsB(i))) fi->SetF(i);
    }

  if(vpp.collapseShortEdge)
  {
    float distThr = m.bbox.Diag() * vpp.collapseShortEdgePerc;
    for(FaceIterator fi=outMesh.face.begin();fi!=outMesh.face.end();++fi) if(!fi->IsD())
    {
      for(int i=0;i<3;++i)
        if((Distance(fi->P0(i),fi->P1(i))<distThr) && !fi->IsF(i))
        {
//          printf("Collapsing face %i:%i e%i \n",tri::Index(outMesh,*fi),tri::Index(outMesh,fi->FFp(i)),i);
          if(!fi->V(i)->IsB())
            face::FFEdgeCollapse(outMesh, *fi,i);
          break;
        }
    }
  }

  //******************** END OF CLEANING ****************


  // ******************* star to tri conversion *********
  if(vpp.triangulateRegion)
  {
    for(FaceIterator fi=outMesh.face.begin();fi!=outMesh.face.end();++fi) if(!fi->IsD())
    {
      for(int i=0;i<3;++i)
      {
        bool b0 = fi->V0(i)->IsB();
        bool b1 = fi->V1(i)->IsB();
        if( ((b0  && b1) || (fi->IsF(i) && !b0 && !b1) ) &&
            tri::Index(outMesh,fi->V(i))<seedVec.size())
        {
//          if(b0==b1)
          if(!seedVec[tri::Index(outMesh,fi->V(i))]->IsS())
            face::FFEdgeCollapse(outMesh, *fi,i);
          break;
        }
      }
    }
  }

  // Now a plain conversion of the non faux edges into a polygonal mesh
  std::vector< typename tri::UpdateTopology<MeshType>::PEdge> EdgeVec;
  tri::UpdateTopology<MeshType>::FillUniqueEdgeVector(outMesh,EdgeVec,false);
  tri::UpdateTopology<MeshType>::AllocateEdge(outMesh);

  for(size_t i=0;i<EdgeVec.size();++i)
  {
    size_t e0 = tri::Index(outMesh,EdgeVec[i].v[0]);
    size_t e1 = tri::Index(outMesh,EdgeVec[i].v[1]);
    assert(e0<outPoly.vert.size());
    tri::Allocator<MeshType>::AddEdge(outPoly,&(outPoly.vert[e0]),&(outPoly.vert[e1]));
  }
}

class VoronoiEdge
{
public:
  VertexPointer r0,r1;
  FacePointer f0,f1;
  bool operator == (const VoronoiEdge &ve) const {return ve.r0==r0 && ve.r1==r1; }
  bool operator < (const VoronoiEdge &ve) const { return (ve.r0==r0)?ve.r1<r1:ve.r0<r0; }
  float Len() const { return Distance(vcg::Barycenter(*f0), vcg::Barycenter(*f1)); }
};

static void BuildVoronoiEdgeVec(MeshType &m, std::vector<VoronoiEdge> &edgeVec)
{
  PerVertexPointerHandle sources =  tri::Allocator<MeshType>:: template GetPerVertexAttribute<VertexPointer> (m,"sources");

  edgeVec.clear();
  std::vector<FacePointer> cornerVec;
  std::vector<FacePointer> borderCornerVec;
  GetFaceCornerVec(m,sources,cornerVec,borderCornerVec);
  // Now find all the voronoi edges: each edge (a *face pair) is identified by two voronoi regions
  typedef std::map< std::pair<VertexPointer,VertexPointer>, std::pair<FacePointer,FacePointer> > EdgeMapType;
  EdgeMapType EdgeMap;
  printf("cornerVec.size() %i\n",cornerVec.size());

  for(size_t i=0;i<cornerVec.size();++i)
  {
    for(int j=0;j<3;++j)
    {
      VertexPointer v0 = sources[cornerVec[i]->V0(j)];
      VertexPointer v1 = sources[cornerVec[i]->V1(j)];
      assert(v0!=v1);
      if(v0>v1) std::swap(v1,v0);
      std::pair<VertexPointer,VertexPointer> adjRegion = std::make_pair(v0,v1);
      if(EdgeMap[adjRegion].first==0)
        EdgeMap[adjRegion].first = cornerVec[i];
      else
        EdgeMap[adjRegion].second = cornerVec[i];
    }
  }
  for(size_t i=0;i<borderCornerVec.size();++i)
  {
    VertexPointer v0 = sources[borderCornerVec[i]->V(0)];
    VertexPointer v1 = sources[borderCornerVec[i]->V(1)];
    if(v0==v1) v1 = sources[borderCornerVec[i]->V(2)];
    assert(v0!=v1);
    if(v0>v1) std::swap(v1,v0);
    std::pair<VertexPointer,VertexPointer> adjRegion = std::make_pair(v0,v1);
    if(EdgeMap[adjRegion].first==0)
      EdgeMap[adjRegion].first = borderCornerVec[i];
    else
      EdgeMap[adjRegion].second = borderCornerVec[i];

  }
  typename EdgeMapType::iterator mi;
  for(mi=EdgeMap.begin();mi!=EdgeMap.end();++mi)
  {
    if((*mi).second.first && (*mi).second.second)
    {
      assert((*mi).first.first && (*mi).first.second);
    edgeVec.push_back(VoronoiEdge());
    edgeVec.back().r0 = (*mi).first.first;
    edgeVec.back().r1 = (*mi).first.second;
    edgeVec.back().f0 = (*mi).second.first;
    edgeVec.back().f1 = (*mi).second.second;
    }
  }
}

static void BuildBiasedSeedVec(MeshType &m,
                               DistanceFunctor &df,
                               std::vector<VertexPointer> &seedVec,
                               std::vector<VertexPointer> &frontierVec,
                               std::vector<VertDist> &biasedFrontierVec,
                               VoronoiProcessingParameter &vpp)
{
  biasedFrontierVec.clear();
  if(vpp.unbiasedSeedFlag)
  {
    for(size_t i=0;i<frontierVec.size();++i)
      biasedFrontierVec.push_back(VertDist(frontierVec[i],0));
    assert(biasedFrontierVec.size() == frontierVec.size());
    return;
  }

  std::vector<VoronoiEdge> edgeVec;
  BuildVoronoiEdgeVec(m,edgeVec);
  printf("Found %lu edges on a diagram of %lu seeds\n",edgeVec.size(),seedVec.size());

  std::map<VertexPointer,std::vector<VoronoiEdge *> > SeedToEdgeVecMap;
  std::map< std::pair<VertexPointer,VertexPointer>, VoronoiEdge *> SeedPairToEdgeMap;
  float totalLen=0;
  for(size_t i=0;i<edgeVec.size();++i)
  {
    SeedToEdgeVecMap[edgeVec[i].r0].push_back(&(edgeVec[i]));
    SeedToEdgeVecMap[edgeVec[i].r1].push_back(&(edgeVec[i]));
    SeedPairToEdgeMap[std::make_pair(edgeVec[i].r0, edgeVec[i].r1)]=&(edgeVec[i]);
    assert (edgeVec[i].r0 < edgeVec[i].r1);
    totalLen +=edgeVec[i].Len();
  }

  // compute the perimeter of each region
  std::map <VertexPointer, float> regionPerymeter;
  for(size_t i=0;i<seedVec.size();++i)
  {
    for(size_t j=0;j<SeedToEdgeVecMap[seedVec[i]].size();++j)
    {
      VoronoiEdge *vep = SeedToEdgeVecMap[seedVec[i]][j];
      regionPerymeter[seedVec[i]]+=vep->Len();
    }
    printf("perimeter of region %i is %f\n",i,regionPerymeter[seedVec[i]]);
  }


  PerVertexPointerHandle sources =  tri::Allocator<MeshType>:: template GetPerVertexAttribute<VertexPointer> (m,"sources");
  // The real bias for each edge is (perim)/(edge)
  // each source can belong to two edges max. so the weight is
  std::map<VertexPointer,float> weight;
  std::map<VertexPointer,int> cnt;
  float biasSum = totalLen/5.0f;
  for(FaceIterator fi=m.face.begin();fi!=m.face.end();++fi)
  {
    for(int i=0;i<3;++i)
    {
      VertexPointer s0 = sources[(*fi).V0(i)];
      VertexPointer s1 = sources[(*fi).V1(i)];
      if(s0!=s1)
      {
        if(s0>s1) std::swap(s0,s1);
        VoronoiEdge *ve = SeedPairToEdgeMap[std::make_pair(s0,s1)];
        if(!ve) printf("v %i %i \n",tri::Index(m,s0),tri::Index(m,s1));
        assert(ve);
        float el = ve->Len();
        weight[(*fi).V0(i)] += (regionPerymeter[s0]+biasSum)/(el+biasSum) ;
        weight[(*fi).V1(i)] += (regionPerymeter[s1]+biasSum)/(el+biasSum) ;
        cnt[(*fi).V0(i)]++;
        cnt[(*fi).V1(i)]++;
      }
    }
  }
  for(VertexIterator vi=m.vert.begin();vi!=m.vert.end();++vi)
  {
    if(cnt[&*vi]>0)
    {
//      float bias = weight[&*vi]/float(cnt[&*vi]);
      float bias = weight[&*vi]/float(cnt[&*vi]) + totalLen;
      biasedFrontierVec.push_back(VertDist(&*vi, bias));
    }
  }
  printf("Collected %i frontier vertexes\n",biasedFrontierVec.size());
}


static void DeleteUnreachedRegions(MeshType &m, PerVertexPointerHandle &sources)
{
  tri::UpdateFlags<MeshType>::VertexClearV(m);
  for(size_t i=0;i<m.vert.size();++i)
    if(sources[i]==0) m.vert[i].SetV();

  for(FaceIterator fi=m.face.begin(); fi!=m.face.end();++fi)
    if(fi->V(0)->IsV() || fi->V(1)->IsV() || fi->V(2)->IsV() )
    {
      face::VFDetach(*fi);
      tri::Allocator<MeshType>::DeleteFace(m,*fi);
    }
  //		qDebug("Deleted faces not reached: %i -> %i",int(m.face.size()),m.fn);
  tri::Clean<MeshType>::RemoveUnreferencedVertex(m);
  tri::Allocator<MeshType>::CompactEveryVector(m);
}

/// Let f_p(q) be the squared distance of q from p
/// f_p(q) = (p_x-q_x)^2 + (p_y-q_y)^2 + (p_z-q_z)^2
/// f_p(q) = p_x^2 -2p_xq_x +q_x^2 + ... + p_z^2 -2p_zq_z +q_z^2
///

struct QuadricSumDistance
{
  ScalarType a;
  ScalarType c;
  CoordType b;
  QuadricSumDistance() {a=0; c=0; b[0]=0; b[1]=0; b[2]=0;}
  void AddPoint(CoordType p)
  {
    a+=1;
    assert(c>=0);
    c+=p*p;
    b[0]+= -2.0f*p[0];
    b[1]+= -2.0f*p[1];
    b[2]+= -2.0f*p[2];
  }

  ScalarType Eval(CoordType p) const
  {
    ScalarType d = a*(p*p) + b*p + c;
    assert(d>=0);
    return d;
  }

  CoordType Min() const
  {
    return b * -0.5f;
  }
};

/// Find the new position according to the geodesic rule.
/// For each region, given the frontiers, it chooses the point with the highest distance from the frontier
///
static void QuadricRelax(MeshType &m, std::vector<VertexType *> &seedVec, std::vector<VertexPointer> &frontierVec,
                          std::vector<VertexType *> &newSeeds,
              DistanceFunctor &df, VoronoiProcessingParameter &vpp)
{
  newSeeds.clear();
  typename MeshType::template PerVertexAttributeHandle<VertexPointer> sources;
  sources = tri::Allocator<MeshType>:: template GetPerVertexAttribute<VertexPointer> (m,"sources");
  QuadricSumDistance dz;
  std::vector<QuadricSumDistance> dVec(m.vert.size(),dz);
  assert(m.vert.size()==m.vn);
  for(VertexIterator vi=m.vert.begin();vi!=m.vert.end();++vi)
  {
    assert(sources[vi]!=0);
    int seedIndex = tri::Index(m,sources[vi]);
    dVec[seedIndex].AddPoint(vi->P());
  }
  // Search the local maxima for each region and use them as new seeds
  std::pair<float,VertexPointer> zz(std::numeric_limits<ScalarType>::max(), static_cast<VertexPointer>(0));
  std::vector< std::pair<float,VertexPointer> > seedMaximaVec(m.vert.size(),zz);
  for(VertexIterator vi=m.vert.begin();vi!=m.vert.end();++vi)
  {
    assert(sources[vi]!=0);
    int seedIndex = tri::Index(m,sources[vi]);
    ScalarType val = dVec[seedIndex].Eval(vi->P());
    vi->Q()=val;
    if(seedMaximaVec[seedIndex].first > val)
    {
      seedMaximaVec[seedIndex].first = val;
      seedMaximaVec[seedIndex].second = &*vi;
    }
  }

  tri::UpdateColor<MeshType>::PerVertexQualityRamp(m);
//  tri::io::ExporterPLY<MeshType>::Save(m,"last.ply",tri::io::Mask::IOM_VERTCOLOR + tri::io::Mask::IOM_VERTQUALITY );

  // update the seedvector with the new maxima (For the vertex not selected)
  for(size_t i=0;i<m.vert.size();++i)
    if(seedMaximaVec[i].second) // only  seeds entries have a non zero pointer
    {
      if(vpp.fixSelectedSeed && sources[seedMaximaVec[i].second]->IsS())
        newSeeds.push_back(sources[seedMaximaVec[i].second]);
      else
        newSeeds.push_back(seedMaximaVec[i].second);
    }
}

/// Find the new position according to the geodesic rule.
/// For each region, given the frontiers, it chooses the point with the highest distance from the frontier
///
static void GeodesicRelax(MeshType &m, std::vector<VertexType *> &seedVec, std::vector<VertexPointer> &frontierVec,
                          std::vector<VertexType *> &newSeeds,
              DistanceFunctor &df, VoronoiProcessingParameter &vpp)
{
  newSeeds.clear();
  typename MeshType::template PerVertexAttributeHandle<VertexPointer> sources;
  sources = tri::Allocator<MeshType>:: template GetPerVertexAttribute<VertexPointer> (m,"sources");

  std::vector<typename tri::Geodesic<MeshType>::VertDist> biasedFrontierVec;
  BuildBiasedSeedVec(m,df,seedVec,frontierVec,biasedFrontierVec,vpp);
  tri::Geodesic<MeshType>::Visit(m,biasedFrontierVec,df);
  tri::UpdateColor<MeshType>::PerVertexQualityRamp(m);
  //    tri::io::ExporterPLY<MeshType>::Save(m,"last.ply",tri::io::Mask::IOM_VERTCOLOR + tri::io::Mask::IOM_VERTQUALITY );

  if(vpp.colorStrategy == VoronoiProcessingParameter::DistanceFromBorder)
    tri::UpdateColor<MeshType>::PerVertexQualityRamp(m);

  // Search the local maxima for each region and use them as new seeds
  std::pair<float,VertexPointer> zz(0.0f,static_cast<VertexPointer>(NULL));
  std::vector< std::pair<float,VertexPointer> > seedMaximaVec(m.vert.size(),zz);
  for(VertexIterator vi=m.vert.begin();vi!=m.vert.end();++vi)
  {
    assert(sources[vi]!=0);
    int seedIndex = tri::Index(m,sources[vi]);
    if(seedMaximaVec[seedIndex].first < (*vi).Q())
    {
      seedMaximaVec[seedIndex].first=(*vi).Q();
      seedMaximaVec[seedIndex].second=&*vi;
    }
  }

  // update the seedvector with the new maxima (For the vertex not selected)
  for(size_t i=0;i<seedMaximaVec.size();++i)
    if(seedMaximaVec[i].second)
    {
      if(vpp.fixSelectedSeed && sources[seedMaximaVec[i].second]->IsS())
        newSeeds.push_back(sources[seedMaximaVec[i].second]);
      else
        newSeeds.push_back(seedMaximaVec[i].second);
    }
}

static void PruneSeedByRegionArea(std::vector<VertexType *> &seedVec,
                                  std::vector< std::pair<float,VertexPointer> > &regionArea,
                                  VoronoiProcessingParameter &vpp)
{
  // Smaller area region are discarded
  Distribution<float> H;
  for(size_t i=0;i<regionArea.size();++i)
    if(regionArea[i].second) H.Add(regionArea[i].first);
  float areaThreshold=0;
  if(vpp.areaThresholdPerc != 0) areaThreshold = H.Percentile(vpp.areaThresholdPerc);
  std::vector<VertexType *> newSeedVec;

  // update the seedvector with the new maxima (For the vertex not selected)
  for(size_t i=0;i<seedVec.size();++i)
  {
    if(regionArea[i].first >= areaThreshold)
      newSeedVec.push_back(seedVec[i]);
  }
  swap(seedVec,newSeedVec);
}

/// \brief Perform a Lloyd relaxation cycle over a mesh
///
///

static void VoronoiRelaxing(MeshType &m, std::vector<VertexType *> &seedVec, int relaxIter, DistanceFunctor &df,
                            VoronoiProcessingParameter &vpp, vcg::CallBackPos *cb=0)
{
  tri::RequireVFAdjacency(m);
  tri::UpdateFlags<MeshType>::FaceBorderFromVF(m);
  tri::UpdateFlags<MeshType>::VertexBorderFromFace(m);
  typename MeshType::template PerVertexAttributeHandle<VertexPointer> sources;
  sources = tri::Allocator<MeshType>:: template GetPerVertexAttribute<VertexPointer> (m,"sources");

  for(int iter=0;iter<relaxIter;++iter)
  {
    if(cb) cb(iter*100/relaxIter,"Voronoi Lloyd Relaxation: First Partitioning");

    // first run: find for each point what is the closest to one of the seeds.
    tri::Geodesic<MeshType>::Compute(m, seedVec, df,std::numeric_limits<ScalarType>::max(),0,&sources);

    if(vpp.colorStrategy == VoronoiProcessingParameter::DistanceFromSeed)
      tri::UpdateColor<MeshType>::PerVertexQualityRamp(m);
    // Delete all the (hopefully) small regions that have not been reached by the seeds;

    if(vpp.deleteUnreachedRegionFlag)  DeleteUnreachedRegions(m,sources);
    std::pair<float,VertexPointer> zz(0.0f,static_cast<VertexPointer>(NULL));
    std::vector< std::pair<float,VertexPointer> > regionArea(m.vert.size(),zz);
    std::vector<VertexPointer> frontierVec;

    GetAreaAndFrontier(m, sources,  regionArea, frontierVec);
    assert(frontierVec.size()>0);

    if(vpp.colorStrategy == VoronoiProcessingParameter::RegionArea) VoronoiAreaColoring(m, seedVec, regionArea);

    //    qDebug("We have found %i regions range (%f %f), avg area is %f, Variance is %f 10perc is %f",(int)seedVec.size(),H.Min(),H.Max(),H.Avg(),H.StandardDeviation(),areaThreshold);

    if(cb) cb(iter*100/relaxIter,"Voronoi Lloyd Relaxation: Searching New Seeds");
    std::vector<VertexPointer> newSeedVec;

    if(vpp.geodesicRelaxFlag)
      GeodesicRelax(m,seedVec, frontierVec, newSeedVec, df,vpp);
    else
      QuadricRelax(m,seedVec,frontierVec, newSeedVec, df,vpp);

    assert(newSeedVec.size() == seedVec.size());
    PruneSeedByRegionArea(newSeedVec,regionArea,vpp);

    for(size_t i=0;i<frontierVec.size();++i)
      frontierVec[i]->C() = Color4b::Gray;
    for(size_t i=0;i<seedVec.size();++i)
      seedVec[i]->C() = Color4b::Black;
    for(size_t i=0;i<newSeedVec.size();++i)
      newSeedVec[i]->C() = Color4b::White;

    swap(newSeedVec,seedVec);
  }
}


// Base vertex voronoi coloring algorithm.
// it assumes VF adjacency. No attempt of computing real geodesic distnace is done. Just a BFS visit starting from the seeds.
static void TopologicalVertexColoring(MeshType &m, std::vector<VertexType *> &seedVec)
{
  std::queue<VertexPointer> VQ;

  tri::UpdateQuality<MeshType>::VertexConstant(m,0);

  for(size_t i=0;i<seedVec.size();++i)
  {
    VQ.push(seedVec[i]);
    seedVec[i]->Q()=i+1;
  }

  while(!VQ.empty())
  {
    VertexPointer vp = VQ.front();
    VQ.pop();

    std::vector<VertexPointer> vertStar;
    vcg::face::VVStarVF<FaceType>(vp,vertStar);
    for(typename std::vector<VertexPointer>::iterator vv = vertStar.begin();vv!=vertStar.end();++vv)
    {
      if((*vv)->Q()==0)
      {
        (*vv)->Q()=vp->Q();
        VQ.push(*vv);
      }
    }
  } // end while(!VQ.empty())

}

// Drastic Simplification algorithm.
// Similar in philosopy to the classic grid clustering but using a voronoi partition instead of the regular grid.
//
// This function assumes that in the mOld mesh,  for each vertex you have a quality that denotes the index of the cluster
// mNew is created by collasping onto a single vertex all the vertices that lies in the same cluster.
// Non degenerate triangles are preserved.

static void VoronoiClustering(MeshType &mOld, MeshType &mNew, std::vector<VertexType *> &seedVec)
{
    std::set<Point3i> clusteredFace;

    FaceIterator fi;
    for(fi=mOld.face.begin();fi!=mOld.face.end();++fi)
    {
        if( (fi->V(0)->Q() != fi->V(1)->Q() ) &&
                (fi->V(0)->Q() != fi->V(2)->Q() ) &&
                (fi->V(1)->Q() != fi->V(2)->Q() )  )
                clusteredFace.insert( Point3i(int(fi->V(0)->Q()), int(fi->V(1)->Q()), int(fi->V(2)->Q())));
    }

    tri::Allocator<MeshType>::AddVertices(mNew,seedVec.size());
    for(size_t i=0;i< seedVec.size();++i)
    mNew.vert[i].ImportData(*(seedVec[i]));

    tri::Allocator<MeshType>::AddFaces(mNew,clusteredFace.size());
    std::set<Point3i>::iterator fsi; ;

    for(fi=mNew.face.begin(),fsi=clusteredFace.begin(); fsi!=clusteredFace.end();++fsi,++fi)
    {
        (*fi).V(0) = & mNew.vert[(int)(fsi->V(0)-1)];
        (*fi).V(1) = & mNew.vert[(int)(fsi->V(1)-1)];
        (*fi).V(2) = & mNew.vert[(int)(fsi->V(2)-1)];
    }
}

}; // end class VoronoiProcessing

} // end namespace tri
} // end namespace vcg
#endif
