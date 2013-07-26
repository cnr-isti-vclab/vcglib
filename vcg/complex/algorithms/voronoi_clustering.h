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
  }
  int colorStrategy;
  float areaThresholdPerc;
  bool deleteUnreachedRegionFlag;
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

static void VoronoiColoring(MeshType &m, std::vector<VertexType *> &seedVec, bool frontierFlag=true)
{
  PerVertexPointerHandle sources =  tri::Allocator<MeshType>:: template GetPerVertexAttribute<VertexPointer> (m,"sources");
  assert(tri::Allocator<MeshType>::IsValidHandle(m,sources));
  tri::Geodesic<MeshType> g;
  VertexPointer farthest;

  if(frontierFlag)
  {
    //static_cast<VertexPointer>(NULL) has been introduced just to avoid an error in the MSVS2010's compiler confusing pointer with int. You could use nullptr to avoid it, but it's not supported by all compilers.
    //The error should have been removed from MSVS2012
    std::pair<float,VertexPointer> zz(0.0f,static_cast<VertexPointer>(NULL));
    std::vector< std::pair<float,VertexPointer> > regionArea(m.vert.size(),zz);
    std::vector<VertexPointer> borderVec;
    std::vector<FacePointer> cornerVec;
    std::vector<FacePointer> borderCornerVec;
    GetAreaAndFrontier(m, sources,  regionArea, borderVec, cornerVec, borderCornerVec);
    tri::Geodesic<MeshType>::Compute(m,borderVec);
  }

  tri::UpdateColor<MeshType>::PerVertexQualityRamp(m);
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

/// Given a mesh with geodesic sources for all vertexes
/// (e.g. for all vertexes we know what is the corresponding voronoi region)
/// we compute Area of all the regions
/// Area is computed only for triangles that fully belong to a given source.

static void GetAreaAndFrontier(MeshType &m, PerVertexPointerHandle &sources,
                               std::vector< std::pair<float,VertexPointer> > &regionArea,
                               std::vector<VertexPointer> &borderVec,
                               std::vector<FacePointer> &cornerVec,
                               std::vector<FacePointer> &borderCornerVec)
{
  tri::UpdateFlags<MeshType>::VertexClearV(m);
  cornerVec.clear();
  borderVec.clear();
  for(FaceIterator fi=m.face.begin();fi!=m.face.end();++fi)
  {
    VertexPointer s0 = sources[(*fi).V(0)];
    VertexPointer s1 = sources[(*fi).V(1)];
    VertexPointer s2 = sources[(*fi).V(2)];
    if((s0 != s1) || (s0 != s2) )
    {
      for(int i=0;i<3;++i)
        borderVec.push_back(fi->V(i));

      if(s1!=s2 && s0!=s1 && s0!=s2) {
        cornerVec.push_back(&*fi);
      }
      else
      {
        for(int i=0;i<3;++i)
        {
          if(sources[(*fi).V0(i)] != sources[(*fi).V1(i)] && fi->IsB(i))
              borderCornerVec.push_back(&*fi);
        }
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


static void ConvertVoronoiDiagramToMesh(MeshType &m, MeshType &outM, MeshType &poly, std::vector<VertexType *> &seedVec,  DistanceFunctor &df, VoronoiProcessingParameter &vpp )
{
  typename MeshType::template PerVertexAttributeHandle<VertexPointer> sources;
  sources = tri::Allocator<MeshType>:: template GetPerVertexAttribute<VertexPointer> (m,"sources");

  tri::Geodesic<MeshType>::Compute(m,seedVec, df,std::numeric_limits<ScalarType>::max(),0,&sources);

  std::map<VertexPointer,int> seedMap;
  for(size_t i=0;i<seedVec.size();++i)
    seedMap[seedVec[i]]=i;

  std::pair<float,VertexPointer> zz(0.0f,VertexPointer(NULL));
  std::vector< std::pair<float,VertexPointer> > regionArea(m.vert.size(),zz);
  std::vector<VertexPointer> borderVec;
  std::vector<FacePointer> cornerVec;
  std::vector<FacePointer> borderCornerVec;
  GetAreaAndFrontier(m, sources,  regionArea, borderVec, cornerVec, borderCornerVec);
  outM.Clear();
  poly.Clear();

  std::map<FacePointer,int> cornerMap;
  for(size_t i=0;i<cornerVec.size();++i)
    cornerMap[cornerVec[i]]=i;

  for(size_t i=0;i<borderCornerVec.size();++i)
    cornerMap[borderCornerVec[i]]=i;

  tri::Allocator<MeshType>::AddVertices(outM,seedVec.size()+cornerVec.size()+borderCornerVec.size());

  for(size_t i=0;i<seedVec.size();++i){
    outM.vert[i].P()=seedVec[i]->P();
    outM.vert[i].C()=Color4b::White;
  }

  int cOff = seedVec.size();
  for(size_t i=0;i<cornerVec.size();++i)
  {
    outM.vert[cOff+i].P()=vcg::Barycenter(*(cornerVec[i]));
    outM.vert[cOff+i].C()=Color4b::Gray;
  }

  int bcOff =seedVec.size()+cornerVec.size();
  for(size_t i=0;i<borderCornerVec.size();++i)
    outM.vert[bcOff+i].P()=vcg::Barycenter(*(borderCornerVec[i]));

  tri::Append<MeshType,MeshType>::MeshCopy(poly,outM);

  // There is a voronoi edge if there are two corner face that share two sources.
  // In such a case we add a pair of triangles with an edge connecting these two corner faces
  // and with the two involved sources
  std::map<std::pair<VertexPointer,VertexPointer>, FacePointer > VoronoiEdge;

  for(size_t i=0;i<cornerVec.size();++i)
  {
    for(int j=0;j<3;++j)
    {
      VertexPointer v0 = sources[cornerVec[i]->V0(j)];
      VertexPointer v1 = sources[cornerVec[i]->V1(j)];
      if(v1<v0) std::swap(v0,v1); assert(v1!=v0);

      if(VoronoiEdge[std::make_pair(v0,v1)] == 0)
        VoronoiEdge[std::make_pair(v0,v1)] = cornerVec[i];
      else
      {
        int otherCorner = cornerMap[VoronoiEdge[std::make_pair(v0,v1)]];
        VertexPointer corner0 = &(outM.vert[cOff+i]);
        VertexPointer corner1 = &(outM.vert[cOff+otherCorner]);
        FaceIterator fi;
        fi = tri::Allocator<MeshType>::AddFace(outM,&(outM.vert[seedMap[v0]]), corner0, corner1);
        fi->SetF(0); fi->SetF(2);
        fi = tri::Allocator<MeshType>::AddFace(outM,&(outM.vert[seedMap[v1]]), corner0, corner1);
        fi->SetF(0); fi->SetF(2);

        tri::Allocator<MeshType>::AddEdge(poly,&(poly.vert[tri::Index(outM,corner0)]),&(poly.vert[tri::Index(outM,corner1)])  );
      }
    }
  }

  // Now build the boundary facets, e.g. the triangles with an edge on the boundary that connects two bordercorner face.
    for(size_t i=0;i<borderCornerVec.size();++i)
    {
        VertexPointer v0 = sources[borderCornerVec[i]->V(0)];
        VertexPointer v1 = sources[borderCornerVec[i]->V(1)];
        if(v1==v0)    v1 = sources[borderCornerVec[i]->V(2)];
        if(v1<v0) std::swap(v0,v1); assert(v1!=v0);

        if(VoronoiEdge[std::make_pair(VertexPointer(0),v0)] == 0)
          VoronoiEdge[std::make_pair(VertexPointer(0),v0)] = borderCornerVec[i];
        else
        {
          int otherCorner = cornerMap[VoronoiEdge[std::make_pair(VertexPointer(0),v0)]];
          VertexPointer corner0 = &(outM.vert[bcOff+i]);
          VertexPointer corner1 = &(outM.vert[bcOff+otherCorner]);
          FaceIterator fi = tri::Allocator<MeshType>::AddFace(outM,&(outM.vert[seedMap[v0]]), corner0, corner1);
          fi->SetF(0);fi->SetF(2);
        }
        if(VoronoiEdge[std::make_pair(VertexPointer(0),v1)] == 0)
          VoronoiEdge[std::make_pair(VertexPointer(0),v1)] = borderCornerVec[i];
        else
        {
          int otherCorner = cornerMap[VoronoiEdge[std::make_pair(VertexPointer(0),v1)]];
          FaceIterator fi=tri::Allocator<MeshType>::AddFaces(outM,1);
          VertexPointer corner0 = &(outM.vert[bcOff+i]);
          VertexPointer corner1 = &(outM.vert[bcOff+otherCorner]);
          fi->V(0) = &(outM.vert[seedMap[v1]]);
          fi->V(1) = corner0;
          fi->V(2) = corner1;
          fi->SetF(0);fi->SetF(2);
        }
        if(VoronoiEdge[std::make_pair(v0,v1)] == 0)
          assert(0);
        else
        {
          int otherCorner = cornerMap[VoronoiEdge[std::make_pair(v0,v1)]];
          FaceIterator fi=tri::Allocator<MeshType>::AddFaces(outM,2);
          VertexPointer corner0 = &(outM.vert[bcOff+i]);
          VertexPointer corner1 = &(outM.vert[cOff+otherCorner]);
          fi->V(0) = &(outM.vert[seedMap[v0]]);
          fi->V(1) = corner0;
          fi->V(2) = corner1;
          fi->SetF(0);fi->SetF(2);
          tri::Allocator<MeshType>::AddEdge(poly,&(poly.vert[tri::Index(outM,corner0)]),&(poly.vert[tri::Index(outM,corner1)])  );

          ++fi;
          fi->V(0) = &(outM.vert[seedMap[v1]]);
          fi->V(1) = corner0;
          fi->V(2) = corner1;
          fi->SetF(0);fi->SetF(2);
          tri::Allocator<MeshType>::AddEdge(poly,&(poly.vert[tri::Index(outM,corner0)]),&(poly.vert[tri::Index(outM,corner1)])  );

        }
    }
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

static void VoronoiRelaxing(MeshType &m, std::vector<VertexType *> &seedVec, int relaxIter, DistanceFunctor &df, VoronoiProcessingParameter &vpp, vcg::CallBackPos *cb=0)
{
  tri::RequireVFAdjacency(m);
  tri::UpdateFlags<MeshType>::FaceBorderFromVF(m);
  typename MeshType::template PerVertexAttributeHandle<VertexPointer> sources;
  sources = tri::Allocator<MeshType>:: template GetPerVertexAttribute<VertexPointer> (m,"sources");

  for(int iter=0;iter<relaxIter;++iter)
  {
    if(cb) cb(iter*100/relaxIter,"Voronoi Lloyd Relaxation: First Partitioning");
    // first run: find for each point what is the closest to one of the seeds.

    tri::Geodesic<MeshType>::Compute(m,seedVec, df,std::numeric_limits<ScalarType>::max(),0,&sources);
    if(vpp.colorStrategy == VoronoiProcessingParameter::DistanceFromSeed)
      tri::UpdateColor<MeshType>::PerVertexQualityRamp(m);
    // Delete all the (hopefully) small regions that have not been reached by the seeds;

    if(vpp.deleteUnreachedRegionFlag)
      DeleteUnreachedRegions(m,sources);
    //static_cast<VertexPointer>(NULL) has been introduced just to avoid an error in the MSVS2010's compiler confusing pointer with int. You could use nullptr to avoid it, but it's not supported by all compilers.
    //The error should have been removed from MSVS2012
    std::pair<float,VertexPointer> zz(0.0f,static_cast<VertexPointer>(NULL));
    std::vector< std::pair<float,VertexPointer> > regionArea(m.vert.size(),zz);
    std::vector<VertexPointer> borderVec;
    std::vector<FacePointer> cornerVec;
    std::vector<FacePointer> borderCornerVec;

    GetAreaAndFrontier(m, sources,  regionArea, borderVec, cornerVec,borderCornerVec);

    // Smaller area region are discarded
    Distribution<float> H;
    for(size_t i=0;i<regionArea.size();++i)
      if(regionArea[i].second) H.Add(regionArea[i].first);

    if(vpp.colorStrategy == VoronoiProcessingParameter::RegionArea)
    {
      float meshArea = tri::Stat<MeshType>::ComputeMeshArea(m);
      float expectedArea = meshArea/float(seedVec.size());
      for(size_t i=0;i<m.vert.size();++i)
          m.vert[i].C()=Color4b::ColorRamp(expectedArea *0.75f ,expectedArea*1.25f, regionArea[tri::Index(m,sources[i])].first);
    }

    float areaThreshold=0;
    if(vpp.areaThresholdPerc != 0) areaThreshold = H.Percentile(vpp.areaThresholdPerc);

//    qDebug("We have found %i regions range (%f %f), avg area is %f, Variance is %f 10perc is %f",(int)seedVec.size(),H.Min(),H.Max(),H.Avg(),H.StandardDeviation(),areaThreshold);

    if(cb) cb(iter*100/relaxIter,"Voronoi Lloyd Relaxation: Searching New Seeds");

    tri::Geodesic<MeshType>::Compute(m,borderVec,df);

    if(vpp.colorStrategy == VoronoiProcessingParameter::DistanceFromBorder)
      tri::UpdateColor<MeshType>::PerVertexQualityRamp(m);

    // Search the local maxima for each region and use them as new seeds
    std::vector< std::pair<float,VertexPointer> > seedMaxima(m.vert.size(),zz);
    for(VertexIterator vi=m.vert.begin();vi!=m.vert.end();++vi)
    {
      assert(sources[vi]!=0);
      int seedIndex = tri::Index(m,sources[vi]);
      if(seedMaxima[seedIndex].first < (*vi).Q())
      {
        seedMaxima[seedIndex].first=(*vi).Q();
        seedMaxima[seedIndex].second=&*vi;
      }
    }
    std::vector<VertexPointer> newSeeds;
    for(size_t i=0;i<seedMaxima.size();++i)
      if(seedMaxima[i].second)
      {
        seedMaxima[i].second->C() = Color4b::Gray;
        if(regionArea[i].first >= areaThreshold)
          newSeeds.push_back(seedMaxima[i].second);
      }


    for(size_t i=0;i<borderVec.size();++i)
      borderVec[i]->C() = Color4b::Gray;

    for(size_t i=0;i<cornerVec.size();++i)
      for(int j=0;j<3;++j)
      cornerVec[i]->V(j)->C() = Color4b::Green;

    for(size_t i=0;i<seedVec.size();++i)
      seedVec[i]->C() = Color4b::Black;

    swap(newSeeds,seedVec);

    for(size_t i=0;i<seedVec.size();++i)
      seedVec[i]->C() = Color4b::White;
  }

//  tri::Allocator<MeshType>::DeletePerVertexAttribute (m,"sources");
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
