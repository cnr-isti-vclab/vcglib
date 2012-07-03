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

//#include <vcg/simplex/face/topology.h>
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
		typedef typename MeshType::VertexType			VertexType;
		
		ClusteringSampler()
		{
			sampleVec=0;
		}
		ClusteringSampler(std::vector<VertexType *> *_vec)
		{
			sampleVec = _vec;
		}
		std::vector<VertexType *> *sampleVec;
		
		void AddVert(const VertexType &p) 
		{
			sampleVec->push_back((VertexType *)(&p));
		}
	}; // end class ClusteringSampler



template <class MeshType >
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

static void ComputePerVertexSources(MeshType &m, std::vector<VertexType *> &seedVec)
{
  tri::Geo<MeshType> g;
  VertexPointer farthest;
  tri::Allocator<MeshType>::DeletePerVertexAttribute(m,"sources"); // delete any conflicting handle regardless of the type...
  PerVertexPointerHandle vertexSources =  tri::Allocator<MeshType>:: template AddPerVertexAttribute<VertexPointer> (m,"sources");
  tri::Allocator<MeshType>::DeletePerFaceAttribute(m,"sources"); // delete any conflicting handle regardless of the type...
  PerFacePointerHandle faceSources =  tri::Allocator<MeshType>:: template AddPerFaceAttribute<VertexPointer> (m,"sources");
  assert(tri::Allocator<MeshType>::IsValidHandle(m,vertexSources));
  g.FarthestVertex(m,seedVec,farthest,std::numeric_limits<ScalarType>::max(),&vertexSources);
}

static void VoronoiColoring(MeshType &m, std::vector<VertexType *> &seedVec, bool frontierFlag=true)
{
  PerVertexPointerHandle sources =  tri::Allocator<MeshType>:: template GetPerVertexAttribute<VertexPointer> (m,"sources");
  assert(tri::Allocator<MeshType>::IsValidHandle(m,sources));
  tri::Geo<MeshType> g;
  VertexPointer farthest;

		if(frontierFlag)
		{
				std::pair<float,VertexPointer> zz(0,0);
				std::vector< std::pair<float,VertexPointer> > regionArea(m.vert.size(),zz);
				std::vector<VertexPointer> borderVec;
				GetAreaAndFrontier(m, sources,  regionArea, borderVec);
				g.FarthestVertex(m,borderVec,farthest);
		}

		tri::UpdateColor<MeshType>::VertexQualityRamp(m);
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
  PerFacePointerHandle sources =  tri::Allocator<MeshType>:: template GetPerFaceAttribute<VertexPointer> (m,"sources");
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

// find the vertexes of frontier faces
// and compute Area of all the regions
static void GetAreaAndFrontier(MeshType &m, PerVertexPointerHandle &sources,
		std::vector< std::pair<float,VertexPointer> > &regionArea,
		std::vector<VertexPointer> &borderVec)
{
		tri::UpdateFlags<MeshType>::VertexClearV(m);
		for(FaceIterator fi=m.face.begin();fi!=m.face.end();++fi)
		{
			if(	sources[(*fi).V(0)] != sources[(*fi).V(1)]  ||
					sources[(*fi).V(0)] != sources[(*fi).V(2)] )
						{
							for(int i=0;i<3;++i)
								{
									(*fi).V(i)->SetV();	
									(*fi).V(i)->C() = Color4b::Black;
								}
						}
				else // the face belongs to a single region; accumulate area;
				{
					if(sources[(*fi).V(0)] != 0)
					{
					   int seedIndex = sources[(*fi).V(0)] - &*m.vert.begin();
					   regionArea[seedIndex].first+=DoubleArea(*fi);
					   regionArea[seedIndex].second=sources[(*fi).V(0)];
					}
				}
		}
				
		// Collect the frontier vertexes
		borderVec.clear();
		for(VertexIterator vi=m.vert.begin();vi!=m.vert.end();++vi) 
			if((*vi).IsV()) borderVec.push_back(&*vi);
}

static void VoronoiRelaxing(MeshType &m, std::vector<VertexType *> &seedVec, int relaxIter, int /*percentileClamping*/, vcg::CallBackPos *cb=0)
{			
	for(int iter=0;iter<relaxIter;++iter)
	{
		if(cb) cb(iter*100/relaxIter,"Voronoi Lloyd Relaxation: First Partitioning");
		tri::Geo<MeshType> g;
    VertexPointer farthest;
		// first run: find for each point what is the closest to one of the seeds.
		typename MeshType::template PerVertexAttributeHandle<VertexPointer> sources;
		sources = tri::Allocator<MeshType>:: template AddPerVertexAttribute<VertexPointer> (m,"sources");
		
    g.FarthestVertex(m,seedVec,farthest,std::numeric_limits<ScalarType>::max(),&sources);
		
		std::pair<float,VertexPointer> zz(0,0);
		std::vector< std::pair<float,VertexPointer> > regionArea(m.vert.size(),zz);
		std::vector<VertexPointer> borderVec;

		GetAreaAndFrontier(m, sources,  regionArea, borderVec);
		
		// Smaller area region are discarded
		Distribution<float> H;
    for(size_t i=0;i<regionArea.size();++i)
			if(regionArea[i].second) H.Add(regionArea[i].first);
			
		float areaThreshold;
		if(iter==0) areaThreshold = H.Percentile(.1f);
		else areaThreshold = H.Percentile(.001f);
	//qDebug("We have found %i regions range (%f %f), avg area is %f, Variance is %f 10perc is %f",(int)seedVec.size(),H.Min(),H.Max(),H.Avg(),H.StandardDeviation(),areaThreshold);
  
		if(cb) cb(iter*100/relaxIter,"Voronoi Lloyd Relaxation: Searching New Seeds");
			
    g.FarthestVertex(m,borderVec,farthest);
        tri::UpdateColor<MeshType>::VertexQualityRamp(m);

		// Search the local maxima for each region and use them as new seeds	
		std::vector< std::pair<float,VertexPointer> > seedMaxima(m.vert.size(),zz);
		for(VertexIterator vi=m.vert.begin();vi!=m.vert.end();++vi) 
		{
			int seedIndex = sources[vi] - &*m.vert.begin();
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
		
		tri::UpdateColor<MeshType>::VertexQualityRamp(m);
    for(size_t i=0;i<seedVec.size();++i)
			seedVec[i]->C() = Color4b::Black;
		
    for(size_t i=0;i<borderVec.size();++i)
			borderVec[i]->C() = Color4b::Gray;
		
		swap(newSeeds,seedVec);
		
    for(size_t i=0;i<seedVec.size();++i)
			seedVec[i]->C() = Color4b::White;

		tri::Allocator<MeshType>::DeletePerVertexAttribute (m,"sources");
		
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
	mNew.vert[i].ImportLocal(*(seedVec[i]));																										

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
