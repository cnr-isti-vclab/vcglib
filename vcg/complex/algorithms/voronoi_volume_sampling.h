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
#ifndef __VORONOI_VOLUME_SAMPLING_H
#define __VORONOI_VOLUME_SAMPLING_H
#include <vcg/complex/algorithms/voronoi_processing.h>
#include <vcg/complex/algorithms/create/marching_cubes.h>
#include <vcg/complex/algorithms/create/mc_trivial_walker.h>
#include <vcg/complex/algorithms/point_sampling.h>


namespace vcg
{
namespace tri
{

template< class MeshType>
class VoronoiVolumeSampling
{
public:
  typedef typename tri::VoronoiProcessing<MeshType>::QuadricSumDistance QuadricSumDistance;
  typedef typename MeshType::ScalarType ScalarType;
  typedef typename MeshType::BoxType BoxType;
  typedef typename MeshType::VertexIterator VertexIterator;
  typedef typename MeshType::VertexPointer VertexPointer;
  typedef typename MeshType::CoordType CoordType;
  typedef typename MeshType::FacePointer FacePointer;
  typedef typename vcg::GridStaticPtr<typename MeshType::FaceType, ScalarType> GridType;

  typedef SimpleVolume<SimpleVoxel<ScalarType> >                              VVSVolume;
  typedef typename vcg::tri::TrivialWalker<MeshType,VVSVolume>               VVSWalker;
  typedef typename vcg::tri::MarchingCubes<MeshType, VVSWalker>              VVSMarchingCubes;

  class Param
  {
  public:
    Param()
    {
      elemType=1;
      isoThr=0.1;
      surfFlag=false;
      voxelSide=0;
    }

    int elemType; // the type of element
    ScalarType isoThr;
    ScalarType voxelSide;
    bool surfFlag; 
  };

  VoronoiVolumeSampling(MeshType &_baseMesh)
    :surfTree(0),seedTree(0),baseMesh(_baseMesh),cb(0),restrictedRelaxationFlag(false)
  {
   tri::RequirePerFaceMark(baseMesh);
   tri::UpdateBounding<MeshType>::Box(baseMesh);
   tri::UpdateNormal<MeshType>::PerFaceNormalized(baseMesh);
  }
  
  KdTree<ScalarType>  *surfTree; // used for fast inside query 
  KdTree<ScalarType>  *seedTree; // used to accumulate barycenter in relaxation
  KdTree<ScalarType>  *seedDomainTree; // used to accumulate barycenter in relaxation
  
  typename KdTree<ScalarType>::PriorityQueue pq;
  GridType surfGrid; // used for fast inside query
  typedef FaceTmark<MeshType> MarkerFace;
  MarkerFace mf;
  vcg::face::PointDistanceBaseFunctor<ScalarType> PDistFunct;
      
  MeshType &baseMesh;
  MeshType seedMesh;
  MeshType poissonSurfaceMesh;
  ScalarType poissonRadiusSurface;
  MeshType montecarloVolumeMesh; // we use this mesh as volume evaluator
  MeshType seedDomainMesh;       // where we choose the seeds (by default is the montecarlo volume mesh)
  vcg::CallBackPos *cb;
  math::MarsenneTwisterRNG rng;
  bool restrictedRelaxationFlag;

  
  // Build up the needed structure for efficient point in mesh search. 
  // It uses a poisson disk sampling of the surface plus a
  // kdtree to speed up query point closest on surface for points far from surface. 
  // It initializes the surfGrid, surfTree and poissonSurfaceMesh members   
  void Init(ScalarType _poissonRadiusSurface=0)
  {
    MeshType montecarloSurfaceMesh;
    
    if(_poissonRadiusSurface==0) poissonRadiusSurface = baseMesh.bbox.Diag()/50.0f;
    else poissonRadiusSurface = _poissonRadiusSurface;
    ScalarType meshArea = Stat<MeshType>::ComputeMeshArea(baseMesh);
    int MontecarloSurfSampleNum = 10 * meshArea / (poissonRadiusSurface*poissonRadiusSurface);
    tri::MeshSampler<MeshType> sampler(montecarloSurfaceMesh);
    tri::SurfaceSampling<MeshType,tri::MeshSampler<MeshType> >::SamplingRandomGenerator()=rng;    
    tri::SurfaceSampling<MeshType,tri::MeshSampler<MeshType> >::Montecarlo(baseMesh, sampler, MontecarloSurfSampleNum);
    montecarloSurfaceMesh.bbox = baseMesh.bbox; // we want the same bounding box
    poissonSurfaceMesh.Clear();
    tri::MeshSampler<MeshType> mps(poissonSurfaceMesh);
    typename tri::SurfaceSampling<MeshType,tri::MeshSampler<MeshType> >::PoissonDiskParam pp;
    pp.geodesicDistanceFlag=false;
    
    tri::SurfaceSampling<MeshType,tri::MeshSampler<MeshType> >::PoissonDiskPruning(mps, montecarloSurfaceMesh, poissonRadiusSurface,pp);
    vcg::tri::UpdateBounding<MeshType>::Box(poissonSurfaceMesh);

    printf("Surface Sampling radius %f - montecarlo %ivn - Poisson %ivn\n",poissonRadiusSurface,montecarloSurfaceMesh.vn,poissonSurfaceMesh.vn);
    VertexConstDataWrapper<MeshType> ww(poissonSurfaceMesh);
    if(surfTree) delete surfTree;
    surfTree = new  KdTree<ScalarType>(ww);

    surfGrid.SetWithRadius(baseMesh.face.begin(),baseMesh.face.end(),poissonRadiusSurface);
    mf.SetMesh(&baseMesh);
  }

  // Compute the signed distance from the surface exploting both a kdtree and a ugrid 
  // for a query point p first we use the kdtree with a good poisson sampling of the surface;
  // to get the nearest point on the surface, then if the point is far from the surface we can use the point point distance, while if it is near (e.g. less than 3*poisson radius) we rely on point face distance with a grid.
ScalarType DistanceFromSurface(const CoordType &q, CoordType &closestP)
{
  ScalarType squaredDist;
  unsigned int ind;
  surfTree->doQueryClosest(q,ind,squaredDist);
  ScalarType dist = sqrt(squaredDist);
  if( dist > 3.0f*poissonRadiusSurface)
  {
//    CoordType dir = surfTree->getNeighbor(0) - p;
    CoordType dir = this->poissonSurfaceMesh.vert[ind].P() - q;
    const CoordType &surfN = this->poissonSurfaceMesh.vert[ind].N();
    if(dir* surfN > 0) dist= -dist;
    closestP=this->poissonSurfaceMesh.vert[ind].P();
    return dist;
  }

  ScalarType _maxDist = this->poissonRadiusSurface*3.0f;
  dist=_maxDist;
  FacePointer f=surfGrid.GetClosest(PDistFunct,mf,q,_maxDist,dist,closestP);
  assert(f);
  assert (dist >=0);
  CoordType dir = closestP - q;
  if(dir*f->cN() > 0) dist = -dist;

  return dist;
}


ScalarType DistanceFromVoronoiSeed(const CoordType &p_point)
{
  ScalarType squaredDist;
  unsigned int ind;
  seedTree->doQueryClosest(p_point,ind,squaredDist);
  return math::Sqrt(squaredDist);
}

ScalarType DistanceFromVoronoiFace(const CoordType &p_point)
{

  seedTree->doQueryK(p_point,2,pq);

  std::vector<std::pair<ScalarType, CoordType> > closeSeedVec;
  CoordType p0= this->seedMesh.vert[pq.getIndex(0)].P();
  CoordType p1= this->seedMesh.vert[pq.getIndex(1)].P();
  Plane3<ScalarType> pl; pl.Init((p0+p1)/2.0f,p0-p1);
  return fabs(SignedDistancePlanePoint(pl,p_point));
}

/*
 * Function: scaffolding
 * ----------------------------
 * calculates the distance between the point P and the line R
 * (intersection of the plane P01 P02)
 *
 *   p_point: point to calculate
 *   p_tree: KdTree of the mesh of point
 *   p_m: Mesh of points ( surface and inside )
 *
 *   returns: distance between the point P and the line R
 */

ScalarType DistanceFromVoronoiInternalEdge(const CoordType &p_point)
{
  seedTree->doQueryK(p_point,3,pq);
  std::vector<std::pair<ScalarType, CoordType> > closeSeedVec;
  CoordType p0= this->seedMesh.vert[pq.getIndex(0)].P();
  CoordType p1= this->seedMesh.vert[pq.getIndex(1)].P();
  CoordType p2= this->seedMesh.vert[pq.getIndex(2)].P();

    Plane3<ScalarType>  pl01; pl01.Init((p0+p1)/2.0f,p0-p1);
    Plane3<ScalarType>  pl02; pl02.Init((p0+p2)/2.0f,p0-p2);
    Line3<ScalarType>   voroLine;

    // Calculating the line R that intersect the planes pl01 and pl02
    vcg::IntersectionPlanePlane(pl01,pl02,voroLine);
    // Calculating the distance k between the point p_point and the line R.
    CoordType closestPt;
    ScalarType closestDist;
    vcg::LinePointDistance(voroLine,p_point,closestPt, closestDist);

    return closestDist;
}

ScalarType DistanceFromVoronoiSurfaceEdge(const CoordType &p_point, const CoordType &surfPt)
{
  seedTree->doQueryK(p_point,3,pq);
  pq.sort();
  assert(pq.getWeight(0) <= pq.getWeight(1));
  
  CoordType p0= this->seedMesh.vert[pq.getIndex(0)].P();
  CoordType p1= this->seedMesh.vert[pq.getIndex(1)].P();
  CoordType p2= this->seedMesh.vert[pq.getIndex(2)].P();

  Plane3<ScalarType>  pl01; pl01.Init((p0+p1)/2.0f,p0-p1);
  Plane3<ScalarType>  pl02; pl02.Init((p0+p2)/2.0f,p0-p2);
  Plane3<ScalarType>  pl12; pl12.Init((p1+p2)/2.0f,p1-p2);
  Line3<ScalarType>   voroLine;
  
    
    // Calculating the line R that intersect the planes pl01 and pl02
    vcg::IntersectionPlanePlane(pl01,pl02,voroLine);
    // Calculating the distance k between the point p_point and the line R.
    CoordType closestPt;
    ScalarType voroLineDist;
    vcg::LinePointDistance(voroLine,p_point,closestPt, voroLineDist);

    Plane3<ScalarType> plSurf; plSurf.Init(surfPt, surfPt - p_point);
    Line3<ScalarType>   surfLine;
    // Calculating the line R that intersect the planes pl01 and pl02
    
    ScalarType surfLineDist;
    vcg::IntersectionPlanePlane(pl01,plSurf,surfLine);
    vcg::LinePointDistance(surfLine,p_point,closestPt, surfLineDist);    
    
    return min(voroLineDist,surfLineDist);
}


ScalarType DistanceFromVoronoiCorner(const CoordType &p_point)
{
  seedTree->doQueryK(p_point,4,pq);
  std::vector<std::pair<ScalarType, CoordType> > closeSeedVec;
  CoordType p0= this->seedMesh.vert[pq.getIndex(0)].P();
  CoordType p1= this->seedMesh.vert[pq.getIndex(1)].P();
  CoordType p2= this->seedMesh.vert[pq.getIndex(2)].P();
  CoordType p3= this->seedMesh.vert[pq.getIndex(3)].P();

    Plane3<ScalarType>  pl01; pl01.Init((p0+p1)/2.0f,p0-p1);
    Plane3<ScalarType>  pl02; pl02.Init((p0+p2)/2.0f,p0-p2);
    Plane3<ScalarType>  pl03; pl03.Init((p0+p3)/2.0f,p0-p3);
    Line3<ScalarType>   voroLine;
    
    // Calculating the line R that intersect the planes pl01 and pl02
    vcg::IntersectionPlanePlane(pl01,pl02,voroLine);
    CoordType intersectionPt;
    vcg::IntersectionLinePlane(voroLine,pl03,intersectionPt);
    
    return vcg::Distance(p_point,intersectionPt);
}


void BarycentricRelaxVoronoiSamples(int relaxStep)
{
  bool changed=false;
  assert(montecarloVolumeMesh.vn > seedMesh.vn*20);
  int i;
  for(i=0;i<relaxStep;++i)
  {
    
    std::vector<std::pair<int,CoordType> > sumVec(seedMesh.vn,std::make_pair(0,CoordType(0,0,0)));
    
    // First accumulate for each seed the coord of all the samples that are closest to him.
    for(typename MeshType::VertexIterator vi=montecarloVolumeMesh.vert.begin();vi!=montecarloVolumeMesh.vert.end();++vi)
    {
      unsigned int seedInd;
      ScalarType sqdist;
      seedTree->doQueryClosest(vi->P(),seedInd,sqdist);
      sumVec[seedInd].first++;
      sumVec[seedInd].second+=vi->cP();
    }

    changed=false;
    for(size_t i=0;i<seedMesh.vert.size();++i)
    {
      if(sumVec[i].first == 0) tri::Allocator<MeshType>::DeleteVertex(seedMesh,seedMesh.vert[i]);
      else
      {
        CoordType prevP = seedMesh.vert[i].P();
        seedMesh.vert[i].P() = sumVec[i].second /ScalarType(sumVec[i].first);
        seedMesh.vert[i].Q() = sumVec[i].first;
        if(restrictedRelaxationFlag)
        { 
          unsigned int seedInd;
          ScalarType sqdist;
          seedDomainTree->doQueryClosest(seedMesh.vert[i].P(),seedInd,sqdist);
          seedMesh.vert[i].P() = seedDomainMesh.vert[seedInd].P();
        }        
        if(prevP != seedMesh.vert[i].P()) changed = true;
      }
    }
    tri::Allocator<MeshType>::CompactVertexVector(seedMesh);

    // Kdtree for the seeds must be rebuilt at the end of each step;
    VertexConstDataWrapper<MeshType> vdw(seedMesh);
    delete seedTree;
    seedTree = new KdTree<ScalarType>(vdw);
    if(!changed)
      break;
  }
//  qDebug("performed %i relax step on %i",i,relaxStep);
}

// Given a volumetric sampling of the mesh, and a set of seeds
void QuadricRelaxVoronoiSamples(int relaxStep)
{
  bool changed=false;
  assert(montecarloVolumeMesh.vn > seedMesh.vn*20);
  int i;
  for(i=0;i<relaxStep;++i)
  {
    QuadricSumDistance dz;
    std::vector<QuadricSumDistance> dVec(montecarloVolumeMesh.vert.size(),dz);
    tri::UpdateQuality<MeshType>::VertexConstant(seedMesh,0);

    // Each voronoi region has a quadric representing the sum of the squared distances of all the points of its region.
    // First Loop:
    // For each point of the volume add its distance to the quadric of its region.
    for(typename MeshType::VertexIterator vi=montecarloVolumeMesh.vert.begin();vi!=montecarloVolumeMesh.vert.end();++vi)
    {
      unsigned int seedInd;
      ScalarType sqdist;
      seedTree->doQueryClosest(vi->P(),seedInd,sqdist);
      dVec[seedInd].AddPoint(vi->P());
      seedMesh.vert[seedInd].Q() +=1;
    }

    // Second Loop:  for each region we search in the seed domain the point that has minimal squared distance from all other points in that region.
    // We do that evaluating the quadric in each point
    std::vector< std::pair<ScalarType,int> > seedMinimaVec(seedMesh.vert.size(),std::make_pair(std::numeric_limits<ScalarType>::max(),-1 ));
    for(typename MeshType::VertexIterator vi=seedDomainMesh.vert.begin();vi!=seedDomainMesh.vert.end();++vi)
    {
      unsigned int seedInd;
      ScalarType sqdist;
      seedTree->doQueryClosest(vi->P(),seedInd,sqdist);

      ScalarType val = dVec[seedInd].Eval(vi->P());
      if(val < seedMinimaVec[seedInd].first)
      {
        seedMinimaVec[seedInd].first = val;
        seedMinimaVec[seedInd].second = tri::Index(seedDomainMesh,*vi);
      }
    }
    changed=false;
    for(int i=0;i<seedMesh.vert.size();++i)
    {
      CoordType prevP = seedMesh.vert[i].P() ;
      if(seedMinimaVec[i].second == -1) tri::Allocator<MeshType>::DeleteVertex(seedMesh,seedMesh.vert[i]);
      seedMesh.vert[i].P() = seedDomainMesh.vert[seedMinimaVec[i].second].P();
      if(prevP != seedMesh.vert[i].P()) changed = true;
    }
    tri::Allocator<MeshType>::CompactVertexVector(seedMesh);

    // Kdtree for the seeds must be rebuilt at the end of each step;
    VertexConstDataWrapper<MeshType> vdw(seedMesh);
    delete seedTree;
    seedTree = new KdTree<ScalarType>(vdw);
    if(!changed)
      break;
  }
//  qDebug("performed %i relax step on %i",i,relaxStep);
}


ScalarType ImplicitFunction(const CoordType &p, Param &pp)
{
  CoordType closest;
  ScalarType surfDist = this->DistanceFromSurface(p,closest);
  
  ScalarType elemDist;
  switch(pp.elemType)
  {
  case 0: elemDist = DistanceFromVoronoiSeed(p) - pp.isoThr; break;
  case 1: elemDist = DistanceFromVoronoiSurfaceEdge(p,closest) - pp.isoThr; break;
  case 2: elemDist = DistanceFromVoronoiFace(p) - pp.isoThr; break;
  case 3: elemDist = DistanceFromVoronoiCorner(p) - pp.isoThr; break;
  case 4: elemDist = DistanceFromVoronoiInternalEdge(p) - pp.isoThr; break;
  default: assert(0);
  }
  ScalarType val;
  if(pp.surfFlag)
    val = std::max(-elemDist,surfDist);
  else
    val = std::max(elemDist,surfDist);
  
  return val;
}

/*
 * Function: BuildScaffoldingMesh
 * ----------------------------
 * Build a mesh that is the scaffolding of the original mesh.
 * uses an implicit function and a voronoi3d diagram consisting of the set of inside and
 * surface points of the original mesh m
 *
 *   m: original mesh
 *   surVertex: mesh of surface points
 *   PruningPoisson: mesh of inside and surface points, it's the voronoi3d diagram
 *   n_voxel: number of voxels for the greater side
 */
void BuildScaffoldingMesh(MeshType &scaffoldingMesh, Param &pp)
{
  VVSVolume    volume;
  int         sizeX = (baseMesh.bbox.DimX() / pp.voxelSide)+1;
  int         sizeY = (baseMesh.bbox.DimY() / pp.voxelSide)+1;
  int         sizeZ = (baseMesh.bbox.DimZ() / pp.voxelSide)+1;
  
  int t0=clock();
  BoxType bb = BoxType::Construct(baseMesh.bbox);
  bb.Offset(pp.voxelSide+pp.isoThr*2.0f);
  volume.Init(Point3i(sizeX,sizeY,sizeZ),bb);
  int cnt=0;
  for(ScalarType i=0;i<sizeX;i++)
    for(ScalarType j=0;j<sizeY;j++)
      for(ScalarType k=0;k<sizeZ;k++)
      {
        CoordType p;
        volume.IPiToPf(Point3i(i,j,k),p);
        volume.Val(i,j,k) = ImplicitFunction(p,pp);
        cnt++;
      }
  int t1=clock();
  VVSWalker    walker;
  VVSMarchingCubes	 mc(scaffoldingMesh, walker);
  walker.template BuildMesh <VVSMarchingCubes>(scaffoldingMesh, volume, mc,0);
  int t2=clock();
  printf("Fill Volume (%3i %3i %3i) %5.2f\n", sizeX,sizeY,sizeZ,float(t1-t0)/CLOCKS_PER_SEC);
  printf("Marching %i tris %5.2f\n", scaffoldingMesh.fn,float(t2-t1)/CLOCKS_PER_SEC);
}
 /**
  * @brief Compute an evaulation of the thickness as distance from the medial axis.
  * It starts from a montecarlo volume sampling and try to search for the samples that can be part of the medial axis.
  * It use a sampled representation of the surface. A volume sample is considered part 
  * of the medial axis if there are at least two points that are (almost) the same minimal distance to that point.
  *
  * 
  */
 void ThicknessEvaluator(float distThr, int smoothSize, int smoothIter, MeshType *skelM=0)
 {
   tri::UpdateQuality<MeshType>::VertexConstant(poissonSurfaceMesh,0);
   std::vector<VertexPointer> medialSrc(poissonSurfaceMesh.vert.size(),0);
   for(VertexIterator vi=montecarloVolumeMesh.vert.begin(); vi!=montecarloVolumeMesh.vert.end(); ++vi)
    {
     unsigned int ind;
     ScalarType sqdist;
     this->surfTree->doQueryClosest(vi->P(),ind,sqdist);
     VertexPointer vp = &poissonSurfaceMesh.vert[ind];
     ScalarType minDist = math::Sqrt(sqdist);
     if(vp->Q() < minDist) 
     {
       std::vector<unsigned int> indVec;
       std::vector<ScalarType> sqDistVec;
       
       this->surfTree->doQueryDist( vi->P(), minDist*distThr,indVec,sqDistVec);
       if(indVec.size()>1)
       {
         for(size_t i=0;i<indVec.size();++i)
         {
           VertexPointer vp = &poissonSurfaceMesh.vert[indVec[i]];
           //ScalarType dist = math::Sqrt(sqDistVec[i]);
           if(vp->Q() < minDist) {
             vp->Q()=minDist;
             medialSrc[indVec[i]]=&*vi;
           }             
         }
       }       
     }
   }
   // Now collect the vertexes of the volume mesh that are on the medial surface 
   if(skelM)
   {
     tri::UpdateFlags<MeshType>::VertexClearV(montecarloVolumeMesh);
     for(size_t i=0;i<medialSrc.size();++i)
       medialSrc[i]->SetV();
     for(VertexIterator vi=montecarloVolumeMesh.vert.begin(); vi!=montecarloVolumeMesh.vert.end(); ++vi)
       if(vi->IsV()) tri::Allocator<MeshType>::AddVertex(*skelM,vi->P());
     printf("Generated a medial surf of %i vertexes\n",skelM->vn);
   }
   
   
   tri::Smooth<MeshType>::PointCloudQualityMedian(poissonSurfaceMesh);
   tri::Smooth<MeshType>::PointCloudQualityAverage(poissonSurfaceMesh,smoothSize,smoothIter);
   tri::UpdateColor<MeshType>::PerVertexQualityRamp(poissonSurfaceMesh);
   tri::RedetailSampler<MeshType> rs;
   rs.init(&poissonSurfaceMesh);
   rs.dist_upper_bound = poissonSurfaceMesh.bbox.Diag()*0.05 ;
   rs.qualityFlag = true;
   tri::SurfaceSampling<MeshType, RedetailSampler<MeshType> >::VertexUniform(baseMesh, rs, baseMesh.vn, false);
 }

 void RefineSkeletonVolume(MeshType &skelMesh)
 {
   CoordType closestP;
   int trialNum=0;
   for(int i=0;i<skelMesh.vn;++i)
    {
        CoordType point = math::GeneratePointInBox3Uniform(rng,baseMesh.bbox);
        trialNum++;
        ScalarType d = this->DistanceFromSurface(point, closestP);
        if(d<0){
          vcg::tri::Allocator<MeshType>::AddVertex(montecarloVolumeMesh,point);
          montecarloVolumeMesh.vert.back().Q() = fabs(d);
        }
    }
 }
 

 void BuildMontecarloSampling(int montecarloSampleNum)
 {
   montecarloVolumeMesh.Clear();
   
   int trialNum=0;
   CoordType closest;
   while(montecarloVolumeMesh.vn < montecarloSampleNum)
    {
        CoordType point = math::GeneratePointInBox3Uniform(rng,baseMesh.bbox);
        trialNum++;
        ScalarType d = this->DistanceFromSurface(point,closest);
        if(d<0){
          vcg::tri::Allocator<MeshType>::AddVertex(montecarloVolumeMesh,point);
          montecarloVolumeMesh.vert.back().Q() = fabs(d);
        }
        if(cb && (montecarloVolumeMesh.vn%1000)==0)
          cb((100*montecarloVolumeMesh.vn)/montecarloSampleNum,"Montecarlo Sampling...");
    }
   printf("Made %i Trials to get %i samples\n",trialNum,montecarloSampleNum);
   tri::UpdateBounding<MeshType>::Box(montecarloVolumeMesh);
 }

/*
 * Function: BuildVolumeSampling
 * ----------------------------
 * Build a Poisson-Disk Point cloud that cover all the space of the original mesh m
 *
 */
 void BuildVolumeSampling(int montecarloSampleNum, int poissonSampleNum, ScalarType &poissonRadius, int randSeed)
 {
   if(montecarloSampleNum >0) 
     this->BuildMontecarloSampling(montecarloSampleNum);
   if(seedDomainMesh.vn == 0) 
     tri::Append<MeshType,MeshType>::MeshCopy(seedDomainMesh,montecarloVolumeMesh);     
   
    vector<VertexPointer>  pruningVec;    
    if(poissonRadius ==0 && poissonSampleNum!=0)
      tri::PoissonPruningExact(seedDomainMesh,pruningVec,poissonRadius,poissonSampleNum,0.04,10,randSeed);
    else
      tri::PoissonPruning(seedDomainMesh,pruningVec,poissonRadius,randSeed);

    std::vector<CoordType> seedPts(pruningVec.size());
    for(size_t i=0;i<pruningVec.size();++i)
      seedPts[i]=pruningVec[i]->P();
    tri::BuildMeshFromCoordVector(this->seedMesh,seedPts);
    // Kdtree must be rebuilt at the end of each step;
    VertexConstDataWrapper<MeshType> vdw(seedMesh);
    if(seedTree) delete seedTree;
    seedTree = new KdTree<ScalarType>(vdw);
    
    VertexConstDataWrapper<MeshType> vdw2(seedDomainMesh);
    if(seedDomainTree) delete seedTree;
    seedDomainTree = new KdTree<ScalarType>(vdw);
}

}; // end class


} // end namespace vcg
} // end namespace vcg
#endif // VORONOI_VOLUME_SAMPLING_H
