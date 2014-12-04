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

  typedef SimpleVolume<SimpleVoxel<ScalarType> >                              MyVolume;
  typedef typename vcg::tri::TrivialWalker<MeshType,MyVolume>               MyWalker;
  typedef typename vcg::tri::MarchingCubes<MeshType, MyWalker>              MyMarchingCubes;

  VoronoiVolumeSampling(MeshType &_baseMesh, MeshType &_seedMesh)
    :seedTree(0),surfTree(0),baseMesh(_baseMesh),seedMesh(_seedMesh)
  {

  }

  KdTree<ScalarType>  *seedTree;
  KdTree<ScalarType>  *surfTree;
  typename KdTree<ScalarType>::PriorityQueue pq;
  GridType surfGrid;
  typedef FaceTmark<MeshType> MarkerFace;
  MarkerFace mf;
  vcg::face::PointDistanceBaseFunctor<ScalarType> PDistFunct;

  MeshType &baseMesh;
  MeshType &seedMesh;
  MeshType poissonSurfaceMesh;
  ScalarType poissonRadiusSurface;
  MeshType montecarloVolumeMesh;

  void Init(ScalarType radius=0)
  {
    MeshType montecarloSurfaceMesh;

    if(radius==0) poissonRadiusSurface = baseMesh.bbox.Diag()/50.0f;
    else poissonRadiusSurface = radius;
    ScalarType meshArea = Stat<MeshType>::ComputeMeshArea(baseMesh);
    int MontecarloSampleNum = 10 * meshArea / (radius*radius);
    tri::MeshSampler<MeshType> sampler(montecarloSurfaceMesh);
    tri::SurfaceSampling<MeshType,tri::MeshSampler<CMeshO> >::Montecarlo(baseMesh, sampler, MontecarloSampleNum);
    montecarloSurfaceMesh.bbox = baseMesh.bbox; // we want the same bounding box
    poissonSurfaceMesh.Clear();
    tri::MeshSampler<MeshType> mps(poissonSurfaceMesh);
    typename tri::SurfaceSampling<MeshType,tri::MeshSampler<MeshType> >::PoissonDiskParam pp;
    pp.geodesicDistanceFlag=false;
    tri::SurfaceSampling<MeshType,tri::MeshSampler<MeshType> >::PoissonDiskPruning(mps, montecarloSurfaceMesh, poissonRadiusSurface,pp);
    vcg::tri::UpdateBounding<MeshType>::Box(poissonSurfaceMesh);

    qDebug("Surface Sampling radius %f - montecarlo %ivn - Poisson %ivn",poissonRadiusSurface,montecarloSurfaceMesh.vn,poissonSurfaceMesh.vn);
    VertexConstDataWrapper<MeshType> ww(poissonSurfaceMesh);
    if(surfTree) delete surfTree;
    surfTree = new  KdTree<ScalarType>(ww);

    surfGrid.SetWithRadius(baseMesh.face.begin(),baseMesh.face.end(),poissonRadiusSurface);
    mf.SetMesh(&baseMesh);
  }

// Compute the signed distance from the surface
ScalarType DistanceFromSurface(CoordType &p)
{
  ScalarType squaredDist;
  unsigned int ind;
  surfTree->doQueryClosest(p,ind,squaredDist);
  ScalarType dist = sqrt(squaredDist);
  if( dist > 3.0f*poissonRadiusSurface)
  {
//    CoordType dir = surfTree->getNeighbor(0) - p;
    CoordType dir = this->poissonSurfaceMesh.vert[ind].P() - p;
    const CoordType &surfN = this->poissonSurfaceMesh.vert[ind].N();
    if(dir* surfN > 0) dist= -dist;
    return dist;
  }

  ScalarType _maxDist = this->poissonRadiusSurface*3.0f;
  dist=_maxDist;
  CoordType _closestPt;
  FacePointer f=surfGrid.GetClosest(PDistFunct,mf,p,_maxDist,dist,_closestPt);
  assert(f);
  assert (dist >=0);
  CoordType dir = _closestPt - p;
  if(dir*f->cN() > 0) dist = -dist;

  return dist;
}


ScalarType DistanceFromVoronoiSeed(CoordType p_point)
{
  ScalarType squaredDist;
  unsigned int ind;
  surfTree->doQueryClosest(p_point,ind,squaredDist);
  return math::Sqrt(squaredDist);
}

ScalarType DistanceFromVoronoiFace(CoordType p_point)
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

 ScalarType DistanceFromVoronoiEdge(CoordType p_point)
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

void BarycentricRelaxVoronoiSamples(int relaxStep)
{
  bool changed=false;
  assert(montecarloVolumeMesh.vn > seedMesh.vn*20);
  int i;
  for(i=0;i<relaxStep;++i)
  {
    std::vector<std::pair<int,CoordType> > sumVec(seedMesh.vn,std::make_pair(0,CoordType(0,0,0)));
    for(typename MeshType::VertexIterator vi=montecarloVolumeMesh.vert.begin();vi!=montecarloVolumeMesh.vert.end();++vi)
    {
      unsigned int seedInd;
      ScalarType sqdist;
      seedTree->doQueryClosest(vi->P(),seedInd,sqdist);
      sumVec[seedInd].first++;
      sumVec[seedInd].second+=vi->cP();
    }

    changed=false;
    for(int i=0;i<seedMesh.vert.size();++i)
    {
      if(sumVec[i].first == 0) tri::Allocator<MeshType>::DeleteVertex(seedMesh,seedMesh.vert[i]);
      else
      {
        CoordType prevP = seedMesh.vert[i].P();
        seedMesh.vert[i].P() = sumVec[i].second /ScalarType(sumVec[i].first);
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
  qDebug("performed %i relax step on %i",i,relaxStep);
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

    // Each region has a quadric representing the sum of the squared distances of all the points of its region.
    // First Loop:
    // For each point of the volume add its distance to the quadric of its region.
    for(typename MeshType::VertexIterator vi=montecarloVolumeMesh.vert.begin();vi!=montecarloVolumeMesh.vert.end();++vi)
    {
      unsigned int seedInd;
      ScalarType sqdist;
      seedTree->doQueryClosest(vi->P(),seedInd,sqdist);
      dVec[seedInd].AddPoint(vi->P());
    }

    // Second Loop: Search for each region the point that has minimal squared distance from all other points in that region.
    // We do that evaluating the quadric in each point
    std::vector< std::pair<ScalarType,int> > seedMinimaVec(seedMesh.vert.size(),std::make_pair(std::numeric_limits<ScalarType>::max(),-1 ));

    for(typename MeshType::VertexIterator vi=montecarloVolumeMesh.vert.begin();vi!=montecarloVolumeMesh.vert.end();++vi)
    {
      unsigned int seedInd;
      ScalarType sqdist;
      seedTree->doQueryClosest(vi->P(),seedInd,sqdist);

      ScalarType val = dVec[seedInd].Eval(vi->P());
      if(val < seedMinimaVec[seedInd].first)
      {
        seedMinimaVec[seedInd].first = val;
        seedMinimaVec[seedInd].second = tri::Index(montecarloVolumeMesh,*vi);
      }
    }
    changed=false;
    for(int i=0;i<seedMesh.vert.size();++i)
    {
      CoordType prevP = seedMesh.vert[i].P() ;
      if(seedMinimaVec[i].second == -1) tri::Allocator<MeshType>::DeleteVertex(seedMesh,seedMesh.vert[i]);
      seedMesh.vert[i].P() = montecarloVolumeMesh.vert[seedMinimaVec[i].second].P();
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
  qDebug("performed %i relax step on %i",i,relaxStep);
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
 void BuildScaffoldingMesh(MeshType &scaffoldingMesh, int volumeSide, ScalarType isoThr,int elemEnum, bool surfFlag)
{
   printf("Scaffolding of the mesh \n");
    MyVolume    volume;
    ScalarType       max = math::Max(baseMesh.bbox.DimX(),baseMesh.bbox.DimY(),baseMesh.bbox.DimZ());
    ScalarType       voxel = max / volumeSide;
    int         sizeX = (baseMesh.bbox.DimX() / voxel)+1;
    int         sizeY = (baseMesh.bbox.DimY() / voxel)+1;
    int         sizeZ = (baseMesh.bbox.DimZ() / voxel)+1;

    // Kdtree
//    seedTree->setMaxNofNeighbors(4);

    BoxType bb = BoxType::Construct(baseMesh.bbox);
    bb.Offset(baseMesh.bbox.Diag()*0.04f);
    volume.Init(Point3i(sizeX,sizeY,sizeZ),bb);

    qDebug("Init Volume of %i %i %i",sizeX,sizeY,sizeZ);
   int cnt=0;
   ScalarType offset= volume.voxel.Norm()*isoThr;
    for(ScalarType i=0;i<sizeX;i++)
        for(ScalarType j=0;j<sizeY;j++)
          for(ScalarType k=0;k<sizeZ;k++)
          {
            // check if the point is inside the mesh
            CoordType p;
            volume.IPiToPf(Point3i(i,j,k),p);
            ScalarType surfDist = this->DistanceFromSurface(p);

            ScalarType elemDist;
            switch(elemEnum)
            {
            case 0: elemDist = DistanceFromVoronoiSeed(p) - offset; break;
            case 1: elemDist = DistanceFromVoronoiEdge(p) - offset; break;
            case 2: elemDist = DistanceFromVoronoiFace(p) - offset; break;
            default: assert(0);
            }

            ScalarType val;
            if(surfFlag)
              val = std::max(-elemDist,surfDist);
            else
              val = std::max(elemDist,surfDist);
            volume.Val(i,j,k) = val;
            cnt++;

          }

    // MARCHING CUBES
    qDebug("voxel out %i on %i",cnt,sizeX*sizeY*sizeZ);
    MyWalker    walker;
    MyMarchingCubes	 mc(scaffoldingMesh, walker);
    walker.template BuildMesh <MyMarchingCubes>(scaffoldingMesh, volume, mc,0);
}
 /**
  * @brief
  * start from the montecarlo.
  * Write onto the poisson surface sampling the maximum distance from a vertex inside.
  *
  */
 void ThicknessEvaluator()
 {
//   surfTree->setMaxNofNeighbors(1);
   tri::UpdateQuality<MeshType>::VertexConstant(poissonSurfaceMesh,0);
   for(VertexIterator vi=montecarloVolumeMesh.vert.begin(); vi!=montecarloVolumeMesh.vert.end(); ++vi)
    {
     unsigned int ind;
     ScalarType sqdist;
     this->surfTree->doQueryClosest(vi->P(),ind,sqdist);
     VertexPointer vp = &poissonSurfaceMesh.vert[ind];
     ScalarType dist = math::Sqrt(sqdist);
     if(vp->Q() < dist) vp->Q()=dist;
   }
   tri::UpdateColor<MeshType>::PerVertexQualityRamp(poissonSurfaceMesh);
 }



/*
 * Function: BuildVolumeSampling
 * ----------------------------
 * Build a Poisson-Disk Point cloud that cover all the space of the original mesh m
 *
 */
 void BuildVolumeSampling(int montecarloSampleNum, int seedNum, ScalarType &poissonRadius, vcg::CallBackPos *cb=0)
 {
   montecarloVolumeMesh.Clear();
    math::SubtractiveRingRNG rng;
//    surfTree->setMaxNofNeighbors(1);

    while(montecarloVolumeMesh.vn < montecarloSampleNum)
    {
        CoordType point = math::GeneratePointInBox3Uniform(rng,baseMesh.bbox);
        ScalarType d = this->DistanceFromSurface(point);
        if(d<0){
          vcg::tri::Allocator<MeshType>::AddVertex(montecarloVolumeMesh,point);
          montecarloVolumeMesh.vert.back().Q() = fabs(d);
        }
        if(cb && (montecarloVolumeMesh.vn%1000)==0)
          cb((100*montecarloVolumeMesh.vn)/montecarloSampleNum,"Montecarlo Sampling...");
    }

    vector<VertexPointer>  pruningVec;
    tri::UpdateBounding<MeshType>::Box(montecarloVolumeMesh);
    if(poissonRadius ==0 && seedNum!=0)
      tri::PoissonPruningExact(montecarloVolumeMesh,pruningVec,poissonRadius,seedNum);
    else
      tri::PoissonPruning(montecarloVolumeMesh,pruningVec,poissonRadius,seedNum);

    std::vector<CoordType> seedPts(pruningVec.size());
    for(size_t i=0;i<pruningVec.size();++i)
      seedPts[i]=pruningVec[i]->P();
    tri::Build(this->seedMesh,seedPts);
    // Kdtree must be rebuilt at the end of each step;
    VertexConstDataWrapper<MeshType> vdw(seedMesh);
    if(seedTree) delete seedTree;
    seedTree = new KdTree<ScalarType>(vdw);
}

}; // end class


} // end namespace vcg
} // end namespace vcg
#endif // VORONOI_VOLUME_SAMPLING_H
