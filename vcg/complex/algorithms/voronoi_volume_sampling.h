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
  typedef SimpleVolume<SimpleVoxel>                              MyVolume;
  typedef typename vcg::tri::TrivialWalker<MeshType,MyVolume>               MyWalker;
  typedef typename vcg::tri::MarchingCubes<MeshType, MyWalker>              MyMarchingCubes;
  typedef typename vcg::GridStaticPtr<typename MeshType::FaceType> GridType;
  typedef typename MeshType::VertexIterator VertexIterator;
  typedef typename MeshType::VertexPointer VertexPointer;
  typedef typename MeshType::CoordType CoordType;
  typedef typename MeshType::FacePointer FacePointer;

  VoronoiVolumeSampling(MeshType &_baseMesh, MeshType &_seedMesh)
    :baseMesh(_baseMesh),seedMesh(_seedMesh),seedTree(0),surfTree(0)
  {

  }

  KdTree<float>  *seedTree;
  KdTree<float>  *surfTree;
  GridType surfGrid;
  typedef FaceTmark<MeshType> MarkerFace;
  MarkerFace mf;
  vcg::face::PointDistanceBaseFunctor<float> PDistFunct;

  MeshType &baseMesh;
  MeshType &seedMesh;
  MeshType poissonSurfaceMesh;
  float poissonRadiusSurface;
  MeshType montecarloVolumeMesh;

  void Init(float radius=0)
  {
    MeshType montecarloSurfaceMesh;

    if(radius==0) poissonRadiusSurface = baseMesh.bbox.Diag()/50.0f;
    else poissonRadiusSurface = radius;
    float meshArea = Stat<MeshType>::ComputeMeshArea(baseMesh);
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
    surfTree = new  KdTree<float>(ww);
    surfTree->setMaxNofNeighbors(1);

    surfGrid.SetWithRadius(baseMesh.face.begin(),baseMesh.face.end(),poissonRadiusSurface);
    mf.SetMesh(&baseMesh);
  }

// Compute the signed distance from the surface
float DistanceFromSurface(Point3f &p)
{
  surfTree->doQueryK(p);
  float dist = sqrtf(surfTree->getNeighborSquaredDistance(0));
  if( dist > 3.0f*poissonRadiusSurface)
  {
    Point3f dir = surfTree->getNeighbor(0) - p;
    const Point3f &surfN = this->poissonSurfaceMesh.vert[surfTree->getNeighborId(0)].N();
    if(dir* surfN > 0) dist= -dist;
    return dist;
  }

  float _maxDist = this->poissonRadiusSurface*3.0f;
  dist=_maxDist;
  Point3f _closestPt;
  FacePointer f=surfGrid.GetClosest(PDistFunct,mf,p,_maxDist,dist,_closestPt);
  assert(f);
  assert (dist >=0);
  Point3f dir = _closestPt - p;
  if(dir*f->cN() > 0) dist = -dist;

  return dist;
}


float DistanceFromVoronoiSeed(Point3f p_point)
{
  // Calculating the closest point to p_point
  seedTree->doQueryK(p_point);
  float minD = seedTree->getNeighborSquaredDistance(0);
  for(int i=1;i<seedTree->getNofFoundNeighbors();++i)
    minD=std::min(minD,seedTree->getNeighborSquaredDistance(i));
  return sqrtf(minD);
}

float DistanceFromVoronoiFace(Point3f p_point)
{
  seedTree->doQueryK(p_point);

  std::vector<std::pair<float, Point3f> > closeSeedVec;
  for(int i=0;i<seedTree->getNofFoundNeighbors();++i)
     closeSeedVec.push_back(std::make_pair(seedTree->getNeighborSquaredDistance(i),seedTree->getNeighbor(i)));

  std::sort(closeSeedVec.begin(),closeSeedVec.end());
  Point3f p0=closeSeedVec[0].second;
  Point3f p1=closeSeedVec[1].second;
  Plane3f pl; pl.Init((p0+p1)/2.0f,p0-p1);
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

 float DistanceFromVoronoiEdge(Point3f p_point)
{

  seedTree->doQueryK(p_point);
  std::vector<std::pair<float, Point3f> > closeSeedVec;
  for(int i=0;i<seedTree->getNofFoundNeighbors();++i)
     closeSeedVec.push_back(std::make_pair(seedTree->getNeighborSquaredDistance(i),seedTree->getNeighbor(i)));

  std::sort(closeSeedVec.begin(),closeSeedVec.end());
  Point3f p0=closeSeedVec[0].second;
  Point3f p1=closeSeedVec[1].second;
  Point3f p2=closeSeedVec[2].second;


    Plane3f          pl01; pl01.Init((p0+p1)/2.0f,p0-p1);
    Plane3f          pl02; pl02.Init((p0+p2)/2.0f,p0-p2);
    Line3f           voroLine;

    // Calculating the line R that intersect the planes po1 and p02
    vcg::IntersectionPlanePlane(pl01,pl02,voroLine);
    // Calculating the distance k between the point p_point and the line R.
    Point3f closestPt;
    float closestDist;
    vcg::LinePointDistance(voroLine,p_point,closestPt, closestDist);

    return closestDist;
}


 void RelaxVoronoiSamples(int relaxStep)
{
  bool changed=false;
  assert(montecarloVolumeMesh.vn > seedMesh.vn*20);
  int i;
  for(i=0;i<relaxStep;++i)
  {
    seedTree->setMaxNofNeighbors(1);
    QuadricSumDistance dz;
    std::vector<QuadricSumDistance> dVec(montecarloVolumeMesh.vert.size(),dz);

    for(typename MeshType::VertexIterator vi=montecarloVolumeMesh.vert.begin();vi!=montecarloVolumeMesh.vert.end();++vi)
    {
      seedTree->doQueryK(vi->P());
      int seedIndex = seedTree->getNeighborId(0);
      dVec[seedIndex].AddPoint(vi->P());
    }

    // Search the local maxima for each region and use them as new seeds
    std::vector< std::pair<float,int> > seedMaximaVec(seedMesh.vert.size(),std::make_pair(std::numeric_limits<float>::max(),-1 ));

    for(typename MeshType::VertexIterator vi=montecarloVolumeMesh.vert.begin();vi!=montecarloVolumeMesh.vert.end();++vi)
    {
      seedTree->doQueryK(vi->P());
      int seedIndex = seedTree->getNeighborId(0);
      float val = dVec[seedIndex].Eval(vi->P());
      if(val < seedMaximaVec[seedIndex].first)
      {
        seedMaximaVec[seedIndex].first = val;
        seedMaximaVec[seedIndex].second = tri::Index(montecarloVolumeMesh,*vi);
      }
    }
    changed=false;
    for(int i=0;i<seedMesh.vert.size();++i)
    {
      Point3f prevP = seedMesh.vert[i].P() ;
      seedMesh.vert[i].P() = montecarloVolumeMesh.vert[seedMaximaVec[i].second].P();
      if(prevP != seedMesh.vert[i].P()) changed = true;
    }

    // Kdtree must be rebuilt at the end of each step;
    VertexConstDataWrapper<MeshType> vdw(seedMesh);
    delete seedTree;
    seedTree = new KdTree<float>(vdw);
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
 void BuildScaffoldingMesh(MeshType &scaffoldingMesh, int volumeSide, float isoThr,int elemEnum, bool surfFlag)
{
   printf("Scaffolding of the mesh \n");
    MyVolume    volume;
    float       max = math::Max(baseMesh.bbox.DimX(),baseMesh.bbox.DimY(),baseMesh.bbox.DimZ());
    float       voxel = max / volumeSide;
    int         sizeX = (baseMesh.bbox.DimX() / voxel)+1;
    int         sizeY = (baseMesh.bbox.DimY() / voxel)+1;
    int         sizeZ = (baseMesh.bbox.DimZ() / voxel)+1;

    // Kdtree
    seedTree->setMaxNofNeighbors(4);

    volume.bbox=baseMesh.bbox;
    volume.bbox.Offset(baseMesh.bbox.Diag()*0.04f);
    volume.siz = Point3i(sizeX,sizeY,sizeZ);
    volume.ComputeDimAndVoxel();
    volume.Init(Point3i(sizeX,sizeY,sizeZ));

    qDebug("Init Volume of %i %i %i",sizeX,sizeY,sizeZ);
   int cnt=0;
   float offset= volume.voxel.Norm()*isoThr;
    for(float i=0;i<sizeX;i++)
        for(float j=0;j<sizeY;j++)
          for(float k=0;k<sizeZ;k++)
          {
            // check if the point is inside the mesh
            Point3f p;
            volume.IPiToPf(Point3i(i,j,k),p);
            float surfDist = this->DistanceFromSurface(p);

            float elemDist;
            switch(elemEnum)
            {
            case 0: elemDist = DistanceFromVoronoiSeed(p) - offset; break;
            case 1: elemDist = DistanceFromVoronoiEdge(p) - offset; break;
            case 2: elemDist = DistanceFromVoronoiFace(p) - offset; break;
            default: assert(0);
            }

            float val;
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
   surfTree->setMaxNofNeighbors(1);
   tri::UpdateQuality<MeshType>::VertexConstant(poissonSurfaceMesh,0);
   for(VertexIterator vi=montecarloVolumeMesh.vert.begin(); vi!=montecarloVolumeMesh.vert.end(); ++vi)
    {
     this->surfTree->doQueryK(vi->P());
     VertexPointer vp = &poissonSurfaceMesh.vert[surfTree->getNeighborId(0)];
     float dist = sqrt(surfTree->getNeighborSquaredDistance(0));
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
 void BuildVolumeSampling(int montecarloSampleNum, int seedNum, float &poissonRadius, vcg::CallBackPos *cb=0)
 {
   montecarloVolumeMesh.Clear();
    math::SubtractiveRingRNG rng;
    surfTree->setMaxNofNeighbors(1);

    while(montecarloVolumeMesh.vn < montecarloSampleNum)
    {
        Point3f point = math::GeneratePointInBox3Uniform(rng,baseMesh.bbox);
        float d = this->DistanceFromSurface(point);
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
    tri::Build(this->seedMesh,pruningVec);
    // Kdtree must be rebuilt at the end of each step;
    VertexConstDataWrapper<MeshType> vdw(seedMesh);
    if(seedTree) delete seedTree;
    seedTree = new KdTree<float>(vdw);
}

}; // end class


} // end namespace vcg
} // end namespace vcg
#endif // VORONOI_VOLUME_SAMPLING_H
