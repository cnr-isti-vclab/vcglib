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
#ifndef VORONOI_ATLAS_H
#define VORONOI_ATLAS_H

#include<vcg/complex/algorithms/parametrization/poisson_solver.h>
#include<vcg/complex/algorithms/parametrization/uv_utils.h>
#include<vcg/complex/algorithms/parametrization/distortion.h>
#include<vcg/space/outline2_packer.h>
#include<vcg/space/rasterized_outline2_packer.h>
#include<vcg/complex/algorithms/update/texture.h>
#include<vcg/complex/algorithms/point_sampling.h>
#include<vcg/complex/algorithms/voronoi_processing.h>

//#include<wrap/qt/outline2_rasterizer.h>

namespace vcg {
namespace tri {

template <class MeshType>
class VoronoiAtlas
{
//private:
public:
  class VoroEdge;
  class VoroFace;
  class VoroVertex;
  struct VoroUsedTypes : public UsedTypes<	Use<VoroVertex>   ::template AsVertexType,
                                          Use<VoroEdge>     ::template AsEdgeType,
                                          Use<VoroFace>     ::template AsFaceType>{};

  class VoroVertex  : public Vertex< VoroUsedTypes, vertex::Coord3f, vertex::Normal3f, vertex::TexCoord2f, vertex::VFAdj , vertex::Qualityf, vertex::Color4b, vertex::BitFlags  >{};
  class VoroFace    : public Face<  VoroUsedTypes, face::VertexRef, face::BitFlags, face::FFAdj ,face::VFAdj , face::CurvatureDirf,face::WedgeTexCoord2f> {};
  class VoroEdge    : public Edge< VoroUsedTypes>{};
  class VoroMesh    : public tri::TriMesh< std::vector<VoroVertex>, std::vector<VoroFace> , std::vector<VoroEdge>  > {};

  typedef typename VoroMesh::FaceIterator FaceIterator;
  typedef typename VoroMesh::VertexType VertexType;
  typedef typename VoroMesh::FaceType FaceType;

  static void CollectUVBorder(VoroMesh *rm, std::vector<Point2f> &uvBorder)
  {
    tri::UpdateTopology<VoroMesh>::FaceFace(*rm);
    tri::UpdateFlags<VoroMesh>::FaceClearV(*rm);
    for(FaceIterator fi=rm->face.begin();fi!=rm->face.end();++fi)
    {
      for(short j=0;j<3;++j)
        if(face::IsBorder(*fi,j) && !(fi->IsV()))
        {
          face::Pos<FaceType> pp(&*fi,j,fi->V(j));
          assert(pp.IsBorder());
          face::Pos<FaceType> startPos = pp;
          do
          {
            uvBorder.push_back( pp.F()->WT(pp.VInd()).P() );
            pp.F()->SetV();
            pp.NextB();
          } while(pp != startPos);
        }
    }
  }

 // take a mesh and rescale its uv so that they are in the 0..1 range
 static void RegularizeTexArea(VoroMesh &m)
  {
    float areaTex=0;
    float areaGeo=0;

    vcg::Box2f UVBox = tri::UV_Utils<VoroMesh>::PerWedgeUVBox(m);
    for(FaceIterator fi=m.face.begin();fi!=m.face.end();++fi)
    {
      areaTex+= fabs((fi->WT(1).P() - fi->WT(0).P()) ^ (fi->WT(2).P() - fi->WT(0).P())) ;
      areaGeo+= DoubleArea(*fi);
    }

    float ratio = sqrt(areaGeo/areaTex);

    for(FaceIterator fi=m.face.begin();fi!=m.face.end();++fi)
    {
      for(short j=0;j<3;++j)
        fi->WT(j).P() = (fi->WT(j).P()-UVBox.min) *ratio;
    }
  }


public:
 struct VoronoiAtlasParam
 {
   VoronoiAtlasParam()
   {
     maxIterNum = 5;
     sampleNum=10;
     overlap=false;
   }

   struct Stat
   {
     void clear() { iterNum=totalTime=unwrapTime=voronoiTime=samplingTime=0;}
     int totalTime;
     int unwrapTime;
     int voronoiTime;
     int samplingTime;

     int regionNum;
     int iterNum;
   };

   int sampleNum;
   bool overlap;
   Stat vas;
   int maxIterNum;
   CallBackPos *cb=vcg::CErrCallBackPos;
 };

 // Main parametrization function:
 // it takes a startMesh, copy it and

 //static void Build( MeshType &startMesh, MeshType &paraMesh, VoronoiAtlasParam &pp,short fccoef,short fnLim)
  static void Build( MeshType &startMesh, MeshType &paraMesh, VoronoiAtlasParam &pp)
  {

    pp.vas.clear();

    int t0=clock();
    short fccoef = 10,short fnLim = 50;

  VoroMesh m;  // the mesh used for the processing is a copy of the passed one.
  tri::Append<VoroMesh, MeshType>::Mesh(m, startMesh);
  tri::Clean<VoroMesh>::RemoveUnreferencedVertex(m);
  tri::Allocator<VoroMesh>::CompactVertexVector(m);
  tri::Allocator<VoroMesh>::CompactFaceVector(m);

  tri::UpdateBounding<VoroMesh>::Box(m);
  std::vector<VoroMesh *> meshRegionVec;
  std::vector< std::vector<Point2f> > uvBorders;

  std::vector<Point3f> PoissonSamples;
  int st0,st1,st2,tp0,tp1,selCnt,foldedCnt;
  float diskRadius;
  std::vector<VertexType *> seedVec;
  std::vector<VoroMesh *> badRegionVec;
  std::vector<Point2f> uvBorder;
  //std::unique_ptr<VoroMesh> rm;
  //VoroMesh *rm;

  // Main processing loop
  do
  {
    printf("ITERATION %i sampling mesh of %i with %i *\n",pp.vas.iterNum,m.fn,pp.sampleNum);

    st0=clock();

    PoissonSamples.clear();
    diskRadius=0;
    tri::PoissonSampling(m,PoissonSamples,pp.sampleNum,diskRadius);

    st1=clock();

    pp.vas.samplingTime+= st1-st0;
    pp.cb(50,StrFormat("Sampling created a new mesh of %lu points\n",PoissonSamples.size()));

    printf("Sampling created a new mesh of %lu points\n",PoissonSamples.size());

    EuclideanDistance<VoroMesh> edFunc;

    seedVec.clear();

    tri::VoronoiProcessing<VoroMesh>::SeedToVertexConversion(m,PoissonSamples,seedVec);
    tri::UpdateTopology<VoroMesh>::VertexFace(m);
    tri::VoronoiProcessing<VoroMesh>::ComputePerVertexSources(m,seedVec,edFunc);
    tri::VoronoiProcessing<VoroMesh>::FaceAssociateRegion(m);
    tri::VoronoiProcessing<VoroMesh>::VoronoiColoring(m,true);

    badRegionVec.clear();

    st2=clock();

    pp.vas.voronoiTime+=st2-st1;

    printf("Voronoi prepr.1 time = %i\n",pp.vas.voronoiTime);
    printf("seedVec.size() = %i\n",seedVec.size());

    for(size_t i=0; i<seedVec.size();++i)
    {

      //rm = new VoroMesh();
      //rm.reset(new VoroMesh());
      VoroMesh *rm = new VoroMesh();

      selCnt = tri::VoronoiProcessing<VoroMesh>::FaceSelectAssociateRegion(m,seedVec[i]);

      pp.cb(50,StrFormat("Region %i of %i faces\n",i,selCnt));

      printf("Region %i of %i faces\n",i,selCnt);

      if(selCnt==0) continue;
      assert(selCnt>0);

      if(pp.overlap){
      tri::UpdateSelection<VoroMesh>::VertexFromFaceLoose(m);
      tri::UpdateSelection<VoroMesh>::FaceFromVertexLoose(m);
      }

      tri::Append<VoroMesh,VoroMesh>::Mesh(*rm, m, true);
      tp0=clock();
      tri::PoissonSolver<VoroMesh> PS(*rm);
      tri::UpdateBounding<VoroMesh>::Box(*rm);

      if(PS.IsFeasible())
      {
        PS.Init();
        PS.FixDefaultVertices();
        PS.SolvePoisson(false);
        tri::UpdateTexture<VoroMesh>::WedgeTexFromVertexTex(*rm);
        RegularizeTexArea(*rm);

        uvBorder.clear();
        //CollectUVBorder(rm.get(),uvBorder);
        //meshRegionVec.push_back(rm.get());
        CollectUVBorder(rm,uvBorder);
        meshRegionVec.push_back(rm);
        uvBorders.push_back(uvBorder);
        foldedCnt = tri::Distortion<VoroMesh,false>::FoldedNum(*rm);

        if( foldedCnt > rm->fn/fccoef)
        {
            //badRegionVec.push_back(rm.get());
            badRegionVec.push_back(rm);
          printf("-- region %i Parametrized but with %i fold on %i\n",i,foldedCnt,rm->fn);
        }
        else printf("-- region %i Parametrized\n",i);
      }
      else
      {
        printf("-- region %i is NOT homeomorphic to a disk\n",i);
        //badRegionVec.push_back(rm.get());
        badRegionVec.push_back(rm);
      }

      tp1=clock();
      pp.vas.unwrapTime +=tp1-tp0;
      ++pp.vas.iterNum;

      printf("unwrapTime = %i,",pp.vas.unwrapTime);
      //printf("1. rm = %i\n",rm.get());
      //printf("1. rm = %i\n",rm);
      //Why if next line is uncommented next subsequent operation fails with bad vector length:  tri::Append<VoroMesh,VoroMesh>::Mesh(m, *badRegionVec[i], false);
      //delete rm;

    }

    printf("\n -- Completed (%i bad regions) -- \n", badRegionVec.size());

    //rm = new VoroMesh();
    VoroMesh *rm = new VoroMesh();
    //rm.reset(new VoroMesh());

    tri::VoronoiProcessing<VoroMesh>::FaceSelectAssociateRegion(m,0);

    tri::Append<VoroMesh,VoroMesh>::Mesh(*rm, m, true);

    if(rm->fn>0)
    {
      printf("ACH - unreached faces %i fn\n",rm->fn);
      //badRegionVec.push_back(rm.get());
      badRegionVec.push_back(rm);
    }

    m.Clear();
    pp.sampleNum = fccoef;

    printf("badRegionVec.size() = %i\n",badRegionVec.size());

    if(!badRegionVec.empty())
    {
      for(size_t i=0;i<badRegionVec.size();++i)
        {
          printf("badRegionVec[i]->fn = %i\n",badRegionVec[i]->fn);
        if(badRegionVec[i]->fn>fnLim)
          {
            tri::Append<VoroMesh,VoroMesh>::Mesh(m, *badRegionVec[i], false);
          }
      }

      tri::Clean<VoroMesh>::RemoveDuplicateFace(m);
      tri::Clean<VoroMesh>::RemoveUnreferencedVertex(m);
      tri::Allocator<VoroMesh>::CompactVertexVector(m);
      tri::Allocator<VoroMesh>::CompactFaceVector(m);
    }

    //printf("rm = %i\n",rm.get());
    printf("rm = %i\n",rm);

    delete rm;
  } while (m.fn>0);

  std::vector<Similarity2f> trVec;
  Point2f finalSize;
  //PolyPacker<float>::WritePolyVec(uvBorders,"borders.poly");
  PolyPacker<float>::PackAsObjectOrientedRect(uvBorders,Point2i(1024,1024),trVec,finalSize);
//  RasterizedOutline2Packer<float,QtOutline2Rasterizer>::Parameters prp;
//  RasterizedOutline2Packer<float,QtOutline2Rasterizer>::Pack(uvBorders,Point2i(1024,1024),trVec,prp);

  // loop again over all the patches
  pp.vas.regionNum=meshRegionVec.size();
  printf("meshRegionVec.size() = %i\n",meshRegionVec.size());

  //VoroMesh *rm;
  for(size_t i=0; i<meshRegionVec.size();++i)
  {
    //rm.reset(meshRegionVec[i]);

    printf("rm->face.size() = %i||",meshRegionVec[i]->face.size());

    for(FaceIterator fi=meshRegionVec[i]->face.begin();fi!=meshRegionVec[i]->face.end();++fi)
    {
      for(short j=0;j<3;++j)
      {
        Point2f pp(fi->WT(j).U(),fi->WT(j).V());
        Point2f newpp=trVec[i]*pp;
        fi->WT(j).U()=newpp[0]/1024.0f;
        fi->WT(j).V()=newpp[1]/1024.0f;
      }
    }
    //tri::Append<MeshType,VoroMesh>::Mesh(paraMesh, *rm, false);
    tri::Append<MeshType,VoroMesh>::Mesh(paraMesh, *(meshRegionVec[i]), false);

  }
  int t2=clock();
  pp.vas.totalTime=t2-t0;

  printf("\n pp.vas.totalTime = %i\n",pp.vas.totalTime);
}
}; //end


} // end namespace vcg
} // end namespace tri


#endif // VORONOI_ATLAS_H

