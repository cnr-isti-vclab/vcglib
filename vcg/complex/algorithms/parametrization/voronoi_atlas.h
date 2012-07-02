#ifndef VORONOI_ATLAS_H
#define VORONOI_ATLAS_H
#include<vcg/complex/algorithms/parametrization/poisson_solver.h>
#include<vcg/complex/algorithms/parametrization/uv_utils.h>
#include<vcg/complex/algorithms/parametrization/distortion.h>
#include<vcg/space/poly_packer.h>
#include<vcg/complex/algorithms/update/texture.h>
#include<vcg/complex/algorithms/point_sampling.h>
#include<vcg/complex/algorithms/voronoi_clustering.h>

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
  class VoroFace    : public Face<  VoroUsedTypes, face::VertexRef, face::BitFlags, face::FFAdj ,face::VFAdj , face::WedgeTexCoord2f> {};
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
      for(int j=0;j<3;++j)
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
      for(int j=0;j<3;++j)
        fi->WT(j).P() = (fi->WT(j).P()-UVBox.min) *ratio;
    }
  }


public:
 // Main parametrization function:
 // it takes a startMesh, copy it and


  static void Build( MeshType &startMesh, MeshType &paraMesh, int sampleNum, bool overlap)
  {
  VoroMesh m;  // the mesh used for the processing is a copy of the passed one.
  tri::Append<VoroMesh, MeshType>::Mesh(m, startMesh);
  tri::Clean<VoroMesh>::RemoveUnreferencedVertex(m);
  tri::Allocator<VoroMesh>::CompactVertexVector(m);
  tri::Allocator<VoroMesh>::CompactFaceVector(m);

  tri::UpdateBounding<VoroMesh>::Box(m);
  std::vector<VoroMesh *> meshRegionVec;
  std::vector< std::vector<Point2f> > uvBorders;
  do
  {
    std::vector<Point3f> PoissonSamples;
    float diskRadius=0;
    tri::PoissonSampling(m,PoissonSamples,sampleNum,diskRadius);
    printf("Sampling created a new mesh of %lu points\n",PoissonSamples.size());
    std::vector<VertexType *> seedVec;
    tri::VoronoiProcessing<VoroMesh>::SeedToVertexConversion(m,PoissonSamples,seedVec);
    tri::UpdateTopology<VoroMesh>::VertexFace(m);
    tri::VoronoiProcessing<VoroMesh>::ComputePerVertexSources(m,seedVec);
    tri::VoronoiProcessing<VoroMesh>::FaceAssociateRegion(m);
    tri::VoronoiProcessing<VoroMesh>::VoronoiColoring(m,seedVec,true);
    tri::io::ExporterPLY<VoroMesh>::Save(m,"dd.ply",tri::io::Mask::IOM_VERTCOLOR);

    std::vector<VoroMesh *> badRegionVec;

    for(size_t i=0; i<seedVec.size();++i)
    {
      VoroMesh *rm = new VoroMesh();
      int selCnt = tri::VoronoiProcessing<VoroMesh>::FaceSelectAssociateRegion(m,seedVec[i]);
      assert(selCnt>0);
      if(overlap){
      tri::UpdateSelection<VoroMesh>::VertexFromFaceLoose(m);
      tri::UpdateSelection<VoroMesh>::FaceFromVertexLoose(m);
      }
      tri::Append<VoroMesh,VoroMesh>::Mesh(*rm, m, true);
      char buf[100]; sprintf(buf,"reg%02i.ply",i);
      tri::io::ExporterPLY<VoroMesh>::Save(*rm,buf,tri::io::Mask::IOM_VERTCOLOR|tri::io::Mask::IOM_WEDGTEXCOORD );

      tri::PoissonSolver<VoroMesh> PS(*rm);
      if(PS.IsFeaseable())
      {
        PS.Init();
        PS.FixDefaultVertices();
        PS.SolvePoisson(false);
        tri::UpdateTexture<VoroMesh>::WedgeTexFromVertexTex(*rm);
        RegularizeTexArea(*rm);

        std::vector<Point2f> uvBorder;
        CollectUVBorder(rm,uvBorder);
        meshRegionVec.push_back(rm);
        uvBorders.push_back(uvBorder);
      } else
      {
        qDebug("ACH - mesh %i is NOT homeomorphic to a disk\n",i);
        badRegionVec.push_back(rm);
      }
    }

    VoroMesh *rm = new VoroMesh();
    tri::VoronoiProcessing<VoroMesh>::FaceSelectAssociateRegion(m,0);
    tri::Append<VoroMesh,VoroMesh>::Mesh(*rm, m, true);

    if(rm->fn>0)
    {
      qDebug("ACH - unreached faces %i fn\n",rm->fn);
      badRegionVec.push_back(rm);
    }
    m.Clear();
    sampleNum = 10;
    if(!badRegionVec.empty())
    {
      for(size_t i=0;i<badRegionVec.size();++i)
        if(badRegionVec[i]->fn>10)
          tri::Append<VoroMesh,VoroMesh>::Mesh(m, *badRegionVec[i], false);

//      tri::io::ExporterPLY<VoroMesh>::Save(m,"buf.ply",tri::io::Mask::IOM_VERTCOLOR|tri::io::Mask::IOM_WEDGTEXCOORD );

      tri::Clean<VoroMesh>::RemoveDuplicateFace(m);
      tri::Clean<VoroMesh>::RemoveUnreferencedVertex(m);
      tri::UpdateNormals<VoroMesh>::PerVertexPerFace(m);
      tri::Allocator<VoroMesh>::CompactVertexVector(m);
      tri::Allocator<VoroMesh>::CompactFaceVector(m);
      qDebug("Still %i faces (from %i regions) to process\n",m.fn,badRegionVec.size());
    }
  } while (m.fn>0);
//  tri::io::ExporterPLY<VoroMesh>::Save(m,"vorocolor.ply",tri::io::Mask::IOM_VERTCOLOR);

  std::vector<Similarity2f> trVec;
  Point2f finalSize;
  PolyPacker<float>::PackAsObjectOrientedRect(uvBorders,Point2f(1024.0f,1024.0f),trVec,finalSize);
  // loop again over all the patches
  for(size_t i=0; i<meshRegionVec.size();++i)
  {
    VoroMesh *rm = meshRegionVec[i];
    for(FaceIterator fi=rm->face.begin();fi!=rm->face.end();++fi)
    {
      for(int j=0;j<3;++j)
      {
        Point2f pp(fi->WT(j).U(),fi->WT(j).V());
        Point2f newpp=trVec[i]*pp;
        fi->WT(j).U()=newpp[0]/1024.0f;
        fi->WT(j).V()=newpp[1]/1024.0f;
      }
    }

    char buf[32]; sprintf(buf,"region_aa_%03i.ply",i);
//    tri::io::ExporterPLY<VoroMesh>::Save(*rm,buf,tri::io::Mask::IOM_VERTCOLOR|tri::io::Mask::IOM_WEDGTEXCOORD );
    tri::Append<MeshType,VoroMesh>::Mesh(paraMesh, *rm, false);
  }
}
};


} // end namespace vcg
} // end namespace tri


#endif // VORONOI_ATLAS_H
