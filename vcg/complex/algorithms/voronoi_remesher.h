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

#ifndef _VCGLIB_VORONOI_REMESHER_H
#define _VCGLIB_VORONOI_REMESHER_H

#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/update/topology.h>
#include <vcg/complex/algorithms/refine.h>
#include <vcg/complex/algorithms/clean.h>
#include <vcg/complex/algorithms/voronoi_processing.h>
#include <vcg/complex/algorithms/point_sampling.h>
#include <vcg/complex/algorithms/crease_cut.h>
#include <vcg/complex/algorithms/curve_on_manifold.h>

#include <memory>
#include <string>
#include <vector>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <cmath>
#include <array>
#include <utility>

//#define DEBUG_VORO 1
//#include <wrap/io_trimesh/export.h>
//#include <QString>

namespace vcg {
namespace tri {

class VoroEdgeMeshAux
{
    class EmEdgeType;
    class EmVertexType;
    class EUsedTypes : public vcg::UsedTypes<vcg::Use<EmVertexType>::AsVertexType,
                                             vcg::Use<EmEdgeType>::AsEdgeType> {};
    class EmVertexType : public vcg::Vertex<EUsedTypes
            , vcg::vertex::Coord3d
            , vcg::vertex::BitFlags
            , vcg::vertex::VEAdj> {};
    class EmEdgeType   : public vcg::Edge<EUsedTypes
            , vcg::edge::VertexRef
            , vcg::edge::BitFlags
            , vcg::edge::EEAdj
            , vcg::edge::VEAdj> {};
public:
    class EdgeMeshType : public vcg::tri::TriMesh<std::vector<EmVertexType>, std::vector<EmEdgeType> >
    {
    public:
      ~EdgeMeshType()
      {
        this->Clear();
        this->ClearAttributes();
      }
    };
};

template <class MeshType>
class Remesher
{
public:
  typedef Remesher                  ThisType;
  
  typedef MeshType                      Mesh;
  typedef typename Mesh::ScalarType     ScalarType;
  typedef typename Mesh::CoordType      CoordType;
  typedef typename Mesh::FaceType       FaceType;
  typedef typename Mesh::FacePointer    FacePointer;
  typedef typename Mesh::VertexType     VertexType;
  typedef typename Mesh::VertexPointer  VertexPointer;
  typedef typename Mesh::FaceIterator   FaceIterator;
  typedef typename Mesh::VertexIterator VertexIterator;
  
  typedef std::shared_ptr<Mesh>         MeshPtr;
  
protected:
  typedef face::Pos<FaceType>                            PosType;

  typedef typename VoroEdgeMeshAux::EdgeMeshType EdgeMeshType;

  /// \brief splitCC split the provided mesh into connected components.
  /// \param mesh the inputMesh.
  /// \return the vector of connected components (meshes) for the input model
  /// (if the input mesh is a single connected component returns an empty vector).
  ///
  inline static std::vector<MeshPtr> splitCC(MeshType & mesh)
  {
    std::vector<MeshPtr> ret;
    
    // find the connected components
    std::vector<std::pair<int, typename MeshType::FacePointer> > CCV;
    Clean<MeshType>::ConnectedComponents(mesh, CCV);
    
    if (CCV.size() == 1)
      return ret;
    for(size_t i=0; i<CCV.size(); ++i)
    {
      UpdateSelection<MeshType>::Clear(mesh);
      CCV[i].second->SetS();
      UpdateSelection<MeshType>::FaceConnectedFF(mesh);
      ret.push_back(std::make_shared<MeshType>());
      Append<MeshType, MeshType>::MeshCopy(*(ret.back()), mesh, true);
    }
    
    return ret;
  }
  
public:
  static const int VoroRelaxationStep = 20;
  
  ///
  /// \brief Remesh the main function that remeshes a mesh preserving creases.
  /// \param original the mesh
  /// \param samplingRadius is the sampling ragius for remeshing
  /// \param borderCreaseAngleDeg is the angle treshold for preserving corner points on the mesh boundary
  /// \return the remeshed mesh
  ///
  static inline MeshPtr Remesh(Mesh & original, const ScalarType samplingRadius, const ScalarType borderCreaseAngleDeg = 0.0)
  {
    UpdateTopology<Mesh>::FaceFace(original);
    UpdateFlags<Mesh>::FaceBorderFromFF(original);
    UpdateFlags<Mesh>::VertexBorderFromFaceAdj(original);

    RequireFFAdjacency(original);
	RequireVFAdjacency(original);
    
    if (Clean<Mesh>::CountNonManifoldEdgeFF(original) > 0)
    {
      std::cout << "Input mesh has non manifold edges" << std::endl;
      return nullptr;
    }
    
    // for closed watertight mesh try to split
    if (Clean<Mesh>::CountHoles(original) < 1)
    {
      CreaseCut<Mesh>(original, vcg::math::ToRad(borderCreaseAngleDeg));
      Allocator<Mesh>::CompactEveryVector(original);
      UpdateTopology<Mesh>::FaceFace(original);
      UpdateFlags<Mesh>::FaceBorderFromFF(original);
      UpdateFlags<Mesh>::VertexBorderFromFaceAdj(original);
#ifdef DEBUG_VORO
      io::Exporter<Mesh>::Save(original, "creaseSplit.ply", 0);
#endif
    }
    
    // One CC
    std::vector<MeshPtr> ccs = splitCC(original);
    if (ccs.empty())
      return RemeshOneCC(original, samplingRadius, borderCreaseAngleDeg);
    
    
    // Multiple CCs
    std::cout << "Remeshing " << ccs.size() << " components" << std::endl;
    for (size_t i=0; i<ccs.size(); i++)
    {
      std::cout << "Remeshing component " << (i+1) << "/" << ccs.size() << std::endl;
      ccs[i] = RemeshOneCC(*ccs[i], samplingRadius, borderCreaseAngleDeg, i);
    }
    
    MeshPtr ret = std::make_shared<Mesh>();
    for (MeshPtr & mesh : ccs)
    {
      Append<Mesh,Mesh>::Mesh(*ret, *mesh);
    }
    Clean<Mesh>::RemoveDuplicateVertex(*ret, true);
    return ret;
  }
  
  ///
  /// \brief RemeshOneCC the function that remeshes a single connected component mesh preserving its boundary (consistently for eventually adjacent meshes).
  /// \param original the mesh
  /// \param samplingRadius is the sampling ragius for remeshing
  /// \param borderCreaseAngleDeg is the angle treshold for preserving corner points on the mesh boundary
  /// \return the remeshed mesh
  ///
  static inline MeshPtr RemeshOneCC(Mesh & original, const ScalarType samplingRadius, const ScalarType borderCreaseAngleDeg = 0.0, int idx = 0)
  {
    RequireCompactness(original);
    RequirePerFaceFlags(original);
    
    UpdateTopology<Mesh>::FaceFace(original);
    UpdateFlags<Mesh>::FaceBorderFromFF(original);
    UpdateFlags<Mesh>::VertexBorderFromFaceAdj(original);
    
    
    // Resample border
    Mesh poissonEdgeMesh;
    {
      typedef typename EdgeMeshType::CoordType Coord;
      
      EdgeMeshType em;
      //			ThisType::ExtractMeshSides(original, em);
      
      ThisType::ExtractMeshBorders(original, em);
      
      // wtf we should close the loops
      Clean<EdgeMeshType>::RemoveDuplicateVertex(em);
      Allocator<EdgeMeshType>::CompactVertexVector(em);
      Allocator<EdgeMeshType>::CompactEdgeVector(em);
      
#ifdef DEBUG_VORO
      io::ExporterOBJ<EdgeMeshType>::Save(em, QString("edgeMesh_%1.obj").arg(idx).toStdString().c_str(), io::Mask::IOM_EDGEINDEX);
#endif
      
      // eventually split on 'creases'
      if (borderCreaseAngleDeg > 0.0)
      {
        UpdateFlags<EdgeMeshType>::VertexClearS(em);
        UpdateFlags<EdgeMeshType>::VertexClearV(em);
        Clean<EdgeMeshType>::SelectCreaseVertexOnEdgeMesh(em, vcg::math::ToRad(borderCreaseAngleDeg));
        std::cout << Clean<EdgeMeshType>::SplitSelectedVertexOnEdgeMesh(em) << " splits" << std::endl;
      }
#ifdef DEBUG_VORO
      io::ExporterOBJ<EdgeMeshType>::Save(em, QString("edgeMesh_split_%1.obj").arg(idx).toStdString().c_str(), io::Mask::IOM_EDGEINDEX);
#endif
      
      // Samples vector
      std::vector<Coord> borderSamples;
      TrivialSampler<EdgeMeshType> ps(borderSamples);
      
      // uniform sampling
      UpdateTopology<EdgeMeshType>::EdgeEdge(em);
      SurfaceSampling<EdgeMeshType>::EdgeMeshUniform(em, ps, samplingRadius, true);
      BuildMeshFromCoordVector(poissonEdgeMesh, borderSamples);
      UpdateBounding<Mesh>::Box(poissonEdgeMesh);
      
      // remove duplicate vertices
      Clean<Mesh>::RemoveDuplicateVertex(poissonEdgeMesh, false);
      Allocator<Mesh>::CompactVertexVector(poissonEdgeMesh);
      
      // select all vertices (to mark them fixed)
      UpdateFlags<Mesh>::VertexSetS(poissonEdgeMesh);
      
#ifdef DEBUG_VORO
      //			// temp remove
      //			UpdateColor<Mesh>::PerVertexConstant(poissonEdgeMesh, vcg::Color4b::Gray);
      
      //			typedef typename vcg::SpatialHashTable<VertexType, ScalarType> HashVertexGrid;
      //			HashVertexGrid HG;
      //			HG.Set(poissonEdgeMesh.vert.begin(),poissonEdgeMesh.vert.end());
      //			for (size_t i=0; i<creases.size(); i++)
      //			{
      //				const float dist_upper_bound=FLT_MAX;
      //				ScalarType dist;
      //				VertexType * vp = GetClosestVertex<MeshType,HashVertexGrid>(poissonEdgeMesh, HG, creases[i], dist_upper_bound, dist);
      //				assert(vp);
      //				vp->C() = vcg::Color4b::Red;
      //			}
      io::ExporterPLY<MeshType>::Save(poissonEdgeMesh, QString("borderMesh_%1.ply").arg(idx).toStdString().c_str(), io::Mask::IOM_VERTCOLOR);
#endif
    }
    
    typedef VoronoiProcessing<Mesh>            Voronoi;
    typedef TrivialSampler<Mesh>               BaseSampler;
    typedef SurfaceSampling<Mesh, BaseSampler> SurfaceSampler;
    typedef SurfaceSampling<Mesh, FixSampler>  SurfaceFixSampler;
    
    // copy original mesh
    Mesh baseMesh;
    Append<Mesh, Mesh>::MeshCopy(baseMesh, original, false, true);
    
    // refine to obtain a base mesh
    VoronoiProcessingParameter vpp;
    vpp.refinementRatio = 4.0f;
    Voronoi::PreprocessForVoronoi(baseMesh, samplingRadius, vpp);
    
    // Poisson sampling preserving border
    Mesh poissonMesh;
    std::vector<CoordType> seedPointVec;
    std::vector<bool>      seedFixedVec;
    FixSampler fix_sampler(seedPointVec, seedFixedVec);
    
    // montecarlo sampler
    std::vector<CoordType> sampleVec;
    BaseSampler mps(sampleVec);
    
    // NOTE in order to make results always the same the random sampling generator is reinitialized for
    // for each patch resampling
    SurfaceSampler::SamplingRandomGenerator().initialize(5489u);
    
    // Montecarlo oversampling
    Mesh montecarloMesh;
    int poissonCount = SurfaceSampler::ComputePoissonSampleNum(original, samplingRadius) * 0.7;
    
    std::cout << "poisson Count: " << poissonCount << std::endl;
    if (poissonCount <= 0)
    {
      // no need for internal sampling
      Append<Mesh, Mesh>::MeshCopy(poissonMesh, poissonEdgeMesh);
      for (auto vi = poissonEdgeMesh.vert.begin(); vi != poissonEdgeMesh.vert.end(); vi++)
      {
        fix_sampler.AddVert(*vi);
      }
    }
    else
    {
      // Montecarlo poisson sampling
      SurfaceSampler::MontecarloPoisson(original, mps, poissonCount * 20);
      BuildMeshFromCoordVector(montecarloMesh,sampleVec);
      
#ifdef DEBUG_VORO
      io::ExporterPLY<MeshType>::Save(montecarloMesh, QString("montecarloMesh_%1.ply").arg(idx).toStdString().c_str());
#endif
      
      // Poisson disk pruning initialized with edges
      typename SurfaceFixSampler::PoissonDiskParam pp;
      pp.preGenMesh = &poissonEdgeMesh;
      pp.preGenFlag = true;
      SurfaceFixSampler::PoissonDiskPruning(fix_sampler, montecarloMesh, samplingRadius, pp);
      
#ifdef DEBUG_VORO
      BuildMeshFromCoordVector(poissonMesh,seedPointVec);
      io::ExporterPLY<MeshType>::Save(poissonMesh, QString("poissonMesh_%1.ply").arg(idx).toStdString().c_str());
#endif
    }
    
    
    std::cout << "poisson samples " << seedPointVec.size() << std::endl;
    
    // restricted relaxation with fixed points
    vpp.seedPerturbationProbability = 0.0f;
    // TODO check preserveFixedSeed flag (NO)
    Voronoi::RestrictedVoronoiRelaxing(baseMesh, seedPointVec, seedFixedVec, VoroRelaxationStep, vpp);
    
#ifdef DEBUG_VORO
    BuildMeshFromCoordVector(poissonMesh,seedPointVec);
    io::ExporterPLY<MeshType>::Save(poissonMesh, QString("relaxedMesh_%1.ply").arg(idx).toStdString().c_str());
#endif
    
    // FAIL?
    MeshPtr finalMeshPtr = std::make_shared<Mesh>();
    std::vector<VertexType *> seedVertexVec;
    //		Voronoi::SeedToVertexConversion(baseMesh, seedPointVec, seedVertexVec, false);
    ThisType::SeedToFixedBorderVertexConversion(baseMesh, seedPointVec, seedFixedVec, seedVertexVec);
    EuclideanDistance<Mesh> dd;
    std::cout << "BEGIN compute vertex sources (basemesh vn:" << baseMesh.VN() << " fn:" << baseMesh.FN() << ")" << std::endl;
    
    Voronoi::ComputePerVertexSources(baseMesh, seedVertexVec, dd);
    std::cout << "END   compute vertex sources" << std::endl;
    //		Voronoi::ConvertDelaunayTriangulationToMesh(baseMesh, *finalMeshPtr, seedVertexVec, false); // traditional
    ThisType::ConvertDelaunayTriangulationExtendedToMesh(baseMesh, *finalMeshPtr, seedVertexVec); // border-preserving
    
#ifdef DEBUG_VORO
    io::ExporterPLY<MeshType>::Save(*finalMeshPtr, QString("voroMesh_%1.ply").arg(idx).toStdString().c_str());
    io::ExporterPLY<MeshType>::Save(baseMesh, QString("baseMesh_%1.ply").arg(idx).toStdString().c_str(), io::Mask::IOM_VERTCOLOR);
#endif
    
    return finalMeshPtr;
  }
  
protected:
  static inline void ExtractMeshBorders(Mesh & mesh, EdgeMeshType & sides)
  {
    RequireFFAdjacency(mesh);
    
    // clean the edge mesh containing the borders
    sides.Clear();
    
    // gather into separate vertices lists
    std::vector<std::vector<VertexType *> > edges;
    
    for (auto fi = mesh.face.begin(); fi != mesh.face.end(); fi++)
    {
      for (int e=0; e<fi->VN(); e++)
      {
        if (vcg::face::IsBorder(*fi, e))
        {
          std::vector<VertexType *> tmp;
          tmp.push_back(fi->V(e));
          tmp.push_back(fi->V((e+1)%fi->VN()));
          edges.push_back(tmp);
        }
      }
    }
    
    // convert to edge mesh
    for (auto & e : edges)
    {
      assert(e.size() >= 2);
      
      std::vector<typename EdgeMeshType::VertexType *> newVtx;
      
      // insert new vertices and store their pointer
      auto vi = Allocator<EdgeMeshType>::AddVertices(sides, e.size());
      for (const auto & v : e)
      {
        vi->ImportData(*v);
        newVtx.push_back(&(*vi++));
      }
      
      auto ei = Allocator<EdgeMeshType>::AddEdges(sides, e.size() - 1);
      for (int i=0; i<static_cast<int>(e.size() - 1); i++)
      {
        ei->V(0) = newVtx[i];
        ei->V(1) = newVtx[i+1];
        ei++;
      }
    }
    
    Clean<EdgeMeshType>::RemoveDuplicateVertex(sides);
  }
  
  
  ///
  /// \brief ExtractMeshSides
  /// \param mesh the mesh (topology already computed)
  /// \param sides the edge mesh filled with the extracted borders
  ///
  static inline void ExtractMeshSides(Mesh & mesh, EdgeMeshType & sides)
  {
    // TODO change this.... maybe wrong
    RequireFFAdjacency(mesh);
    RequireVFAdjacency(mesh);
    
    // clean the edge mesh containing the borders
    sides.Clear();
    
    // find a border edge
    assert(Clean<Mesh>::CountHoles(mesh) >= 1);
    PosType pos;
    for (auto fi = mesh.face.begin(); fi != mesh.face.end() && pos.IsNull(); fi++)
    {
      for (int e=0; e<fi->VN(); e++)
      {
        if (vcg::face::IsBorder(*fi, e))
        {
          pos = PosType(&(*fi), e);
          break;
        }
      }
    }
    assert(!pos.IsNull());
    assert(pos.IsBorder());
    
    // navigate to a corner
    VertexType * v = pos.V();
    do
    {
      pos.NextB();
    } while(!pos.V()->IsV() && pos.V() != v);
    // if it's a loop mark the initial point as a corner
    pos.V()->SetV();
    v = pos.V();
    
    // gather into separate vertices lists
    std::vector<std::vector<VertexType *> > edges;
    std::vector<VertexType *> edgePtrVec;
    do
    {
      edgePtrVec.push_back(pos.V());
      pos.NextB();
      if (pos.V()->IsV())
      {
        edgePtrVec.push_back(pos.V());
        edges.push_back(edgePtrVec);
        edgePtrVec.clear();
      }
    } while (pos.V() != v);
    
    // convert to edge mesh
    for (auto & e : edges)
    {
      assert(e.size() >= 2);
      
      std::vector<typename EdgeMeshType::VertexType *> newVtx;
      
      // insert new vertices and store their pointer
      auto vi = Allocator<EdgeMeshType>::AddVertices(sides, e.size());
      for (const auto & v : e)
      {
        vi->ImportData(*v);
        newVtx.push_back(&(*vi++));
      }
      
      auto ei = Allocator<EdgeMeshType>::AddEdges(sides, e.size() - 1);
      for (int i=0; i<static_cast<int>(e.size() - 1); i++)
      {
        ei->V(0) = newVtx[i];
        ei->V(1) = newVtx[i+1];
        ei++;
      }
    }
  }
  
  static void SeedToFixedBorderVertexConversion(MeshType & m,
                                                const std::vector<CoordType> & seedPVec,
                                                const std::vector<bool> & seedFixed,
                                                std::vector<VertexType *> & seedVVec)
  {
    // TODO mark here all seeds (cross-border)

    typedef typename vcg::SpatialHashTable<VertexType, ScalarType> HashVertexGrid;
    seedVVec.clear();
    
    UpdateTopology<MeshType>::FaceFace(m);
    UpdateFlags<MeshType>::VertexBorderFromFaceAdj(m);
    
    typename MeshType::BoxType bbox = m.bbox;
    bbox.Offset(bbox.Diag()/4.0);
    
    // internal vertices grid
    HashVertexGrid HG;
    HG.Set(m.vert.begin(),m.vert.end(), bbox);
    
    // boundary vertices grid
    MeshType borderMesh;
    HashVertexGrid borderHG;
    {
      // get border vertices and build another mesh
      std::vector<CoordType> borderPts;
      for (auto vit=m.vert.begin(); vit!=m.vert.end(); vit++)
      {
        if (!vit->IsD() && vit->IsB())
          borderPts.push_back(vit->cP());
      }
      BuildMeshFromCoordVector(borderMesh,borderPts);
      borderMesh.bbox = m.bbox;
      borderHG.Set(borderMesh.vert.begin(), borderMesh.vert.end(), bbox);
    }
    
    const float dist_upper_bound=m.bbox.Diag()/4.0;
    VertexType * vp = NULL;
    
    for( size_t i = 0; i < seedPVec.size(); i++)
    {
      const CoordType & p = seedPVec[i];
      const bool fixed    = seedFixed[i];
      
      if (!fixed)
      {
        ScalarType dist;
        vp = GetClosestVertex<MeshType,HashVertexGrid>(m, HG, p, dist_upper_bound, dist);
      }
      else
      {
        vp = NULL;
        
        ScalarType dist;
        VertexType * borderVp = GetClosestVertex<MeshType,HashVertexGrid>(borderMesh, borderHG, p, dist_upper_bound, dist);
        
        if (borderVp)
        {
          vp = GetClosestVertex<MeshType,HashVertexGrid>(m, HG, borderVp->cP(), dist_upper_bound, dist);
        }
      }
      
      if (vp)
      {
        seedVVec.push_back(vp);
      }
    }
  }
  
  static void ConvertDelaunayTriangulationExtendedToMesh(MeshType &m,
                                                         MeshType &outMesh,
                                                         std::vector<VertexType *> &seedVec)
  {
    typedef VoronoiProcessing<MeshType> Voronoi;
    
    RequirePerVertexAttribute(m ,"sources");
    RequireCompactness(m);
    RequireVFAdjacency(m);
    
    auto sources = Allocator<MeshType>::template GetPerVertexAttribute<VertexPointer> (m,"sources");
    
    outMesh.Clear();
    UpdateTopology<MeshType>::FaceFace(m);
    UpdateFlags<MeshType>::FaceBorderFromFF(m);
    
    std::map<VertexPointer, int> seedMap;  // It says if a given vertex of m is a seed (and its index in seedVec)
    Voronoi::BuildSeedMap(m, seedVec, seedMap);
    
    std::vector<FacePointer> innerCornerVec,   // Faces adjacent to three different regions
        borderCornerVec;  // Faces that are on the border and adjacent to at least two regions.
    Voronoi::GetFaceCornerVec(m, sources, innerCornerVec, borderCornerVec);
    
    // First add all the needed vertices: seeds and corners
    
    for(size_t i=0;i<seedVec.size();++i)
    {
      Allocator<MeshType>::AddVertex(outMesh, seedVec[i]->P(), vcg::Color4b::White);
    }
    
    // Now just add one face for each inner corner
    for(size_t i=0; i<innerCornerVec.size(); ++i)
    {
      VertexPointer s0 = sources[innerCornerVec[i]->V(0)];
      VertexPointer s1 = sources[innerCornerVec[i]->V(1)];
      VertexPointer s2 = sources[innerCornerVec[i]->V(2)];
      assert ( (s0!=s1) && (s0!=s2) && (s1!=s2) );
      VertexPointer v0 = & outMesh.vert[seedMap[s0]];
      VertexPointer v1 = & outMesh.vert[seedMap[s1]];
      VertexPointer v2 = & outMesh.vert[seedMap[s2]];
      Allocator<MeshType>::AddFace(outMesh,  v0, v1, v2 );
    }
    
    // Now loop around the borders and find the missing delaunay triangles
    // select border seed vertices only and pick one
    UpdateFlags<Mesh>::VertexBorderFromFaceAdj(m);
    UpdateFlags<Mesh>::VertexClearS(m);
    UpdateFlags<Mesh>::VertexClearV(m);
    
    std::vector<VertexPointer> borderSeeds;
    for (auto & s : seedVec)
    {
      if (s->IsB())
      {
        s->SetS();
        borderSeeds.emplace_back(s);
      }
    }
    
    for (VertexPointer startBorderVertex : borderSeeds)
    {
      if (startBorderVertex->IsV())
      {
        continue;
      }
      
      // unvisited border seed found
      
      // put the pos on the border
      PosType pos(startBorderVertex->VFp(), startBorderVertex->VFi());
      do {
        pos.NextE();
      } while (!pos.IsBorder() || (pos.VInd() != pos.E()));
      
      // check all border edges between each consecutive border seeds pair
      do {
        std::vector<VertexPointer> edgeVoroVertices(1, sources[pos.V()]);
        //	among all sources found
        do {
          pos.NextB();
          VertexPointer source = sources[pos.V()];
          if (edgeVoroVertices.empty() || edgeVoroVertices.back() != source)
          {
            edgeVoroVertices.push_back(source);
          }
        } while (!pos.V()->IsS());
        
        pos.V()->SetV();
        
        //				assert(edgeVoroVertices.size() >= 2);
        
        // add face if 3 different voronoi regions are crossed by the edge
        if (edgeVoroVertices.size() == 3)
        {
          VertexPointer v0 = & outMesh.vert[seedMap[edgeVoroVertices[0]]];
          VertexPointer v1 = & outMesh.vert[seedMap[edgeVoroVertices[1]]];
          VertexPointer v2 = & outMesh.vert[seedMap[edgeVoroVertices[2]]];
          Allocator<MeshType>::AddFace(outMesh, v0,v1,v2);
        }
        
      } while ((pos.V() != startBorderVertex));
    }
    
    
    Clean<MeshType>::RemoveUnreferencedVertex(outMesh);
    Allocator<MeshType>::CompactVertexVector(outMesh);
  }
  
  ///
  /// \brief The FixSampler class is used with poisson disk pruning to preserve selected vertices and
  /// keep an auxiliary vector indicating wether the sample is fixed or not
  ///
  class FixSampler
  {
  public:
    typedef typename MeshType::CoordType  CoordType;
    typedef typename MeshType::VertexType VertexType;
    
    FixSampler(std::vector<CoordType> & samples,
               std::vector<bool>      & fixed)
      : sampleVec(samples)
      , fixedVec (fixed)
    {
      reset();
    }
    
    void reset()
    {
      sampleVec.clear();
      fixedVec .clear();
    }
    
    void AddVert(const VertexType &p)
    {
      sampleVec.push_back(p.cP());
      fixedVec .push_back(p.IsS());
    }
    
  private:
    std::vector<CoordType> & sampleVec;
    std::vector<bool>      & fixedVec;
  };

};
} // end namespace tri
} // end namespace vcg

#endif // _VCGLIB_VORONOI_REMESHER_H
