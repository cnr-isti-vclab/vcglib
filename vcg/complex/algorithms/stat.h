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

#ifndef __VCGLIB_TRIMESH_STAT
#define __VCGLIB_TRIMESH_STAT

// Standard headers
// VCG headers

#include <vcg/math/histogram.h>
#include <vcg/complex/algorithms/closest.h>
#include <vcg/space/index/grid_static_ptr.h>
#include <vcg/complex/algorithms/inertia.h>
#include <vcg/space/polygon3.h>


namespace vcg {
namespace tri{
template <class StatMeshType>
class Stat
{
public:
  typedef StatMeshType MeshType;
  typedef typename MeshType::ScalarType			ScalarType;
  typedef typename MeshType::VertexType     VertexType;
  typedef typename MeshType::VertexPointer  VertexPointer;
  typedef typename MeshType::VertexIterator VertexIterator;
  typedef typename MeshType::ConstVertexIterator ConstVertexIterator;
  typedef typename MeshType::EdgeType       EdgeType;
  typedef typename MeshType::EdgeIterator   EdgeIterator;
  typedef typename MeshType::FaceType       FaceType;
  typedef typename MeshType::FacePointer    FacePointer;
  typedef typename MeshType::FaceIterator   FaceIterator;
  typedef typename MeshType::FaceContainer  FaceContainer;
  typedef typename MeshType::TetraType      TetraType;
  typedef typename MeshType::TetraPointer   TetraPointer;
  typedef typename MeshType::TetraIterator  TetraIterator;
  typedef typename MeshType::TetraContainer TetraContainer;
  typedef typename vcg::Box3<ScalarType>  Box3Type;

  static void ComputePerVertexQualityMinMax(MeshType & m, ScalarType &minV, ScalarType &maxV)
  {
    std::pair<ScalarType, ScalarType> pp = ComputePerVertexQualityMinMax(m);
    
    minV=pp.first;
    maxV=pp.second;
  }
  static std::pair<ScalarType, ScalarType> ComputePerVertexQualityMinMax(MeshType & m)
  {
//    assert(0);
    tri::RequirePerVertexQuality(m);
    typename MeshType::template PerMeshAttributeHandle  < std::pair<ScalarType, ScalarType> > mmqH;
    mmqH = tri::Allocator<MeshType>::template GetPerMeshAttribute <std::pair<ScalarType, ScalarType> >(m,"minmaxQ");

    std::pair<ScalarType, ScalarType> minmax = std::make_pair(std::numeric_limits<ScalarType>::max(), -std::numeric_limits<ScalarType>::max());

    for(VertexIterator vi = m.vert.begin(); vi != m.vert.end(); ++vi)
      if(!(*vi).IsD())
      {
        if( (*vi).Q() < minmax.first)  minmax.first  = (*vi).Q();
        if( (*vi).Q() > minmax.second) minmax.second = (*vi).Q();
      }

    mmqH() = minmax;
    return minmax;
  }


  static void ComputePerFaceQualityMinMax(MeshType & m, ScalarType &minV, ScalarType &maxV)
  {
    std::pair<ScalarType, ScalarType> pp = ComputePerFaceQualityMinMax(m);
    minV=pp.first; 
    maxV=pp.second;
  }

  static std::pair<ScalarType,ScalarType> ComputePerFaceQualityMinMax( MeshType & m)
  {
    tri::RequirePerFaceQuality(m);
    std::pair<ScalarType,ScalarType> minmax = std::make_pair(std::numeric_limits<ScalarType>::max(),-std::numeric_limits<ScalarType>::max());

    FaceIterator fi;
    for(fi = m.face.begin(); fi != m.face.end(); ++fi)
      if(!(*fi).IsD())
      {
        if( (*fi).Q() < minmax.first)  minmax.first  = (*fi).Q();
        if( (*fi).Q() > minmax.second) minmax.second = (*fi).Q();
      }
    return minmax;
  }

  static void ComputePerTetraQualityMinMax(MeshType & m, ScalarType & minQ, ScalarType & maxQ)
  {
    std::pair<ScalarType, ScalarType> minmax = ComputerPerTetraQualityMinMax(m);

    minQ = minmax.first;
    maxQ = minmax.second;
  }

  static std::pair<ScalarType, ScalarType> ComputePerTetraQualityMinMax(MeshType & m)
  {
    tri::RequirePerTetraQuality(m);
	std::pair<ScalarType, ScalarType> minmax = std::make_pair(std::numeric_limits<ScalarType>::max(), std::numeric_limits<ScalarType>::min());

    ForEachTetra(m, [&minmax] (TetraType & t) {
      if (t.Q() < minmax.first)  minmax.first  = t.Q();
      if (t.Q() > minmax.second) minmax.second = t.Q();
    });

    return minmax;
  }

  static ScalarType ComputePerTetraQualityAvg(MeshType & m)
  {
    tri::RequirePerTetraQuality(m);
    ScalarType avgQ = 0;

    ForEachTetra(m, [&avgQ] (TetraType & t) {
      avgQ += t.Q();
    });

    return avgQ /= (ScalarType) m.TN();
  }

  static ScalarType ComputePerFaceQualityAvg(MeshType & m)
  {
    tri::RequirePerFaceQuality(m);
    ScalarType AvgQ = 0;

    FaceIterator fi;
    size_t num=0;
    for(fi = m.face.begin(); fi != m.face.end(); ++fi)
    {
      if((*fi).IsD())continue;
        AvgQ+= (*fi).Q();
        num++;
    }
    return (AvgQ/(ScalarType)num);
  }

  static ScalarType ComputePerVertQualityAvg(const MeshType & m)
  {
    tri::RequirePerVertexQuality(m);
    ScalarType AvgQ = 0;

    ConstVertexIterator vi;
    size_t num=0;
    for(vi = m.vert.begin(); vi != m.vert.end(); ++vi)
    {
      if((*vi).IsD())continue;
        AvgQ+= (*vi).cQ();
        num++;
    }
    return (AvgQ/(ScalarType)num);
  }

  static std::pair<ScalarType,ScalarType> ComputePerEdgeQualityMinMax( MeshType & m)
  {
    tri::RequirePerEdgeQuality(m);
    std::pair<ScalarType,ScalarType> minmax = std::make_pair(std::numeric_limits<ScalarType>::max(),-std::numeric_limits<ScalarType>::max());

    EdgeIterator ei;
    for(ei = m.edge.begin(); ei != m.edge.end(); ++ei)
      if(!(*ei).IsD())
      {
        if( (*ei).Q() < minmax.first)  minmax.first =(*ei).Q();
        if( (*ei).Q() > minmax.second) minmax.second=(*ei).Q();
      }
    return minmax;
  }

  /**
  \short compute the pointcloud barycenter.
  E.g. it assume each vertex has a mass. If useQualityAsWeight is true, vertex quality is the mass of the vertices
  */
	static Point3<ScalarType> ComputeCloudBarycenter(const MeshType & m, bool useQualityAsWeight=false) 
	{
		if (useQualityAsWeight)
			tri::RequirePerVertexQuality(m);
	 
		Point3<ScalarType> barycenter(0, 0, 0);
		Point3d accumulator(0.0, 0.0, 0.0);
		double weightSum = 0;
		for (auto vi = m.vert.begin(); vi != m.vert.end(); ++vi) {
			if (!(*vi).IsD()) {
				ScalarType weight = useQualityAsWeight ? (*vi).Q() : 1.0f;
				accumulator[0] += (double)((*vi).P()[0] * weight);
				accumulator[1] += (double)((*vi).P()[1] * weight);
				accumulator[2] += (double)((*vi).P()[2] * weight);
				weightSum += weight;
			}
		}
		barycenter[0] = (ScalarType)(accumulator[0] / weightSum);
		barycenter[1] = (ScalarType)(accumulator[1] / weightSum);
		barycenter[2] = (ScalarType)(accumulator[2] / weightSum);
		return barycenter;
  }

  /**
    \short compute the barycenter of the surface thin-shell.
    E.g. it assume a 'empty' model where all the mass is located on the surface and compute the barycenter of that thinshell.
    Works for any triangulated model (no problem with open, nonmanifold selfintersecting models).
    Useful for computing the barycenter of 2D planar figures.
    */
  static Point3<ScalarType> ComputeShellBarycenter(MeshType & m)
  {
    Point3<ScalarType> barycenter(0,0,0);
    ScalarType areaSum=0;
    FaceIterator fi;
    for(fi = m.face.begin(); fi != m.face.end(); ++fi)
      if(!(*fi).IsD())
      {
        ScalarType area=DoubleArea(*fi);
        barycenter += Barycenter(*fi)*area;
        areaSum+=area;
      }
    return barycenter/areaSum;
  }

  static ScalarType ComputeTetraMeshVolume(MeshType & m)
  {
    ScalarType V = 0;

    ForEachTetra(m, [&V] (TetraType & t) {
      V += Tetra::ComputeVolume(t);
    });

    return V;
  }

  static ScalarType ComputeMeshVolume(const MeshType & m)
  {
    Inertia<MeshType> I(m);
    return I.Mass();
  }

  static ScalarType ComputeMeshArea(const MeshType & m)
  {
    ScalarType area=0;

    for(auto fi = m.face.begin(); fi != m.face.end(); ++fi)
      if(!(*fi).IsD())
        area += DoubleArea(*fi);

    return area/ScalarType(2.0);
  }

  static ScalarType ComputePolyMeshArea(MeshType & m)
  {
    ScalarType area=0;

    for(FaceIterator fi = m.face.begin(); fi != m.face.end(); ++fi)
      if(!(*fi).IsD())
        area += PolyArea(*fi);

    return area;
  }

	static ScalarType ComputeBorderLength(MeshType & m, bool computeFFTopology = true)
	{
		RequireFFAdjacency(m);
		ScalarType sum = 0;
		if (computeFFTopology)
		{
			tri::UpdateTopology<MeshType>::FaceFace(m);
		}
		ForEachFace(m, [&](FaceType &f) {
			for (int k=0; k<f.VN(); k++)
				if (face::IsBorder(f, k))
				{
					sum += Distance(f.cP0(k), f.cP1(k));
				}
		});
		return sum;
	}

  static void ComputePerVertexQualityDistribution(MeshType & m, Distribution<ScalarType> & h, bool selectionOnly = false)    // V1.0
  {
    tri::RequirePerVertexQuality(m);
    for(VertexIterator vi = m.vert.begin(); vi != m.vert.end(); ++vi)
      if(!(*vi).IsD() &&  ((!selectionOnly) || (*vi).IsS()) )
      {
        assert(!math::IsNAN((*vi).Q()) && "You should never try to compute Histogram with Invalid Floating points numbers (NaN)");
        h.Add((*vi).Q());
      }
  }

  static void ComputePerFaceQualityDistribution( MeshType & m,  Distribution<typename MeshType::ScalarType> &h,
                                                 bool selectionOnly = false)    // V1.0
  {
    tri::RequirePerFaceQuality(m);
    for(FaceIterator fi = m.face.begin(); fi != m.face.end(); ++fi)
      if(!(*fi).IsD() &&  ((!selectionOnly) || (*fi).IsS()) )
      {
        assert(!math::IsNAN((*fi).Q()) && "You should never try to compute Histogram with Invalid Floating points numbers (NaN)");
        h.Add((*fi).Q());
      }
  }

  static void ComputePerTetraQualityDistribution(MeshType & m, Distribution<ScalarType> & h, bool selectionOnly = false)
  {
    tri::RequirePerTetraQuality(m);
    ForEachTetra(m, [&] (TetraType & t) {
      if (!selectionOnly || t.IsS())
      {
        assert(!math::IsNAN(t.Q()) && "You should never try to compute Histogram with Invalid Floating points numbers (NaN)");
        h.Add(t.Q());
      }
    });
  }

  static void ComputePerTetraQualityHistogram(MeshType & m, Histogram<ScalarType> & h, bool selectionOnly = false, int HistSize = 10000)
  {
    tri::RequirePerTetraQuality(m);
	std::pair<ScalarType, ScalarType> minmax = tri::Stat<MeshType>::ComputePerTetraQualityMinMax(m);
    
    h.Clear();
    h.SetRange(minmax.first, minmax.second, HistSize);

    ForEachTetra(m, [&] (TetraType & t) {
      if (!selectionOnly || t.IsS())
      {
        assert(!math::IsNAN(t.Q()) && "You should never try to compute Histogram with Invalid Floating points numbers (NaN)");
        h.Add(t.Q());
      }
    });
  }

  static void ComputePerFaceQualityHistogram( MeshType & m, Histogram<ScalarType> &h, bool selectionOnly=false,int HistSize=10000 )
  {
    tri::RequirePerFaceQuality(m);
    std::pair<ScalarType, ScalarType> minmax = tri::Stat<MeshType>::ComputePerFaceQualityMinMax(m);
    h.Clear();
    h.SetRange( minmax.first,minmax.second, HistSize );
    for(FaceIterator fi = m.face.begin(); fi != m.face.end(); ++fi)
      if(!(*fi).IsD() &&  ((!selectionOnly) || (*fi).IsS()) ){
        assert(!math::IsNAN((*fi).Q()) && "You should never try to compute Histogram with Invalid Floating points numbers (NaN)");
        h.Add((*fi).Q());
      }
  }

  static void ComputePerVertexQualityHistogram( MeshType & m, Histogram<ScalarType> &h, bool selectionOnly = false, int HistSize=10000 )    // V1.0
  {
    tri::RequirePerVertexQuality(m);
    std::pair<ScalarType, ScalarType> minmax = ComputePerVertexQualityMinMax(m);

    h.Clear();
    h.SetRange( minmax.first,minmax.second, HistSize);
    for(VertexIterator vi = m.vert.begin(); vi != m.vert.end(); ++vi)
      if(!(*vi).IsD() &&  ((!selectionOnly) || (*vi).IsS()) )
      {
        assert(!math::IsNAN((*vi).Q()) && "You should never try to compute Histogram with Invalid Floating points numbers (NaN)");
        h.Add((*vi).Q());
      }
    // Sanity check; If some very wrong value has happened in the Q value,
    // the histogram is messed. If a significant percentage (20% )of the values are all in a single bin
    // we should try to solve the problem. No easy solution here.
    // We choose to compute the get the 1percentile and 99 percentile values as new mixmax ranges
    // and just to be sure enlarge the Histogram.

    if(h.MaxCount() > HistSize/5)
    {
      std::vector<ScalarType> QV;
      QV.reserve(m.vn);
      for(VertexIterator vi = m.vert.begin(); vi != m.vert.end(); ++vi)
        if(!(*vi).IsD()) QV.push_back((*vi).Q());

      std::nth_element(QV.begin(),QV.begin()+m.vn/100,QV.end());
      ScalarType newmin=*(QV.begin()+m.vn/100);
      std::nth_element(QV.begin(),QV.begin()+m.vn-m.vn/100,QV.end());
      ScalarType newmax=*(QV.begin()+m.vn-m.vn/100);

      h.Clear();
      h.SetRange(newmin, newmax, HistSize*50);
      for(VertexIterator vi = m.vert.begin(); vi != m.vert.end(); ++vi)
        if(!(*vi).IsD() && ((!selectionOnly) || (*vi).IsS()) )
          h.Add((*vi).Q());
    }
  }

  static void ComputeEdgeLengthHistogram(MeshType & m, Histogram<ScalarType> & h)
  {
    assert(m.edge.size()>0);
    h.Clear();
    h.SetRange( 0, m.bbox.Diag(), 10000);
    for(EdgeIterator ei = m.edge.begin(); ei != m.edge.end(); ++ei)
    {
      if(!(*ei).IsD())
      {
        h.Add(Distance<ScalarType>((*ei).V(0)->P(),(*ei).V(1)->P()));
      }
    }
  }

  static ScalarType ComputeEdgeLengthAverage(MeshType & m)
  {
    Histogram<ScalarType> h;
    ComputeEdgeLengthHistogram(m,h);
    return h.Avg();
  }

  static ScalarType ComputeEdgeLengthSum(MeshType & m)
  {
    ScalarType sum=0;
    ForEachEdge(m, [&](EdgeType &e){
      sum+=Distance(e.cP(0),e.cP(1));
    });    
    return sum;
  }
  static void ComputeFaceEdgeLengthDistribution( MeshType & m, Distribution<ScalarType> & h, bool includeFauxEdge=false)
  {
    std::vector< typename tri::UpdateTopology<MeshType>::PEdge > edgeVec;
    tri::UpdateTopology<MeshType>::FillUniqueEdgeVector(m,edgeVec,includeFauxEdge);
    h.Clear();
    tri::UpdateFlags<MeshType>::FaceBorderFromNone(m);
    for(size_t i=0;i<edgeVec.size();++i)
      h.Add(Distance(edgeVec[i].v[0]->P(),edgeVec[i].v[1]->P()));
  }

  static ScalarType ComputeFaceEdgeLengthAverage(MeshType & m, bool selected=false)
  {
    double sum=0;
    for(FaceIterator fi = m.face.begin(); fi != m.face.end(); ++fi)
      if(!(*fi).IsD())
        if(!selected || fi->IsS())
        {
          for(int i=0;i<3;++i)
            sum+=double(Distance(fi->P0(i),fi->P1(i)));
        }
    return sum/(m.fn*3.0);
  }

}; // end class

} //End Namespace tri
} // End Namespace vcg

#endif

