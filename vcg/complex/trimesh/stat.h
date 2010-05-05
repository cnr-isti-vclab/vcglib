/****************************************************************************
* VCGLib                                                            o o     *
* Visual and Computer Graphics Library                            o     o   *
*                                                                _   O  _   *
* Copyright(C) 2006                                                \/)\/    *
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
/****************************************************************************
History

$Log: not supported by cvs2svn $
Revision 1.3  2007/01/11 10:12:19  cignoni
Removed useless and conflicting inclusion of face.h

Revision 1.2  2006/05/25 09:39:09  cignoni
missing std and other gcc detected syntax errors

Revision 1.1  2006/05/21 06:59:13  cignoni
Initial Commit

****************************************************************************/

#ifndef __VCGLIB_TRIMESH_STAT
#define __VCGLIB_TRIMESH_STAT

// Standard headers
// VCG headers

#include <vcg/math/histogram.h>
#include <vcg/simplex/face/pos.h>
#include <vcg/simplex/face/topology.h>
#include <vcg/complex/trimesh/base.h>
#include <vcg/complex/trimesh/closest.h>
#include <vcg/space/index/grid_static_ptr.h>
#include <vcg/complex/trimesh/allocate.h>


namespace vcg {
	namespace tri{
template <class StatMeshType>
class Stat
{
  public:
	typedef StatMeshType MeshType; 
	typedef typename MeshType::VertexType     VertexType;
	typedef typename MeshType::VertexPointer  VertexPointer;
	typedef typename MeshType::VertexIterator VertexIterator;
	typedef	typename MeshType::ScalarType			ScalarType;
	typedef typename MeshType::FaceType       FaceType;
	typedef typename MeshType::FacePointer    FacePointer;
	typedef typename MeshType::FaceIterator   FaceIterator;
	typedef typename MeshType::FaceContainer  FaceContainer;
	typedef typename vcg::Box3<ScalarType>  Box3Type;
	
	static std::pair<float,float> ComputePerVertexQualityMinMax( MeshType & m)    // V1.0
	{
	std::pair<float,float> minmax = std::make_pair(std::numeric_limits<float>::max(),-std::numeric_limits<float>::max());
		 
		VertexIterator vi;
		for(vi = m.vert.begin(); vi != m.vert.end(); ++vi)
			if(!(*vi).IsD()) 
			{
				if( (*vi).Q() < minmax.first) minmax.first=(*vi).Q();
				if( (*vi).Q() > minmax.second) minmax.second=(*vi).Q();
			}
		return minmax;
	}

	static std::pair<float,float> ComputePerFaceQualityMinMax( MeshType & m)    // V1.0
	{
		std::pair<float,float> minmax = std::make_pair(std::numeric_limits<float>::max(),-std::numeric_limits<float>::max());
		
		FaceIterator fi;
		for(fi = m.face.begin(); fi != m.face.end(); ++fi)
			if(!(*fi).IsD()) 
			{
				if( (*fi).Q() < minmax.first)  minmax.first =(*fi).Q();
				if( (*fi).Q() > minmax.second) minmax.second=(*fi).Q();
			}
		return minmax;
	}
	
			static ScalarType ComputeMeshArea(MeshType & m)
			{
				ScalarType area=0;
				
				FaceIterator fi;
        for(fi = m.face.begin(); fi != m.face.end(); ++fi)
						if(!(*fi).IsD()) 
								area += DoubleArea(*fi);
						
				return area/ScalarType(2.0);
			}
			
      static void ComputePerVertexQualityDistribution( MeshType & m, Distribution<float> &h, bool selectionOnly = false)    // V1.0
      {
        VertexIterator vi;
        for(vi = m.vert.begin(); vi != m.vert.end(); ++vi)
            if(!(*vi).IsD() &&  ((!selectionOnly) || (*vi).IsS()) )
              {
                assert(!math::IsNAN((*vi).Q()) && "You should never try to compute Histogram with Invalid Floating points numbers (NaN)");
                h.Add((*vi).Q());
              }
      }

      static void ComputePerFaceQualityDistribution( MeshType & m, Distribution<float> &h, bool selectionOnly = false)    // V1.0
      {
        FaceIterator fi;
        for(fi = m.face.begin(); fi != m.face.end(); ++fi)
            if(!(*fi).IsD() &&  ((!selectionOnly) || (*fi).IsS()) )
              {
                assert(!math::IsNAN((*fi).Q()) && "You should never try to compute Histogram with Invalid Floating points numbers (NaN)");
                h.Add((*fi).Q());
              }
      }

      static void ComputePerFaceQualityHistogram( MeshType & m, Histogramf &h, bool selectionOnly=false )
      {
        FaceIterator fi;
        int HistSize=10000;
        std::pair<float,float> minmax = tri::Stat<MeshType>::ComputePerFaceQualityMinMax(m);
        h.Clear();
        h.SetRange( minmax.first,minmax.second, HistSize );
        for(fi = m.face.begin(); fi != m.face.end(); ++fi)
          if(!(*fi).IsD() &&  ((!selectionOnly) || (*fi).IsS()) ){
            assert(!math::IsNAN((*fi).Q()) && "You should never try to compute Histogram with Invalid Floating points numbers (NaN)");
          h.Add((*fi).Q());
        }
      }

      static void ComputePerVertexQualityHistogram( MeshType & m, Histogramf &h, bool selectionOnly = false)    // V1.0
      {
				VertexIterator vi;
				int HistSize=10000;
				std::pair<float,float> minmax = ComputePerVertexQualityMinMax(m);
					
				h.Clear();
				h.SetRange( minmax.first,minmax.second, HistSize);
				for(vi = m.vert.begin(); vi != m.vert.end(); ++vi)
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
						std::vector<float> QV;
						QV.reserve(m.vn);
						for(vi = m.vert.begin(); vi != m.vert.end(); ++vi)
							if(!(*vi).IsD()) QV.push_back((*vi).Q());

						std::nth_element(QV.begin(),QV.begin()+m.vn/100,QV.end());
						float newmin=*(QV.begin()+m.vn/100);
						std::nth_element(QV.begin(),QV.begin()+m.vn-m.vn/100,QV.end());
						float newmax=*(QV.begin()+m.vn-m.vn/100);
						
						h.Clear();
						h.SetRange(newmin, newmax, HistSize*50);
						for(vi = m.vert.begin(); vi != m.vert.end(); ++vi)
                                                                if(!(*vi).IsD() && ((!selectionOnly) || (*vi).IsS()) )
                                                                    h.Add((*vi).Q());
					}
      }

			static int ComputeEdgeHistogram( MeshType & m, Histogramf &h)    // V1.0
      {
        ScalarType Diag = m.bbox.Diag();
        h.Clear();
        h.SetRange( 0, Diag, 10000);
        FaceIterator fi;
        for(fi = m.face.begin(); fi != m.face.end(); ++fi)
        {
          if(!(*fi).IsD())
          {
            if( !(*fi).V(0)->IsS() && !(*fi).V(1)->IsS()  )
            {
              h.Add(Distance<float>((*fi).V(0)->P(),(*fi).V(1)->P()));
              (*fi).V(0)->SetS();
              (*fi).V(1)->SetS();
            }
            if( !(*fi).V(1)->IsS() && !(*fi).V(2)->IsS())
            {
              h.Add(Distance<float>((*fi).V(1)->P(),(*fi).V(2)->P()));
              (*fi).V(2)->SetS();
              (*fi).V(1)->SetS();
            }
            if( !(*fi).V(2)->IsS() && !(*fi).V(0)->IsS())
            {
              h.Add(Distance<float>((*fi).V(2)->P(),(*fi).V(0)->P()));
              (*fi).V(0)->SetS();
              (*fi).V(2)->SetS();
            }
          }
        }
        VertexIterator vi;
        for(vi = m.vert.begin(); vi != m.vert.end(); ++vi)
          (*vi).ClearS();
      return 0;
	  }
}; // end class
	
	} //End Namespace tri
} // End Namespace vcg

#endif

