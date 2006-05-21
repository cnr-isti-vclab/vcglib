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
****************************************************************************/

#ifndef __VCGLIB_TRIMESH_STAT
#define __VCGLIB_TRIMESH_STAT

// Standard headers
// VCG headers

#include <vcg/math/histogram.h>
#include <vcg/simplex/face/face.h>
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
        std::pair<float,float> minmax = make_pair(std::numeric_limits<float>::max(),-std::numeric_limits<float>::max());
         
        VertexIterator vi;
        for(vi = m.vert.begin(); vi != m.vert.end(); ++vi)
          if(!(*vi).IsD()) 
          {
            if( (*vi).Q() < minmax.first) minmax.first=(*vi).Q();
            if( (*vi).Q() > minmax.second) minmax.second=(*vi).Q();
          }
        return minmax;
      }

      static void ComputePerVertexQualityHistogram( MeshType & m, Histogramf &h)    // V1.0
      {
        h.Clear();
        std::pair<float,float> minmax = ComputePerVertexQualityMinMax(m);
        h.SetRange( minmax.first,minmax.second, 10000);
        VertexIterator vi;
        for(vi = m.vert.begin(); vi != m.vert.end(); ++vi)
          if(!(*vi).IsD()) h.Add((*vi).Q());
      }

			static int ComputeEdgeHistogram( MeshType & m, Histogramf &h)    // V1.0
      {
        ScalarType Diag = m.bbox.Diag();
        h.clear();
        h.SetRange( 0, diagonale, 10000);
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
        for(vi = m.cm.vert.begin(); vi != m.cm.vert.end(); ++vi)
          (*vi).ClearS();
      }
}; // end class
	
	} //End Namespace tri
} // End Namespace vcg

#endif

