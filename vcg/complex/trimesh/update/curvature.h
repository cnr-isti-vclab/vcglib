/****************************************************************************
* VCGLib                                                            o o     *
* Visual and Computer Graphics Library                            o     o   *
*                                                                _   O  _   *
* Copyright(C) 2004                                                \/)\/    *
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
/****************************************************************************/

#ifndef VCGLIB_UPDATE_CURVATURE_
#define VCGLIB_UPDATE_CURVATURE_

#include <vcg/math/base.h>
#include <vcg/simplex/face/pos.h>

namespace vcg {
namespace tri {

/** \addtogroup trimesh */
/*@{*/

/// Management, updating and computation of per-vertex and per-face normals.
/// This class is used to compute or update the normals that can be stored in the vertex or face component of a mesh.
template <class ComputeMeshType>
class UpdateCurvature
{

public:
typedef ComputeMeshType MeshType; 	
typedef typename MeshType::VertexType     VertexType;
typedef typename VertexType::NormalType     NormalType;
typedef typename VertexType::ScalarType ScalarType;
typedef typename MeshType::VertexPointer  VertexPointer;
typedef typename MeshType::VertexIterator VertexIterator;
typedef typename MeshType::FaceType       FaceType;
typedef typename MeshType::FacePointer    FacePointer;
typedef typename MeshType::FaceIterator   FaceIterator;


/** computes the discrete gaussian curvature as proposed in 
Discrete Differential-Geometry Operators for Triangulated 2-Manifolds Mark Meyer,
 Mathieu Desbrun, Peter Schroder, Alan H. Barr VisMath '02, Berlin
*/
static void Gaussian( MeshType &  m){
		assert(m.HasPerVertexQuality());

	MeshType::VertexIterator vi;   // iteratore vertice
	MeshType::FaceIterator fi;     // iteratore facce
	double *area;                  // areamix vector
	int i;												 // index
	double area0, area1, area2;
	double angle0, angle1, angle2; 
	
	//--- Initialization
	area = new double[m.vn];

	//reset the values to 0
	for(vi=m.vert.begin();vi!=m.vert.end();++vi) if(!(*vi).IsD())
		(*vi).Q() = 0.0;

	//--- compute Areamix
	for(fi=m.face.begin();fi!=m.face.end();++fi) if(!(*fi).IsD())
	{
		
		// angles
			 angle0 = math::Abs(Angle(	(*fi).V(1)->P()-(*fi).V(0)->P(),(*fi).V(2)->P()-(*fi).V(0)->P() ));
			 angle1 = math::Abs(Angle(	(*fi).V(0)->P()-(*fi).V(1)->P(),(*fi).V(2)->P()-(*fi).V(1)->P() ));
 			 angle2 = M_PI-(angle0+angle1);
		
		if((angle0 < M_PI/2) || (angle1 < M_PI/2) || (angle2 < M_PI/2))  // triangolo non ottuso
		{ 
			float e01 = SquaredDistance( (*fi).V(1)->P() , (*fi).V(0)->P() );
			float e12 = SquaredDistance( (*fi).V(2)->P() , (*fi).V(1)->P() );
			float e20 = SquaredDistance( (*fi).V(0)->P() , (*fi).V(2)->P() );
			
			// voronoi area v[0]
			area0 = ( e01*(1/tan(angle2)) + e20*(1/tan(angle1)) ) /8;
			// voronoi area v[1]
			area1 = ( e01*(1/tan(angle2)) + e12*(1/tan(angle0)) ) /8;
			// voronoi area v[2]
			area2 = ( e20*(1/tan(angle1)) + e20*(1/tan(angle0)) ) /8;
			
			(*fi).V(0)->Q()  += area0;
			(*fi).V(1)->Q()  += area1;
			(*fi).V(2)->Q()  += area2;
		}
		else // triangolo ottuso
		{  
			(*fi).V(0)->Q() += (*fi).Area() / 3;
			(*fi).V(1)->Q() += (*fi).Area() / 3;
			(*fi).V(2)->Q() += (*fi).Area() / 3;            
		}
	}

	i = 0;
	for(vi=m.vert.begin();vi!=m.vert.end();++vi,++i) if(!(*vi).IsD())
	{
		area[i] = (*vi).Q();
		(*vi).Q() = (float)(2.0 * M_PI);
	}
	
	if(false)
	for(fi=m.face.begin();fi!=m.face.end();++fi) if(!(*fi).IsD())
	{
		for(int i=0;i<3;i++)
		{
			if((*fi).IsBorder(i))
			{
				MeshType::CoordType e1,e2;
				vcg::face::Pos<FaceType> hp(&*fi,i,(*fi).V(i));
				//MeshType::hedgepos_type hp(&*fi,i,(*fi).V(i));
				vcg::face::Pos<FaceType> hp1=hp;
				//MeshType::hedgepos_type hp1=hp;

				hp1.FlipV();
	
				e1= hp1.v->P()-hp.v->P();
				hp1.FlipV();
				hp1.NextB();
				e2= hp1.v->P()-hp.v->P();
				(*fi).V(i)->Q() -=math::Abs(Angle(e1,e2));
			}
	  }
	}
	
	
	for(fi=m.face.begin();fi!=m.face.end();++fi)  if(!(*fi).IsD())
	{
		float angle0 = math::Abs(Angle(
			(*fi).V(1)->P()-(*fi).V(0)->P(),(*fi).V(2)->P()-(*fi).V(0)->P() ));
		float angle1 = math::Abs(Angle(
			(*fi).V(0)->P()-(*fi).V(1)->P(),(*fi).V(2)->P()-(*fi).V(1)->P() ));
		float angle2 = M_PI-(angle0+angle1);
		
		(*fi).V(0)->Q() -= angle0;
		(*fi).V(1)->Q() -= angle1;
		(*fi).V(2)->Q() -= angle2;
	}
	i=0;
	for(vi=m.vert.begin(); vi!=m.vert.end(); ++vi,++i) if(!(*vi).IsD())
	{
		(*vi).Q() /= area[i];
		(*vi).Q()=math::Clamp((*vi).Q(),-0.050f,0.050f);
		
 /*   if ( (*vi).Q() < 0 )
			(*vi).Q() = log( -(*vi).Q() );
		else if( (*vi).Q() > 0 )
			(*vi).Q() = log( (*vi).Q() );*/

	}
	
	//--- DeInit
	
	delete[] area;

}
};
}
}
#endif
