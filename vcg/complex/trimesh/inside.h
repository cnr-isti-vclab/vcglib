
/****************************************************************************
* IDOLib                                                                    *
* Interactive Deformable Objects Library									*
*	http://idolib.sf.net													*	
*																			*
* Copyright(C) 2005															*
* Visual Computing Lab                                                      *
* ISTI - Italian National Research Council                                  *
*                                                                           *
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
Revision 1.2  2007/06/06 14:26:51  pietroni
compiling error resolved



****************************************************************************/
#include <vcg/space/ray3.h>
#include <vcg/space/box3.h>
#include <vcg/space/triangle3.h>

#ifndef VCG_INSIDE
#define VCG_INSIDE

///this static funtion is used to see if one point is inside 
///a triangular mesh or not... First parameter is a spatial indexing 
///structure (eg. a grid) used to perform research operation, 
///initialized with faces of the triangular mesh of type TriMeshType
namespace vcg {
	namespace trimesh {
		template <class FaceSpatialIndexing,class TriMeshType>
		class Inside{

		private:

			typedef typename FaceSpatialIndexing::CoordType CoordType;
			typedef typename FaceSpatialIndexing::ScalarType ScalarType;

		public:



			//return true if the point is inside the mesh 
			static bool Is_Inside(TriMeshType &m,
				FaceSpatialIndexing &_g_mesh,
				const CoordType &test)
			{

				typedef typename TriMeshType::FaceType FaceType;
				typedef typename TriMeshType::ScalarType ScalarType;
				typedef typename TriMeshType::CoordType CoordType;

				const ScalarType EPSILON=0.000001;

				///first test if the element is inside the bb of the mesh
				if (!m.bbox.IsIn(test))
					return false;
				else
				{
					ScalarType dist;
					CoordType Norm,ip,nearest;
					FaceType *f=vcg::trimesh::GetClosestFace<TriMeshType,FaceSpatialIndexing>(m,_g_mesh,test,m.bbox.Diag(),dist,nearest,Norm,ip);
					//any face in the mesh
					assert(f!=NULL);
					//if the point is on the face is considered inside
					if ((test-nearest).Norm()<=EPSILON)
						return true;
					///if the point is in surface of a face then is ok
					if ((ip.V(0)>EPSILON)&&(ip.V(1)>EPSILON)&&(ip.V(2)>EPSILON))
					{
						if ((f->N()*(test-nearest))<0)///in this case normal direction is enought
							return true;
						else
							return false;
					}

					///in this case we are not sure because
					///hit an edge or a vertex then we use a ray that
					///go until the barycenter of found face,
					///then see normal value again
					else
					{
						CoordType bary = vcg::Barycenter< FaceType >(*f);
						vcg::Ray3<ScalarType> r;
						r.Set(test,(bary-test));///set origin and direction
						r.Normalize();
						FaceType *f1=vcg::trimesh::DoRay<TriMeshType,FaceSpatialIndexing>(m,_g_mesh,r,m.bbox.Diag(),dist);
						if (f1==NULL)
							assert(0);
						if ((f1->N()*(test-bary))<0)///in this case normal direction is enought
							return true;
						else
							return false;
					}

				}
			}

		}; // end class
	}
}

#endif