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

#ifndef VCG_PARAM_DISTORSION
#define VCG_PARAM_DISTORSION
#include <vcg/complex/algorithms/parametrization/uv_utils.h>

namespace vcg {
	namespace tri{
		template <class MeshType>
		class Distorsion
		{
			typedef typename MeshType::FaceType FaceType;
			typedef typename MeshType::VertexType VertexType;
			typedef typename MeshType::CoordType CoordType;
			typedef typename MeshType::ScalarType ScalarType;

			static ScalarType Area3D(const FaceType *f)
			{
//				CoordType vp0=f->P(0);
//				CoordType vp1=f->P(1);
//				CoordType vp2=f->P(2);
//				ScalarType Area3D=((vp2-vp0)^(vp1-vp0)).Norm()/2.0;
				return DoubleArea(*f)*(0.5);
			}

			static ScalarType AreaUV(const FaceType *f)
			{
				vcg::Point2<ScalarType> uv0=f->V(0)->T().P();
				vcg::Point2<ScalarType> uv1=f->V(1)->T().P();
				vcg::Point2<ScalarType> uv2=f->V(2)->T().P();
				ScalarType AreaUV=((uv1-uv0)^(uv2-uv0))/2.0;
				return AreaUV;
			}

			static ScalarType EdgeLenght3D(FaceType *f,int e)
			{
				assert((e>=0)&&(e<3));
				ScalarType lenght=(f->P0(e)-f->P1(e)).Norm();
				return (lenght);
			}

			static ScalarType EdgeLenghtUV(FaceType *f,int e)
			{
				assert((e>=0)&&(e<3));
				vcg::Point2<ScalarType> uv0=f->V(e)->T().P();
				vcg::Point2<ScalarType> uv1=f->V((e+1)%3)->T().P();
				ScalarType UVlenght=(uv0-uv1).Norm();
				return (UVlenght);
			}

			static ScalarType Angle3D(const FaceType *f,int e)
			{
				assert((e>=0)&&(e<3));
				CoordType p0=f->P((e+2)%3);
				CoordType p1=f->P(e);
				CoordType p2=f->P((e+1)%3);
				typedef typename CoordType::ScalarType ScalarType;
				CoordType dir0=p2-p1;
				CoordType dir1=p0-p1;
				dir0.Normalize();
				dir1.Normalize();
				ScalarType angle=dir0*dir1;
				return angle;
			}

			static ScalarType AngleUV(const FaceType *f,int e)
			{
				vcg::Point2<ScalarType> uv0=f->V((e+2)%3)->T().P();
				vcg::Point2<ScalarType> uv1=f->V(e)->T().P();
				vcg::Point2<ScalarType> uv2=f->V((e+1)%3)->T().P();
				vcg::Point2<ScalarType> dir0=uv2-uv1;
				vcg::Point2<ScalarType> dir1=uv0-uv1;
				dir0.Normalize();
				dir1.Normalize();
				ScalarType angle=dir0*dir1;
				return angle;
			}

		public:

			///return the variance of angle, normalized
			///in absolute value
			static ScalarType AngleDistorsion(const FaceType *f,int e)
			{
				ScalarType Angle_3D=Angle3D(f,e);
				ScalarType Angle_UV=AngleUV(f,e);
				ScalarType diff=fabs(Angle_3D-Angle_UV);///Angle_3D;
				return diff;
			}
			
			///return the variance of angle, normalized
			///in absolute value
			static ScalarType AngleDistorsion(const FaceType *f)
			{
				ScalarType angleDist=0;
				for (int i=0;i<3;i++)
					angleDist+=AngleDistorsion(f,i);
				return angleDist;
			}

			///return the global scaling factor from 3D to UV
			static ScalarType ScalingFactor(MeshType &m,
				ScalarType &AreaScale,
				ScalarType &EdgeScale)
			{
				ScalarType SumArea3D=0;
				ScalarType SumArea2D=0;
				ScalarType SumEdge3D=0;
				ScalarType SumEdge2D=0;
				for (int i=0;i<m.face.size();i++)
				{
					SumArea3D+=Area3D(&m.face[i]);
					SumArea2D+=AreaUV(&m.face[i]);
					for (int j=0;j<3;j++)
					{
						SumEdge3D+=EdgeLenght3D(&m.face[i],j);
						SumEdge2D+=EdgeLenghtUV(&m.face[i],j);
					}
				}
				AreaScale=SumArea3D/SumArea2D;
				EdgeScale=SumEdge3D/SumEdge2D;
			}

			///return the variance of edge lenght, normalized
			///in absolute value, the scalar EdgeScaleVal may be calculated
			///by using the ScalingFactor function
			static ScalarType EdgeDistorsion(FaceType *f,int e,
				ScalarType EdgeScaleVal)
			{
				ScalarType edgeUV=EdgeLenghtUV(f,e)*EdgeScaleVal;
				ScalarType edge3D=EdgeLenght3D(f,e);
				ScalarType diff=fabs(edge3D-edgeUV)/edge3D;
				return diff;
			}

			///return the variance of area, normalized
			///in absolute value, the scalar AreaScaleVal may be calculated
			///by using the ScalingFactor function
//			static ScalarType AreaDistorsion(FaceType *f,
//											 ScalarType AreaScaleVal)
//			{
//			  ScalarType areaUV=AreaUV(f)*AreaScaleVal;
//			  ScalarType area3D=EdgeLenght3D(f,e);
//			  ScalarType diff=fabs(edge3D-edgeUV)/edge3D;
//			  return diff;
//			}

			///return the number of folded faces
			static bool Folded(const FaceType *f)
			{
				ScalarType areaUV=AreaUV(f);
				/*if (areaUV<0)
					printf("area %5.5f \n",areaUV);*/
				return (areaUV<0);
			}

			static int Folded(const MeshType &m)
			{
				int folded=0;
				for (int i=0;i<m.face.size();i++)
				{
					if (m.face[i].IsD())continue;
					if(Folded(&m.face[i]))folded++;
				}
				return folded;
			}

			static bool GloballyUnFolded(const MeshType &m)
			{
				int num=Folded(m);
				return (num>(m.fn)/2);
			}

			static ScalarType AngleDistorsion(const MeshType &m)
			{
				ScalarType UDdist=0;
				for (int i=0;i<m.face.size();i++)
				{
					if (m.face[i].IsD())continue;
					const FaceType *f=&(m.face[i]);
					UDdist+=AngleDistorsion(f)*Area3D(f);
				}
				return UDdist;
			}
		};
	}
}
#endif
