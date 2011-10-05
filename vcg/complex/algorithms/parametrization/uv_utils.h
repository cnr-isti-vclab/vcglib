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

#ifndef VCG_UV_UTILS
#define VCG_UV_UTILS


namespace vcg {
	namespace tri{

		///calculate the BBox in UV space
		template <class MeshType>
		vcg::Box2<typename MeshType::ScalarType> PerWedgeUVBox(MeshType &m)
		{
			vcg::Box2<typename MeshType::ScalarType> UVBox;
			typename MeshType::FaceIterator fi;
			for (fi=m.face.begin();fi!=m.face.end();fi++)
			{
				if ((*fi).IsD()) continue;
				for (int i=0;i<3;i++)
				UVBox.Add((*fi).WT(i).P());
			}
			return UVBox;
		}

		///transform curvature to UV space
		template <class FaceType>
		vcg::Point2<typename FaceType::ScalarType> Coord3DtoUV(FaceType &f,
														typename FaceType::CoordType dir)
		{
			///then transform to UV
			CoordType bary3d=(f.P(0)+f.P(1)+f.P(2))/3.0;
			vcg::Point2<ScalarType> baryUV=(f.WT(0).P()+f.WT(1).P()+f.WT(2).P())/3.0;
			CoordType dir3d=bary3d+dir;
			CoordType baryCoordsUV;
			vcg::InterpolationParameters<FaceType,ScalarType>(f,dir3d,baryCoordsUV);
			vcg::Point2<ScalarType> dirUV=baryCoordsUV.X()*f.WT(0).P()+
											baryCoordsUV.Y()*f.WT(1).P()+
											baryCoordsUV.Z()*f.WT(2).P()-baryUV;
			dirUV.Normalize();
			return dirUV;
		}

} //End Namespace Tri
} // End Namespace vcg
#endif