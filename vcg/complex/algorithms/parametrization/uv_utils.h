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
} //End Namespace Tri
} // End Namespace vcg
#endif