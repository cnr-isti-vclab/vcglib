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

****************************************************************************/
#ifndef __VCG_TRI_UPDATE_FLAGS
#define __VCG_TRI_UPDATE_FLAGS

namespace vcg {
namespace tri {

template <class UpdateMeshType>
class UpdateFlags
{

public:
typedef UpdateMeshType MeshType; 
typedef typename MeshType::VertexType     VertexType;
typedef typename MeshType::VertexPointer  VertexPointer;
typedef typename MeshType::VertexIterator VertexIterator;
typedef typename MeshType::FaceType       FaceType;
typedef typename MeshType::FacePointer    FacePointer;
typedef typename MeshType::FaceIterator   FaceIterator;


void BorderFF(MeshType &m)
{
	const int BORDERFLAG[3]={FaceType::BORDER0,FaceType::BORDER1,FaceType::BORDER2};
	FaceIterator fi;
	for(fi=face.begin();fi!=face.end();++fi)if(!(*fi).IsD())
		for(int j=0;j<3;++j)
		{
			if(!(*fi).IsManifold(j)) (*fi).SetCF(j);
			else if((*fi).IsBorder(j)) (*fi).SetB(j);
					 else (*fi).ClearB(j);
		}
}



}; // end class

}	// End namespace
}	// End namespace


#endif
