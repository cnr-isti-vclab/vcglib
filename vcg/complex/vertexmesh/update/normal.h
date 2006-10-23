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
#ifndef __VCG_VERTEXMESH_UPDATE_NORMAL
#define __VCG_VERTEXMESH_UPDATE_NORMAL

#include <vcg/space/normal_extrapolation.h>

namespace vcg {
	namespace vertex {

		/** \addtogroup vertexmesh */
		/* @{ */

		/*!
		* This class is used to update the normals of a Vertex mesh
		*/
		template < class VERTEX_CONTAINER > 
		class UpdateNormal
		{
		public:
			typedef 					VERTEX_CONTAINER							VertexContainer;
			typedef typename	VERTEX_CONTAINER::value_type	VertexType;
			typedef typename	VertexType::ScalarType				ScalarType;
			typedef typename	VERTEX_CONTAINER::iterator		VertexIterator;

			/*!
			*/
			static void UpdateNormals(const VertexIterator &begin, const VertexIterator &end, int k)
			{
				vcg::NormalExtrapolation< VertexContainer >::ExtrapolateNormlas(begin, end, k);
			};

		}; //end of class UpdateNormal

		/*! @} */
	}; //end of namespace vertex
}; //end of namespace vcg

#endif //__VCG_VERTEXMESH_UPDATE_NORMAL
