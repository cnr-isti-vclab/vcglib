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

****************************************************************************/

#ifndef __VCG_TRIANGLE2
#define __VCG_TRIANGLE2

#include <vcg\space\point2.h>

namespace vcg {

/** \addtogroup space */
/*@{*/
/** 
		Templated class for storing a generic triangle in a 2D space.
    Note the relation with the Face class of TriMesh complex, both classes provide the P(i) access functions to their points and therefore they share the algorithms on it (e.g. area, normal etc...)
 */
template <class SCALAR_TRIANGLE_TYPE> class Triangle2
{
public:
  typedef SCALAR_TRIANGLE_TYPE ScalarTriangleType;
  typedef ScalarTriangleType ScalarType;
	typedef Point2< ScalarType > CoordType;
	

protected:
	/// Vector of vertex pointer incident in the face
	Point2<ScalarType> _v[3];
public:

	/// Shortcut per accedere ai punti delle facce
	inline CoordType & P0( const int j ) { return _v[j];}
	inline CoordType & P1( const int j ) { return _v[(j+1)%3];}
	inline CoordType & P2( const int j ) { return _v[(j+2)%3];}
	inline const CoordType &  P0( const int j ) const { return _v[j];}
	inline const CoordType &  P1( const int j ) const { return _v[(j+1)%3];}
	inline const CoordType &  P2( const int j ) const { return _v[(j+2)%3];}
	inline const CoordType & cP0( const int j ) const { return _v[j];}
	inline const CoordType & cP1( const int j ) const { return _v[(j+1)%3];}
	inline const CoordType & cP2( const int j ) const { return _v[(j+2)%3];}


}; //end Class

}	 // end namespace
#endif

