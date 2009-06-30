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
#include <vcg/space/triangle3.h>
#include <vcg/space/point2.h>
#include <vcg/space/segment2.h>
#include <float.h>

namespace vcg {

/** \addtogroup space */
/*@{*/
/** 
		Templated class for storing a generic triangle in a 2D space.
    Note the relation with the Face class of TriMesh complex, both classes provide the P(i) access functions to their points and therefore they share the algorithms on it (e.g. area, normal etc...)
 */
template <class SCALAR_TYPE> class Triangle2
{
public:
  typedef SCALAR_TYPE ScalarType;
  typedef Point2< ScalarType > CoordType;
  typedef Triangle2<ScalarType> TriangleType;

protected:
	/// Vector of vertex pointer incident in the face
	Point2<ScalarType> _v[3];
public:

	Triangle2()
	{}

	Triangle2(const CoordType &p0,const CoordType &p1,const CoordType &p2)
	{
		P(0)=p0;
		P(1)=p1;
		P(2)=p2;
	}

	/// Shortcut per accedere ai punti delle facce
	inline CoordType & P( const int j ) { return _v[j];}
	inline CoordType & P0( const int j ) { return _v[j];}
	inline CoordType & P1( const int j ) { return _v[(j+1)%3];}
	inline CoordType & P2( const int j ) { return _v[(j+2)%3];}
	inline const CoordType &  P( const int j ) const { return _v[j];}
	inline const CoordType &  P0( const int j ) const { return _v[j];}
	inline const CoordType &  P1( const int j ) const { return _v[(j+1)%3];}
	inline const CoordType &  P2( const int j ) const { return _v[(j+2)%3];}
	inline const CoordType & cP0( const int j ) const { return _v[j];}
	inline const CoordType & cP1( const int j ) const { return _v[(j+1)%3];}
	inline const CoordType & cP2( const int j ) const { return _v[(j+2)%3];}

/** evaluate barycentric coordinates
	@param bq Point on the face
	@param a barycentric value for V(0)
	@param b barycentric value for V(1)
	@param c barycentric value for V(2)
	@return true se bq appartain to the face, false otherwise
*/
bool InterpolationParameters(const CoordType & bq, ScalarType &a, ScalarType &b, ScalarType &c ) const
{	
	const ScalarType EPSILON = ScalarType(0.0001f);

	ScalarType AreaGlobal=(P(1) - P(0)) ^ (P(2) - P(0));
	ScalarType Area0=((P(2) - P(1)) ^ (bq - P(1)));
	ScalarType Area1=((P(0) - P(2)) ^ (bq - P(2)));
	ScalarType Area2=((P(1) - P(0)) ^ (bq - P(0)));
	//ScalarType AreaGlobal=Area0+Area1+Area2;
	/*if ((Area0>(AreaGlobal+EPSILON))||(Area1>(AreaGlobal+EPSILON))||(Area2>(AreaGlobal+EPSILON)))
		return false;*/
	a=Area0/AreaGlobal;
	b=Area1/AreaGlobal;
	c=Area2/AreaGlobal;

	///test inside/outside
	if(((a>(ScalarType)1+EPSILON)||(b>(ScalarType)1+EPSILON)||(c>(ScalarType)1+EPSILON))||
	  ((a<-EPSILON)||(b<-EPSILON)||(c<-EPSILON)))
		return false;

	///approximation errors
	if(a>1)
		a=(ScalarType)1;
	if(b>1)
		b=(ScalarType)1;
	if(c>1)
		c=(ScalarType)1;
	if(a<0)
		a=(ScalarType)0;
	if(b<0)
		b=(ScalarType)0;
	if(c<0)
		c=(ScalarType)0;

	
	return true;
}

///return the distance to the point q and neighors point p
void PointDistance(const CoordType & q,
				    ScalarType & dist, 
				    CoordType & p ) const
{
	dist=FLT_MAX;
	///find distance to each segment and take minimum
	for (int i=0;i<3;i++)
	{
		vcg::Segment2<float> s=vcg::Segment2<float>(P(i),P((i+1)%3));
		CoordType clos=ClosestPoint<ScalarType>(s,q);
		ScalarType dis_test=(clos-q).Norm();
		if (dis_test<dist)
		{
			dist=dis_test;
			p=clos;
		}
	}
}

}; //end Class


}	 // end namespace
#endif

