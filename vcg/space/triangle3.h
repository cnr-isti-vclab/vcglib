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
#ifndef __VCG_TRIANGLE3
#define __VCG_TRIANGLE3

#include <vcg/space/box3.h>

namespace vcg {

/** \addtogroup space */
/*@{*/
/** 
		Templated class for storing a generic triangle in a 3D space.
    Note the relation with the Face class of TriMesh complex, both classes provide the P(i) access functions to their points and therefore they share the algorithms on it (e.g. area, normal etc...)
 */
template <class ScalarTriangleType> class Triangle3
{
public:
  typedef ScalarTriangleType ScalarType;
	typedef Point3< ScalarType > CoordType;
	/// The bounding box type
	typedef Box3<ScalarType> BoxType;

/*********************************************
    blah
    blah
**/

protected:
	/// Vector of vertex pointer incident in the face
	Point3<ScalarType> _v[3];
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


/** Calcola i coefficienti della combinazione convessa.
	@param bq Punto appartenente alla faccia
	@param a Valore di ritorno per il vertice V(0)
	@param b Valore di ritorno per il vertice V(1)
	@param _c Valore di ritorno per il vertice V(2)
	@return true se bq appartiene alla faccia, false altrimenti
*/
bool InterpolationParameters(const CoordType & bq, ScalarType &a, ScalarType &b, ScalarType &_c ) const
{	
const ScalarType EPSILON = ScalarType(0.000001);


#define x1 (cP(0).x())
#define y1 (cP(0).y())
#define z1 (cP(0).z())
#define x2 (cP(1).x())
#define y2 (cP(1).y())
#define z2 (cP(1).z())
#define x3 (cP(2).x())
#define y3 (cP(2).y())
#define z3 (cP(2).z())
#define px (bq.x())
#define py (bq.y())
#define pz (bq.z())

     ScalarType t1  = px*y2;
     ScalarType t2  = px*y3;
     ScalarType t3  = py*x2;
     ScalarType t4  = py*x3;
     ScalarType t5  = x2*y3;
     ScalarType t6  = x3*y2;
     ScalarType t8  = x1*y2;
     ScalarType t9  = x1*y3;
     ScalarType t10 = y1*x2;
     ScalarType t11 = y1*x3;
     ScalarType t13 = t8-t9-t10+t11+t5-t6;
     if(fabs(t13)>=EPSILON)
	 {
         ScalarType t15 = px*y1;
         ScalarType t16 = py*x1;
         a =  (t1 -t2-t3 +t4+t5-t6 )/t13;
         b = -(t15-t2-t16+t4+t9-t11)/t13;
         _c =  (t15-t1-t16+t3+t8-t10)/t13;
		return true;
     }

     t1  = px*z2;
     t2  = px*z3;
     t3  = pz*x2;
     t4  = pz*x3;
     t5  = x2*z3;
     t6  = x3*z2;
     t8  = x1*z2;
     t9  = x1*z3;
     t10 = z1*x2;
     t11 = z1*x3;
     t13 = t8-t9-t10+t11+t5-t6;
     if(fabs(t13)>=EPSILON)
	 {
		ScalarType t15 = px*z1;
		ScalarType t16 = pz*x1;
		a =  (t1 -t2-t3 +t4+t5-t6 )/t13;
		b = -(t15-t2-t16+t4+t9-t11)/t13;
		_c =  (t15-t1-t16+t3+t8-t10)/t13;
		return true;
     }

     t1  = pz*y2; t2  = pz*y3;
     t3  = py*z2; t4  = py*z3;
     t5  = z2*y3; t6  = z3*y2;
     t8  = z1*y2; t9  = z1*y3;
     t10 = y1*z2; t11 = y1*z3;
     t13 = t8-t9-t10+t11+t5-t6;
     if(fabs(t13)>=EPSILON)
	 {
         ScalarType t15 = pz*y1;
         ScalarType t16 = py*z1;
         a =  (t1 -t2-t3 +t4+t5-t6 )/t13;
         b = -(t15-t2-t16+t4+t9-t11)/t13;
         _c =  (t15-t1-t16+t3+t8-t10)/t13;
		return true;
     }
	 
#undef x1
#undef y1
#undef z1
#undef x2
#undef y2
#undef z2
#undef x3
#undef y3
#undef z3
#undef px
#undef py
#undef pz

     return false;
}




/// Return the _q of the face, the return value is in [0,sqrt(3)/2] = [0 - 0.866.. ]
ScalarType QualityFace( ) const
{
	
	return Quality(P(0), P(1), P(2));
	/*
	CoordType d10 = P(1) - P(0);
	CoordType d20 = P(2) - P(0);
	CoordType d12 = P(1) - P(2);

	CoordType x = d10^d20;

	ScalarType a = Norm( x );		// doppio dell' Area
	ScalarType b;
	
	b = Norm2( d10 );
	ScalarType t = b; 
	t = Norm2( d20 ); if( b<t ) b = t;
	t = Norm2( d12 ); if( b<t ) b = t;

	assert(b!=0.0);

	return a/b;*/

}




}; //end Class

/// Compute a shape quality measure of the triangle composed by points p0,p1,p2
/// It Returns 2*AreaTri/(MaxEdge^2), 
/// the range is range [0.0, 0.866] 
/// e.g. Equilateral triangle sqrt(3)/2, halfsquare: 1/2, ... up to a line that has zero quality.
template<class P3ScalarType>
P3ScalarType Quality( Point3<P3ScalarType> const &p0, Point3<P3ScalarType> const & p1,  Point3<P3ScalarType> const & p2)
{
	Point3<P3ScalarType> d10=p1-p0;
	Point3<P3ScalarType> d20=p2-p0;
	Point3<P3ScalarType> d12=p1-p2;
	Point3<P3ScalarType> x = d10^d20;

	P3ScalarType a = Norm( x );
	if(a==0) return 0; // Area zero triangles have surely quality==0;
	P3ScalarType b = SquaredNorm( d10 );
	P3ScalarType t = b;
	t = SquaredNorm( d20 ); if ( b<t ) b = t;
	t = SquaredNorm( d12 ); if ( b<t ) b = t;
	assert(b!=0.0);
	return a/b;
}

/// Returns the normal to the plane passing through p0,p1,p2
template<class TriangleType>
Point3<typename TriangleType::ScalarType> Normal(const TriangleType &t)
{
	return (( t.P(1) - t.P(0)) ^ (t.P(2) - t.P(0)));
}

/// Like the above, it returns the normal to the plane passing through p0,p1,p2, but normalized.
template<class TriangleType>
Point3<typename TriangleType::ScalarType> NormalizedNormal(const TriangleType &t)
{
	return (( t.P(1) - t.P(0)) ^ (t.P(2) - t.P(0))).Normalize();
}

/// Return the area of the triangle
template<class TriangleType>
typename TriangleType::ScalarType Area(const TriangleType &t) 
{
	return Norm( (t.P(1) - t.P(0)) ^ (t.P(2) - t.P(0)) );
}

template<class TriangleType>
Point3<typename TriangleType::ScalarType> Barycenter(const TriangleType &t) 
{
	return (t.P(0)+t.P(1)+t.P(2))/ScalarType(3.0);
}

template<class TriangleType>
typename TriangleType::ScalarType Perimeter(const TriangleType &t) 
{
	return Distance(t.P(0),t.P(1))+
		     Distance(t.P(1),t.P(2))+
				 Distance(t.P(2),t.P(0));
}


}	 // end namespace


#endif

