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
Revision 1.21  2007/12/02 07:39:19  cignoni
disambiguated sqrt call

Revision 1.20  2007/11/26 14:11:38  ponchio
Added Mean Ratio metric for triangle quality.

Revision 1.19  2007/11/19 17:04:05  ponchio
QualityRadii values fixed.

Revision 1.18  2007/11/18 19:12:54  ponchio
Typo (missing comma).

Revision 1.17  2007/11/16 14:22:35  ponchio
Added qualityRadii: computes inradius /circumradius.
(ok the name is ugly...)

Revision 1.16  2007/10/10 15:11:30  ponchio
Added Circumcenter function.

Revision 1.15  2007/05/10 09:31:15  cignoni
Corrected InterpolationParameters invocation

Revision 1.14  2007/05/04 16:33:27  ganovelli
moved InterpolationParamaters out the class Triangle

Revision 1.13  2007/04/04 23:23:55  pietroni
- corrected and renamed distance to point ( function TrianglePointDistance)

Revision 1.12  2007/01/13 00:25:23  cignoni
Added (Normalized) Normal version templated on three points (instead forcing the creation of a new triangle)

Revision 1.11  2006/10/17 06:51:33  fiorin
In function Barycenter, replaced calls to (the inexistent) cP(i) with P(i)

Revision 1.10  2006/10/10 09:33:47  cignoni
added quality for triangle wrap

Revision 1.9  2006/09/14 08:44:07  ganovelli
changed t.P(*) in t.cP() nella funzione Barycenter

Revision 1.8  2006/06/01 08:38:58  pietroni
added PointDistance function

Revision 1.7  2006/03/01 15:35:09  pietroni
compiled InterspolationParameters function

Revision 1.6  2006/01/22 10:00:56  cignoni
Very Important Change: Area->DoubleArea (and no more Area function)

Revision 1.5  2005/09/23 14:18:27  ganovelli
added constructor

Revision 1.4  2005/04/14 11:35:09  ponchio
*** empty log message ***

Revision 1.3  2004/07/15 13:22:37  cignoni
Added the standard P() access function instead of the shortcut P0()

Revision 1.2  2004/07/15 10:17:42  pietroni
correct access to point funtions call in usage of triangle3 (ex. t.P(0)  in t.P0(0))

Revision 1.1  2004/03/08 01:13:31  cignoni
Initial commit


****************************************************************************/
#ifndef __VCG_TRIANGLE3
#define __VCG_TRIANGLE3

#include <vcg/space/box3.h>
#include <vcg/space/point2.h>
#include <vcg/space/point3.h>
#include <vcg/space/plane3.h>
#include <vcg/space/segment3.h>
#include <vcg/space/triangle2.h>

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
	Triangle3(){}
	Triangle3(const CoordType & c0,const CoordType & c1,const CoordType & c2){_v[0]=c0;_v[1]=c1;_v[2]=c2;}
protected:
	/// Vector of vertex pointer incident in the face
	Point3<ScalarType> _v[3];
public:

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

	bool InterpolationParameters(const CoordType & bq, ScalarType &a, ScalarType &b, ScalarType &_c ) const{
		return InterpolationParameters(*this, bq, a, b,_c ); 
	}

/// Return the _q of the face, the return value is in [0,sqrt(3)/2] = [0 - 0.866.. ]
ScalarType QualityFace( ) const
{
	return Quality(P(0), P(1), P(2));
}

}; //end Class

/// Returns the normal to the plane passing through p0,p1,p2
template<class TriangleType>
typename TriangleType::ScalarType QualityFace(const TriangleType &t)
{
  return Quality(t.cP(0), t.cP(1), t.cP(2));
}

// More robust function to computing barycentric coords of a point inside a triangle.
// it requires the knowledge of what is the direction that is more orthogonal to the face plane. 
// Usually this info can be stored in a bit of the face flags (see updateFlags::FaceProjection(MeshType &m) )
// and accessing the field with
//        if(fp->Flags() & FaceType::NORMX )   axis = 0;
//        else if(fp->Flags() & FaceType::NORMY )   axis = 1;
//        else axis =2;
//        InterpolationParameters(*fp,axis,Point,Bary);
// This direction is used to project the triangle in 2D and solve the problem in 2D where it is well defined.

template<class TriangleType, class ScalarType>
bool InterpolationParameters(const TriangleType t, const int Axis, const Point3<ScalarType> & P,  Point3<ScalarType> & L)
{
  Point2<ScalarType> test;
	typedef Point2<ScalarType> P2;
	if(Axis==0) return InterpolationParameters2( P2(t.P(0)[1],t.P(0)[2]), P2(t.P(1)[1],t.P(1)[2]), P2(t.P(2)[1],t.P(2)[2]), P2(P[1],P[2]), L);
	if(Axis==1) return InterpolationParameters2( P2(t.P(0)[0],t.P(0)[2]), P2(t.P(1)[0],t.P(1)[2]), P2(t.P(2)[0],t.P(2)[2]), P2(P[0],P[2]), L);
	if(Axis==2) return InterpolationParameters2( P2(t.P(0)[0],t.P(0)[1]), P2(t.P(1)[0],t.P(1)[1]), P2(t.P(2)[0],t.P(2)[1]), P2(P[0],P[1]), L);
	return false; 
}
/// Handy Wrapper of the above one that uses the passed normal N to choose the right orientation
template<class TriangleType, class ScalarType>
bool InterpolationParameters(const TriangleType t, const Point3<ScalarType> & N, const Point3<ScalarType> & P,  Point3<ScalarType> & L)
{
	if(N[0]>N[1])
	{
		if(N[0]>N[2])
			return InterpolationParameters(t,0,P,L); /* 0 > 1 ? 2 */
		else 
			return InterpolationParameters(t,2,P,L); /* 2 > 1 ? 2 */
		}
	else 
	{
		if(N[1]>N[2])
			return InterpolationParameters(t,1,P,L); /* 1 > 0 ? 2 */
		else 
			return InterpolationParameters(t,2,P,L); /* 2 > 1 ? 2 */
	}
}

// Function that computes the barycentric coords of a 2D triangle. Used by the above function. 
// Algorithm: simply find a base for the frame of the triangle, assuming v3 as origin (matrix T) invert it and apply to P-v3.

template<class ScalarType>
bool InterpolationParameters2(const Point2<ScalarType> &V1,
													   const Point2<ScalarType> &V2,
														 const Point2<ScalarType> &V3,
														 const Point2<ScalarType> &P, Point3<ScalarType> &L)
{
	vcg::Triangle2<ScalarType> t2=vcg::Triangle2<ScalarType>(V1,V2,V3);
	return (t2.InterpolationParameters(P,L.X(),L.Y(),L.Z() ));
}

//// Function that computes the barycentric coords of a 2D triangle. Used by the above function. 
//// Algorithm: simply find a base for the frame of the triangle, assuming v3 as origin (matrix T) invert it and apply to P-v3.
//
//template<class ScalarType>
//bool InterpolationParameters2(const Point2<ScalarType> &V1,
//													   const Point2<ScalarType> &V2,
//														 const Point2<ScalarType> &V3,
//														 const Point2<ScalarType> &P, Point3<ScalarType> &L)
//{
//	ScalarType T00 = V1[0]-V3[0];  ScalarType T01 = V2[0]-V3[0];  
//	ScalarType T10 = V1[1]-V3[1];  ScalarType T11 = V2[1]-V3[1]; 
//	ScalarType Det = T00 * T11 - T01*T10;
//	if(fabs(Det) < 0.0000001) 
//				return false;
//	
//	ScalarType IT00 =  T11/Det;	ScalarType IT01 = -T01/Det;	
//	ScalarType IT10 = -T10/Det;	ScalarType IT11 =  T00/Det;	
//
//	Point2<ScalarType> Delta = P-V3;
//	
//	L[0] = IT00*Delta[0] + IT01*Delta[1];
//	L[1] = IT10*Delta[0] + IT11*Delta[1]; 
// 
//	if(L[0]<0) L[0]=0;
//	if(L[1]<0) L[1]=0;
//	if(L[0]>1.) L[0]=1;
//	if(L[1]>1.) L[1]=1;
//	
//	L[2] = 1. - L[1] - L[0];
//	if(L[2]<0) L[2]=0;
//	
//	assert(L[2] >= -0.00001);
//
//	return true;
//}


/** Calcola i coefficienti della combinazione convessa.
	@param bq Punto appartenente alla faccia
	@param a Valore di ritorno per il vertice V(0)
	@param b Valore di ritorno per il vertice V(1)
	@param _c Valore di ritorno per il vertice V(2)
	@return true se bq appartiene alla faccia, false altrimenti
*/
template<class TriangleType, class ScalarType>
bool InterpolationParameters(const TriangleType t,const Point3<ScalarType> & N,const Point3<ScalarType> & bq, ScalarType &a, ScalarType &b, ScalarType &c ) 
{	
	Point3<ScalarType> bary;
	bool done= InterpolationParameters(t,N,bq,bary);
	a=bary[0];
	b=bary[1];
	c=bary[2];
	return done;

}


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
  if(b==0) return 0; // Again: area zero triangles have surely quality==0;
	P3ScalarType t = b;
	t = SquaredNorm( d20 ); if ( b<t ) b = t;
	t = SquaredNorm( d12 ); if ( b<t ) b = t;
	return a/b;
}


/// Compute a shape quality measure of the triangle composed by points p0,p1,p2
/// It Returns inradius/circumradius
/// the range is range [0, 1] 
/// e.g. Equilateral triangle 1, halfsquare: 0.81, ... up to a line that has zero quality.
template<class P3ScalarType>
P3ScalarType QualityRadii(Point3<P3ScalarType> const &p0,
													Point3<P3ScalarType> const &p1,
													Point3<P3ScalarType> const &p2) {

	P3ScalarType a=(p1-p0).Norm();
	P3ScalarType b=(p2-p0).Norm();
	P3ScalarType c=(p1-p2).Norm();

	P3ScalarType sum = (a + b + c)*0.5;
	P3ScalarType area2 =  sum*(a+b-sum)*(a+c-sum)*(b+c-sum);
	if(area2 <= 0) return 0;
	//circumradius: (a*b*c)/(4*sqrt(area2))
	//inradius: (a*b*c)/(4*circumradius*sum) => sqrt(area2)/sum;
	return (8*area2)/(a*b*c*sum);
}

/// Compute a shape quality measure of the triangle composed by points p0,p1,p2
/// It Returns mean ratio 2sqrt(a, b)/(a+b) where a+b are the eigenvalues of the M^tM of the
/// transformation matrix into a regular simplex
/// the range is range [0, 1] 
template<class P3ScalarType>
P3ScalarType QualityMeanRatio(Point3<P3ScalarType> const &p0,
													Point3<P3ScalarType> const &p1,
													Point3<P3ScalarType> const &p2) {

	P3ScalarType a=(p1-p0).Norm();
	P3ScalarType b=(p2-p0).Norm();
	P3ScalarType c=(p1-p2).Norm();
	P3ScalarType sum = (a + b + c)*0.5; //semiperimeter
	P3ScalarType area2 =  sum*(a+b-sum)*(a+c-sum)*(b+c-sum);
	if(area2 <= 0) return 0;
	return (4.0*sqrt(3.0)*sqrt(area2))/(a*a + b*b + c*c);
}

/// Returns the normal to the plane passing through p0,p1,p2
template<class TriangleType>
Point3<typename TriangleType::ScalarType> Normal(const TriangleType &t)
{
	return (( t.P(1) - t.P(0)) ^ (t.P(2) - t.P(0)));
}
template<class Point3Type>
Point3Type Normal( Point3Type const &p0, Point3Type const & p1,  Point3Type const & p2)
{
	return (( p1 - p0) ^ (p2 - p0));
}


/// Like the above, it returns the normal to the plane passing through p0,p1,p2, but normalized.
template<class TriangleType>
Point3<typename TriangleType::ScalarType> NormalizedNormal(const TriangleType &t)
{
	return (( t.P(1) - t.P(0)) ^ (t.P(2) - t.P(0))).Normalize();
}
template<class Point3Type>
Point3Type NormalizedNormal( Point3Type const &p0, Point3Type const & p1,  Point3Type const & p2)
{
	return (( p1 - p0) ^ (p2 - p0)).Normalize();
}


/// Handy Wrapper of the above one that calculate the normal on the triangle
template<class TriangleType, class ScalarType>
bool InterpolationParameters(const TriangleType t, const Point3<ScalarType> & P,  Point3<ScalarType> & L)
{
    vcg::Point3<ScalarType> N=vcg::Normal<TriangleType>(t);
        return (InterpolationParameters<TriangleType,ScalarType>(t,N,P,L));
}

/// Return the Double of area of the triangle
// NOTE the old Area function has been removed to intentionally 
// cause compiling error that will help people to check their code...
// A some  people used Area assuming that it returns the double and some not. 
// So please check your codes!!!
// And please DO NOT Insert any Area named function here!

template<class TriangleType>
typename TriangleType::ScalarType DoubleArea(const TriangleType &t) 
{
	return Norm( (t.P(1) - t.P(0)) ^ (t.P(2) - t.P(0)) );
}

template<class TriangleType>
typename TriangleType::ScalarType CosWedge(const TriangleType &t, int k)
{
  typename TriangleType::CoordType 
    e0 = t.P((k+1)%3) - t.P(k),
    e1 = t.P((k+2)%3) - t.P(k);
  return (e0*e1)/(e0.Norm()*e1.Norm());
}

template<class TriangleType>
Point3<typename TriangleType::ScalarType> Barycenter(const TriangleType &t) 
{
	return ((t.P(0)+t.P(1)+t.P(2))/(typename TriangleType::ScalarType) 3.0);
}

template<class TriangleType>
typename TriangleType::ScalarType Perimeter(const TriangleType &t) 
{
	return Distance(t.P(0),t.P(1))+
		     Distance(t.P(1),t.P(2))+
				 Distance(t.P(2),t.P(0));
}

template<class TriangleType>
Point3<typename TriangleType::ScalarType> Circumcenter(const TriangleType &t)
{
   typename TriangleType::ScalarType a2 = (t.P(1) - t.P(2)).SquaredNorm();
   typename TriangleType::ScalarType b2 = (t.P(2) - t.P(0)).SquaredNorm();
   typename TriangleType::ScalarType c2 = (t.P(0) - t.P(1)).SquaredNorm();      
   Point3<typename TriangleType::ScalarType>c = t.P(0)*a2*(-a2 + b2 + c2) + 
                                                t.P(1)*b2*( a2 - b2 + c2) + 
                                                t.P(2)*c2*( a2 + b2 - c2);
   c /= 2*(a2*b2 + a2*c2 + b2*c2) - a2*a2 - b2*b2 - c2*c2;
   return c;
}


}	 // end namespace


#endif

