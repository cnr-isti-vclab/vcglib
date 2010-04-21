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
Revision 1.19  2008/03/05 11:48:08  pietroni
correct 1 error on DistancePoint3Box3 (if the point is inside the box return distance to the nearest face instead of zero)

Revision 1.18  2008/02/03 23:50:51  cignoni
Important Change. Now adding a null bbox to a bbox leave it unchanged (instead of trashing it)

Revision 1.17  2007/07/12 06:41:24  cignoni
added a missing static to the Construct() member

Revision 1.16  2007/03/22 18:12:00  pietroni
*** empty log message ***

Revision 1.15  2007/03/22 17:29:16  pietroni
*** empty log message ***

Revision 1.14  2006/11/13 13:03:45  ponchio
Added GetBBox in Point3 (declaration) the body of the function is in box3.h

Revision 1.13  2006/09/28 13:37:13  m_di_benedetto
Added "const" to Collide()

Revision 1.12  2006/07/06 12:42:10  ganovelli
tolto il tipo enumerato axis

Revision 1.11  2006/06/29 13:22:04  ganovelli
Aggiunte funzioni che determinano il lato magiore (minore)

Revision 1.10  2006/06/09 07:26:37  m_di_benedetto
Added const qualifier to P() method.

Revision 1.9  2005/09/30 15:03:57  cignoni
Added Center-Radius Constructor

Revision 1.8  2004/07/20 23:10:01  cignoni
Uncommented Add of matrix transformed bbox

Revision 1.7  2004/07/15 09:47:55  ganovelli
added function P(int i) to list the box's vertices

Revision 1.6  2004/07/07 23:26:25  cignoni
removed the infamous Inflate. Now only Offset exists

Revision 1.5  2004/03/05 17:51:28  tarini
Errorino "ScalarType" -> "BoxScalarType"

Revision 1.4  2004/03/03 14:32:13  cignoni
Yet another cr lf mismatch

Revision 1.3  2004/02/23 23:44:21  cignoni
cr lf mismatch

Revision 1.2  2004/02/19 15:40:56  cignoni
Added doxygen groups

Revision 1.1  2004/02/13 02:16:22  cignoni
First working release.


****************************************************************************/


#ifndef __VCGLIB_BOX3
#define __VCGLIB_BOX3

#include <vcg/space/point3.h>
#include <vcg/math/matrix44.h>
#include <vcg/space/line3.h>
#include <vcg/space/plane3.h>

namespace vcg {

/** \addtogroup space */
/*@{*/
/** 
Templated class for 3D boxes.
  This is the class for definition of a axis aligned bounding box in 3D space. It is stored just as two Point3
	@param BoxScalarType (template parameter) Specifies the type of scalar used to represent coords.
*/
template <class BoxScalarType>
class Box3
{
public:

	/// The scalar type
	typedef BoxScalarType ScalarType;

	/// min coordinate point
  Point3<BoxScalarType> min;
	/// max coordinate point
	Point3<BoxScalarType> max;
		/// The bounding box constructor
	inline  Box3() { min.X()= 1;max.X()= -1;min.Y()= 1;max.Y()= -1;min.Z()= 1;max.Z()= -1;}
		/// Copy constructor
	inline  Box3( const Box3 & b ) { min=b.min; max=b.max; }
		/// Min Max constructor
	inline  Box3( const Point3<BoxScalarType> & mi, const Point3<BoxScalarType> & ma ) { min = mi; max = ma; }
    /// Point Radius Constructor
  inline Box3(const Point3<BoxScalarType> & center, const BoxScalarType & radius) {
    min = center-Point3<BoxScalarType>(radius,radius,radius);
    max = center+Point3<BoxScalarType>(radius,radius,radius);
  }
		/// The bounding box distructor
	inline ~Box3() { }
		/// Operator to compare two bounding box
	inline bool operator == ( Box3<BoxScalarType> const & p ) const
	{
		return min==p.min && max==p.max;
	}
		/// Operator to dispare two bounding box
	inline bool operator != ( Box3<BoxScalarType> const & p ) const
	{
		return min!=p.min || max!=p.max;
	}
		/** Varia le dimensioni del bounding box scalandole rispetto al parametro scalare.
			@param s Valore scalare che indica di quanto deve variare il bounding box
		*/
	void Offset( const BoxScalarType s )
	{
		Offset( Point3<BoxScalarType> (s,s,s));
	}
		/** Varia le dimensioni del bounding box del valore fornito attraverso il parametro.
			@param delta Point in 3D space
		*/
	void Offset( const Point3<BoxScalarType> & delta )
	{
		min -= delta;
		max += delta;
	}
		/// Initializing the bounding box
	void Set( const Point3<BoxScalarType> & p )
	{
		min = max = p;
	}
		/// Set the bounding box to a null value
	void SetNull()
	{
		min.X()= 1; max.X()= -1;
		min.Y()= 1; max.Y()= -1;
		min.Z()= 1; max.Z()= -1;
	}
		/** Function to add two bounding box
			@param b Il bounding box che si vuole aggiungere
		*/
	void Add( Box3<BoxScalarType> const & b )
	{
		if(b.IsNull()) return; // Adding a null bbox should do nothing
		if(IsNull()) *this=b;
		else
		{
			if(min.X() > b.min.X()) min.X() = b.min.X();
			if(min.Y() > b.min.Y()) min.Y() = b.min.Y();
			if(min.Z() > b.min.Z()) min.Z() = b.min.Z();

			if(max.X() < b.max.X()) max.X() = b.max.X();
			if(max.Y() < b.max.Y()) max.Y() = b.max.Y();
			if(max.Z() < b.max.Z()) max.Z() = b.max.Z();
		}
	}
		/** Funzione per aggiungere un punto al bounding box. Il bounding box viene modificato se il punto
			cade fuori da esso.
			@param p The point 3D
		*/
	void Add( const Point3<BoxScalarType> & p )
	{
		if(IsNull()) Set(p);
		else 
		{
			if(min.X() > p.X()) min.X() = p.X();
			if(min.Y() > p.Y()) min.Y() = p.Y();
			if(min.Z() > p.Z()) min.Z() = p.Z();

			if(max.X() < p.X()) max.X() = p.X();
			if(max.Y() < p.Y()) max.Y() = p.Y();
			if(max.Z() < p.Z()) max.Z() = p.Z();
		}
	}
	
	// Aggiunge ad un box un altro box trasformato secondo la matrice m
	void Add( const Matrix44<BoxScalarType> &m, const Box3<BoxScalarType> & b )
	{
			const Point3<BoxScalarType> &mn= b.min;
			const Point3<BoxScalarType> &mx= b.max;
      Add(m*(Point3<BoxScalarType>(mn[0],mn[1],mn[2])));
			Add(m*(Point3<BoxScalarType>(mx[0],mn[1],mn[2])));
			Add(m*(Point3<BoxScalarType>(mn[0],mx[1],mn[2])));
			Add(m*(Point3<BoxScalarType>(mx[0],mx[1],mn[2])));
			Add(m*(Point3<BoxScalarType>(mn[0],mn[1],mx[2])));
			Add(m*(Point3<BoxScalarType>(mx[0],mn[1],mx[2])));
			Add(m*(Point3<BoxScalarType>(mn[0],mx[1],mx[2])));
			Add(m*(Point3<BoxScalarType>(mx[0],mx[1],mx[2])));
	}
		/** Calcola l'intersezione tra due bounding box. Al bounding box viene assegnato il valore risultante.
			@param b Il bounding box con il quale si vuole effettuare l'intersezione
		*/
	void Intersect( const Box3<BoxScalarType> & b )
	{
		if(min.X() < b.min.X()) min.X() = b.min.X();
		if(min.Y() < b.min.Y()) min.Y() = b.min.Y();
		if(min.Z() < b.min.Z()) min.Z() = b.min.Z();

		if(max.X() > b.max.X()) max.X() = b.max.X();
		if(max.Y() > b.max.Y()) max.Y() = b.max.Y();
		if(max.Z() > b.max.Z()) max.Z() = b.max.Z();

		if(min.X()>max.X() || min.Y()>max.Y() || min.Z()>max.Z()) SetNull();
	}
		/** Trasla il bounding box di un valore definito dal parametro.
			@param p Il bounding box trasla sulla x e sulla y in base alle coordinate del parametro
		*/
	void Translate( const Point3<BoxScalarType> & p )
	{
		min += p;
		max += p;
	}
		/** Verifica se un punto appartiene ad un bounding box.
			@param p The point 3D
			@return True se p appartiene al bounding box, false altrimenti
		*/
	bool IsIn( Point3<BoxScalarType> const & p ) const
	{
		return (
			min.X() <= p.X() && p.X() <= max.X() &&
			min.Y() <= p.Y() && p.Y() <= max.Y() &&
			min.Z() <= p.Z() && p.Z() <= max.Z()
		);
	}
		/** Verifica se un punto appartiene ad un bounding box aperto sul max.
			@param p The point 3D
			@return True se p appartiene al bounding box, false altrimenti
		*/
	bool IsInEx( Point3<BoxScalarType> const & p ) const
	{
		return (
			min.X() <= p.X() && p.X() < max.X() &&
			min.Y() <= p.Y() && p.Y() < max.Y() &&
			min.Z() <= p.Z() && p.Z() < max.Z()
		);
	}
		/** Verifica se due bounding box collidono cioe' se hanno una intersezione non vuota. Per esempio
			due bounding box adiacenti non collidono.
			@param b A bounding box
			@return True se collidoo, false altrimenti
		*/
	/* old version
	bool Collide(Box3<BoxScalarType> const &b)
	{
		Box3<BoxScalarType> bb=*this;
		bb.Intersect(b);
		return bb.IsValid();
	}
	*/
	bool Collide(Box3<BoxScalarType> const &b) const
	{
		return b.min.X()<max.X() && b.max.X()>min.X() &&
			   b.min.Y()<max.Y() && b.max.Y()>min.Y() &&
			   b.min.Z()<max.Z() && b.max.Z()>min.Z() ;
	}
		/** Controlla se il bounding box e' nullo.
			@return True se il bounding box e' nullo, false altrimenti
		*/
	bool IsNull() const { return min.X()>max.X() || min.Y()>max.Y() || min.Z()>max.Z(); }
		/** Controlla se il bounding box e' vuoto.
			@return True se il bounding box e' vuoto, false altrimenti
		*/
	bool IsEmpty() const { return min==max; }
		/// Restituisce la lunghezza della diagonale del bounding box.
	BoxScalarType Diag() const
	{
		return Distance(min,max);
	}
		/// Calcola il quadrato della diagonale del bounding box.
	BoxScalarType SquaredDiag() const
	{
		return SquaredDistance(min,max);
	}
		/// Calcola il centro del bounding box.
	Point3<BoxScalarType> Center() const
	{
		return (min+max)/2;
	}
		/// Compute bounding box size.
	Point3<BoxScalarType> Dim() const
	{
		return (max-min);
	}
	  /// Returns global coords of a local point expressed in [0..1]^3
	Point3<BoxScalarType> LocalToGlobal(Point3<BoxScalarType> const & p) const{
		return Point3<BoxScalarType>( 
			min[0] + p[0]*(max[0]-min[0]), 
			min[1] + p[1]*(max[1]-min[1]),
			min[2] + p[2]*(max[2]-min[2]));
	}
	  /// Returns local coords expressed in [0..1]^3 of a point in 3D
	Point3<BoxScalarType> GlobalToLocal(Point3<BoxScalarType> const & p) const{
		return Point3<BoxScalarType>( 
		  (p[0]-min[0])/(max[0]-min[0]), 
		  (p[1]-min[1])/(max[1]-min[1]), 
		  (p[2]-min[2])/(max[2]-min[2])
			);
	}
		/// Calcola il volume del bounding box.
	BoxScalarType Volume() const
	{
		return (max.X()-min.X())*(max.Y()-min.Y())*(max.Z()-min.Z());
	}
		/// Calcola la dimensione del bounding box sulla x.
	inline BoxScalarType DimX() const { return max.X()-min.X();}
		/// Calcola la dimensione del bounding box sulla y.
	inline BoxScalarType DimY() const { return max.Y()-min.Y();}
		/// Calcola la dimensione del bounding box sulla z.
	inline BoxScalarType DimZ() const { return max.Z()-min.Z();}
		/// Calcola il lato di lunghezza maggiore  
	inline unsigned char MaxDim() const { 
		int i;
		Point3<BoxScalarType> diag = max-min;
		if(diag[0]>diag[1]) i=0; else i=1;
		return (diag[i]>diag[2])? i: 2;
	}
		/// Calcola il lato di lunghezza minore  
	inline unsigned char MinDim() const { 
		int i;
		Point3<BoxScalarType> diag =  max-min;
		if(diag[0]<diag[1]) i=0; else i=1;
		return (diag[i]<diag[2])? i: 2;
	}

	template <class Q>
	inline void Import( const Box3<Q> & b )
	{
		min.Import(b.min);
		max.Import(b.max);
	}

	template <class Q>
	static inline Box3 Construct( const Box3<Q> & b )
	{
    return Box3(Point3<BoxScalarType>::Construct(b.min),Point3<BoxScalarType>::Construct(b.max));
	}
		
		/// gives the ith box vertex in order: (x,y,z),(X,y,z),(x,Y,z),(X,Y,z),(x,y,Z),(X,y,Z),(x,Y,Z),(X,Y,Z)
	Point3<BoxScalarType> P(const int & i) const {
			return Point3<BoxScalarType>(
				min[0]+ (i%2) * DimX(),
				min[1]+ ((i / 2)%2) * DimY(),
				min[2]+ (i>3)* DimZ());
	}
}; // end class definition

template <class T> Box3<T> Point3<T>::GetBBox(Box3<T> &bb) const {
 bb.Set( *this );
}
 

template <class ScalarType> 
ScalarType DistancePoint3Box3(const Point3<ScalarType> &p,
							  const Box3<ScalarType> &b)
{
	///if fall inside return distance to a face
	if (b.IsIn(p))
	{
		const ScalarType dx = std::min<ScalarType>(b.max.X()-p.X(), p.X()-b.min.X());
		const ScalarType dy = std::min<ScalarType>(b.max.Y()-p.Y(), p.Y()-b.min.Y());
		const ScalarType dz = std::min<ScalarType>(b.max.Z()-p.Z(), p.Z()-b.min.Z());

		return std::min<ScalarType>(dx, std::min<ScalarType>(dy, dz));
	}

	{
		ScalarType sq_dist = ScalarType(0);
		for (int i=0; i<3; ++i)
		{
			ScalarType delta = ScalarType(0);
			if       (p[i] < b.min[i])  delta = p[i] - b.min[i];
			else if  (p[i] > b.max[i])  delta = p[i] - b.max[i];
			sq_dist += delta * delta;
		}

    return math::Sqrt(sq_dist);
	}
}


typedef Box3<short>  Box3s;
typedef Box3<int>	 Box3i;
typedef Box3<float>  Box3f;
typedef Box3<double> Box3d;


/*@}*/

} // end namespace
#endif

