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
Revision 1.2  2004/02/19 15:40:56  cignoni
Added doxygen groups

Revision 1.1  2004/02/13 02:16:22  cignoni
First working release.


****************************************************************************/


#ifndef __VCGLIB_BOX3
#define __VCGLIB_BOX3

#include <vcg/space/point3.h>

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
	inline  Box3() { min.x()= 1;max.x()= -1;min.y()= 1;max.y()= -1;min.z()= 1;max.z()= -1;}
		/// Copy constructor
	inline  Box3( const Box3 & b ) { min=b.min; max=b.max; }
		/// Min Max constructor
	inline  Box3( const Point3<BoxScalarType> & mi, const Point3<BoxScalarType> & ma ) { min = mi; max = ma; }
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
	void Inflate( const BoxScalarType s )
	{
		Inflate( (max-min)*s );
	}
		/** Varia le dimensioni del bounding box di (k,k,k) con k = bbox.diag*s
		*/
	void InflateFix( const BoxScalarType s )
	{
		BoxScalarType k = Diag()*s;
		Inflate( Point3<BoxScalarType> (k,k,k));
	}
		/** Varia le dimensioni del bounding box del valore fornito attraverso il parametro.
			@param delta Point in 3D space
		*/
	void Inflate( const Point3<BoxScalarType> & delta )
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
		min.x()= 1; max.x()= -1;
		min.y()= 1; max.y()= -1;
		min.z()= 1; max.z()= -1;
	}
		/** Function to add two bounding box
			@param b Il bounding box che si vuole aggiungere
		*/
	void Add( Box3<BoxScalarType> const & b )
	{
		if(IsNull()) *this=b;
		else
		{
			if(min.x() > b.min.x()) min.x() = b.min.x();
			if(min.y() > b.min.y()) min.y() = b.min.y();
			if(min.z() > b.min.z()) min.z() = b.min.z();

			if(max.x() < b.max.x()) max.x() = b.max.x();
			if(max.y() < b.max.y()) max.y() = b.max.y();
			if(max.z() < b.max.z()) max.z() = b.max.z();
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
			if(min.x() > p.x()) min.x() = p.x();
			if(min.y() > p.y()) min.y() = p.y();
			if(min.z() > p.z()) min.z() = p.z();

			if(max.x() < p.x()) max.x() = p.x();
			if(max.y() < p.y()) max.y() = p.y();
			if(max.z() < p.z()) max.z() = p.z();
		}
	}
	//
	//// Aggiunge ad un box un altro box trasformato secondo la matrice m
	//void Add( const Matrix44<BoxScalarType> &m, const Box3<BoxScalarType> & b )
	//{
	//		const Point3<BoxScalarType> &mn= b.min;
	//		const Point3<BoxScalarType> &mx= b.max;
 //     Add(m.Apply(Point3<BoxScalarType>(mn[0],mn[1],mn[2])));
	//		Add(m.Apply(Point3<BoxScalarType>(mx[0],mn[1],mn[2])));
	//		Add(m.Apply(Point3<BoxScalarType>(mn[0],mx[1],mn[2])));
	//		Add(m.Apply(Point3<BoxScalarType>(mx[0],mx[1],mn[2])));
	//		Add(m.Apply(Point3<BoxScalarType>(mn[0],mn[1],mx[2])));
	//		Add(m.Apply(Point3<BoxScalarType>(mx[0],mn[1],mx[2])));
	//		Add(m.Apply(Point3<BoxScalarType>(mn[0],mx[1],mx[2])));
	//		Add(m.Apply(Point3<BoxScalarType>(mx[0],mx[1],mx[2])));
	//}
		/** Calcola l'intersezione tra due bounding box. Al bounding box viene assegnato il valore risultante.
			@param b Il bounding box con il quale si vuole effettuare l'intersezione
		*/
	void Intersect( const Box3<BoxScalarType> & b )
	{
		if(min.x() < b.min.x()) min.x() = b.min.x();
		if(min.y() < b.min.y()) min.y() = b.min.y();
		if(min.z() < b.min.z()) min.z() = b.min.z();

		if(max.x() > b.max.x()) max.x() = b.max.x();
		if(max.y() > b.max.y()) max.y() = b.max.y();
		if(max.z() > b.max.z()) max.z() = b.max.z();

		if(min.x()>max.x() || min.y()>max.y() || min.z()>max.z()) SetNull();
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
			min.x() <= p.x() && p.x() <= max.x() &&
			min.y() <= p.y() && p.y() <= max.y() &&
			min.z() <= p.z() && p.z() <= max.z()
		);
	}
		/** Verifica se un punto appartiene ad un bounding box aperto sul max.
			@param p The point 3D
			@return True se p appartiene al bounding box, false altrimenti
		*/
	bool IsInEx( Point3<BoxScalarType> const & p ) const
	{
		return (
			min.x() <= p.x() && p.x() < max.x() &&
			min.y() <= p.y() && p.y() < max.y() &&
			min.z() <= p.z() && p.z() < max.z()
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
	bool Collide(Box3<BoxScalarType> const &b)
	{
		return b.min.x()<max.x() && b.max.x()>min.x() &&
			   b.min.y()<max.y() && b.max.y()>min.y() &&
			   b.min.z()<max.z() && b.max.z()>min.z() ;
	}
		/** Controlla se il bounding box e' nullo.
			@return True se il bounding box e' nullo, false altrimenti
		*/
	bool IsNull() const { return min.x()>max.x() || min.y()>max.y() || min.z()>max.z(); }
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
		return (max.x()-min.x())*(max.y()-min.y())*(max.z()-min.z());
	}
		/// Calcola la dimensione del bounding box sulla x.
	inline BoxScalarType DimX() const { return max.x()-min.x();}
		/// Calcola la dimensione del bounding box sulla y.
	inline BoxScalarType DimY() const { return max.y()-min.y();}
		/// Calcola la dimensione del bounding box sulla z.
	inline BoxScalarType DimZ() const { return max.z()-min.z();}

	template <class Q>
	inline void Import( const Box3<Q> & b )
	{
		min.Import(b.min);
		max.Import(b.max);
	}

	template <class Q>
	inline Box3 Construct( const Box3<Q> & b )
	{
    return Box3(Point3<ScalarType>::Construct(b.min),Point3<ScalarType>::Construct(b.max));
	}

}; // end class definition




typedef Box3<short>  Box3s;
typedef Box3<int>	 Box3i;
typedef Box3<float>  Box3f;
typedef Box3<double> Box3d;


/*@}*/

} // end namespace
#endif
