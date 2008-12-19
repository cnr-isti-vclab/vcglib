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
Revision 1.11  2006/10/07 10:02:16  cignoni
Added missing typename for interp.parameters

Revision 1.10  2005/11/30 14:05:04  ponchio
Fixed some UberZ fuynctions and non defined _flags

Revision 1.9  2005/10/14 12:34:55  cignoni
Added ordered constructor that build a edge with unique ordering
among vertices (useful for edge-collapse simplification)

Revision 1.8  2005/10/01 09:22:51  cignoni
Major rewriting of the whole class edge. Removed default flags and nonsense attibutes. Given consistent naming to defines.

Revision 1.7  2005/07/15 15:45:51  ganovelli
template parametere Scalar removed

Revision 1.6  2005/04/14 11:35:09  ponchio
*** empty log message ***

Revision 1.5  2004/10/25 16:25:12  ponchio
inline Set(...)  -> inline void Set(...)

Revision 1.4  2004/10/25 08:21:17  ganovelli
added: constructor,Set and some minor changes.

Revision 1.3  2004/05/10 14:40:28  ganovelli
name of adhacency function updated

Revision 1.2  2004/05/10 14:02:29  ganovelli
created

Revision 1.1  2004/04/26 19:04:23  ganovelli
created

****************************************************************************/
#ifndef __VCGLIB__EDGE_TYPE_BASE
#define __VCGLIB__EDGE_TYPE_BASE

#pragma message("[VCGLIB Warning]  this way to define the simplex edge is DEPRECATED  and no more SUPPORTED") 
#pragma message("[VCGLIB Warning]  use vcg/simplex/edgeplus instead ") 


#include <vcg/space/box3.h>
#include <vcg/space/texcoord2.h>

namespace vcg {

/**
\ingroup segment
    @name segment
		Class Edge.
    This is the base class for definition of a face of the mesh.
		@param SVTYPE (Templete Parameter) Specifies the vertex class type.
 */
template <class EDGENAME, class SVTYPE, class TCTYPE = TexCoord2<float,1> > class EDGE_TYPE
{
public:
	///	The base type of the segment
	typedef EDGE_TYPE BaseEdgeType;
	///	The scalar type derived from the vertex
	typedef typename SVTYPE::ScalarType ScalarType;
	/// The vertex type
	typedef SVTYPE VertexType;
	/// The type of the the vertex coordinate
	typedef Point3< ScalarType > CoordType;
	
	/// The bounding box type
	typedef Box3<ScalarType> BoxType;
	
  /// Default Empty Costructor
  inline EDGE_TYPE(){}

  inline EDGE_TYPE(VertexType* v0,VertexType* v1){v[0]=v0;v[1]=v1;}

  static inline EDGE_TYPE OrderedEdge(VertexType* v0,VertexType* v1){
   if(v0<v1) return EDGE_TYPE(v0,v1);
   else return EDGE_TYPE(v1,v0);
  }



	/// Costructor
  inline void Set(VertexType* v0,VertexType* v1){v[0]=v0;v[1]=v1;}


/***********************************************/
/** @name Vertex Pointer
    blah
    blah
**/
  //@{
protected:
	/// Vector of vertex pointer incident in the face
	VertexType *v[2];
public:
	/** Return the pointer to the j-th vertex of the face.
		@param j Index of the face vertex.
	 */
	inline SVTYPE * & V( const int j )
	{	
		assert( !IsD() );
		assert(j >= 0 && j <  2);
		return v[j];
	}

	inline const SVTYPE * const & V( const int j ) const
	{
		assert( !IsD() );
		assert(j>=0 && j<2);
		return v[j];
	}
	inline const SVTYPE * const & cV( const int j ) const
	{
		assert( !IsD() );
		assert(j>=0 && j<2);
		return v[j];
	}

	// Shortcut per accedere ai punti delle facce
	inline CoordType & P( const int j )
	{	
		assert( !IsD() );
		assert(j>=0 && j<2);
		return v[j]->P();
	}

	inline const CoordType & P( const int j ) const
	{
		assert( !IsD() );
		assert(j>=0 && j<2);
		return v[j]->cP();
	}
	inline const CoordType & cP( const int j ) const
	{
		assert( !IsD() );
		assert(j>=0 && j<2);
		return v[j]->cP();
	}

	/** Return the pointer to the ((j+1)%3)-th vertex of the face.
		@param j Index of the face vertex.
	 */
	inline SVTYPE * & V0( const int j ) { return V(j);}
	inline SVTYPE * & V1( const int j ) { return V((j+1)%2);}
	inline const SVTYPE * const &  V0( const int j ) const { return V(j);}
	inline const SVTYPE * const &  V1( const int j ) const { return V((j+1)%2);}
	inline const SVTYPE * const & cV0( const int j ) const { return cV(j);}
	inline const SVTYPE * const & cV1( const int j ) const { return cV((j+1)%2);}

	/// Shortcut per accedere ai punti delle facce
	inline CoordType & P0( const int j ) { return V(j)->P();}
	inline CoordType & P1( const int j ) { return V((j+1)%2)->P();}
	inline const CoordType &  P0( const int j ) const { return V(j)->P();}
	inline const CoordType &  P1( const int j ) const { return V((j+1)%2)->P();}
	inline const CoordType & cP0( const int j ) const { return cV(j)->P();}
	inline const CoordType & cP1( const int j ) const { return cV((j+1)%2)->P();}

	inline SVTYPE * & UberV( const int j )
	{	
		assert(j>=0 && j<2);
		return v[j];
	}

	inline const SVTYPE * const & UberV( const int j ) const
	{
		assert(j>=0 && j<2);
		return v[j];
	}


  //@}

/***********************************************/
/** @name Normal
    blah
    blah
**/
  //@{

#ifdef __VCGLIB_EDGE_EN
	/// This vector indicates the normal of the face (defines if FACE_N is defined)
protected:
	CoordType _n;
public:
#endif

  /// Return the reference of the normal to the face (if __VCGLIB_EDGE_FN is defined).
	inline CoordType & N()
	{
#ifdef __VCGLIB_EDGE_EN
	return _n;
#else
	assert(0);
	return *(CoordType *)0;
#endif
	}
		/// Return the reference of the normal to the face (if __VCGLIB_EDGE_FN is defined).
	inline const CoordType & N() const
	{
#ifdef __VCGLIB_EDGE_EN
		return _n;
#else
	return *(CoordType *)0;
#endif
	}
	/// Return the reference of the normal to the face (if __VCGLIB_EDGE_FN is defined).
	inline const CoordType cN() const
	{
#ifdef __VCGLIB_EDGE_EN
		return _n;
#else
	return *(CoordType *)0;
#endif
	}

  //@}

/***********************************************/
/** @name Quality
    blah
    blah
**/
  //@{

#ifdef __VCGLIB_EDGE_EQ
protected:
	float _q;
#endif
public:
	float & Q()
	{
#ifdef __VCGLIB_EDGE_EQ
		return _q;
#else
		assert(0);
		return *(float*)(0);
#endif
	}

const float & Q() const
	{
#ifdef __VCGLIB_EDGE_EQ
		return _q;
#else
		assert(0);
		return *(float*)(0);
#endif
	}

  //@}

/***********************************************/
/** @name Colors
    blah
    blah
**/
  //@{
protected:
#ifdef __VCGLIB_EDGE_EC
	Color4b _c;
#endif

public:
	Color4b & C()
	{
#ifdef __VCGLIB_EDGE_EC
		return _c;
#else
		assert(0);
		return *(Color4b*)(0);
#endif
	}

	const Color4b C() const
	{
#ifdef __VCGLIB_EDGE_EC
		return _c;
#else
		return Color4b(Color4b::White);
#endif
	}



  //@}

/***********************************************/
/** @name Adjacency
    blah
    blah
**/
  //@{

protected:
#if defined(__VCGLIB_EDGE_AE)
  /// Vector of face pointer, it's used to indicate the adjacency relations (defines if FACE_A is defined)
	EDGENAME   *ee[2];				// edge adiacenti
	/// Index of the face in the arrival face 
	char zs[2];									
#endif

#ifdef __VCGLIB_EDGE_AV
	///Vettore di puntatori a edge, utilizzato per indicare le adiacenze vertice faccia
	EDGENAME *ev[2];
	char zv[2];
#endif

public:




	/** Return the pointer to the j-th adjacent edge.
	    @param j Index of the edge.
	 */
	inline EDGENAME * & EEp( const int j )
	{
		assert( !IsD() );
		assert(j>=0 && j<2);
#if defined(__VCGLIB_EDGE_AE)
		  return ee[j];
#else 
		assert(0);
    return *(EDGENAME **)(0);;
#endif
	}

	inline const EDGENAME * const & EEp( const int j ) const
	{
		assert( !IsD() );
		assert(j>=0 && j<2);
#if defined(__VCGLIB_EDGE_AE)
		  return ee[j];
#else
		  assert(0);
		  return (EDGENAME *)0;
#endif
	}
	inline EDGENAME * & EEp1( const int j ) { return EEp((j+1)%2);}
	inline const EDGENAME * const&  EEp1( const int j ) const { return EEp((j+1)%2);}

/** Return the pointer to the j-th adjacent face.
	    @param j Index of the edge.
	 */
	inline EDGENAME * & UberEEp( const int j )
	{
		assert(j>=0 && j<2);
#if defined(__VCGLIB_EDGE_AE)
		return ee[j];
#else 
		assert(0); // if you stop here you are probably trying to use FF topology in a face without it
    return *(EDGENAME **)(0);
#endif
	}

	inline const EDGENAME * const & UberEEp( const int j ) const
	{
		assert(j>=0 && j<2);
#if defined(__VCGLIB_EDGE_AE)
		return ee[j];
#else
		assert(0); // if you stop here you are probably trying to use FF topology in a face without it
    return *(EDGENAME **)(0);
#endif
	}
	

	inline EDGENAME * & VEp( const int j )
	{
		assert( !IsD() );
		assert(j>=0 && j<2);
#ifdef __VCGLIB_EDGE_AV
		return ev[j];
#else
		assert(0); // you are probably trying to use VF topology in a vertex without it
    return *(EDGENAME **)(0);
#endif
	}

	inline const EDGENAME * const & VEp( const int j ) const
	{
		assert( !IsD() );
		assert(j>=0 && j<2);
#ifdef __VCGLIB_EDGE_AV
		return ev[j];
#else
		assert(0);
    return *(EDGENAME **)(0);
#endif
	}


	/** Return the index that the face have in the j-th adjacent face.
	    @param j Index of the edge.
	 */
	inline char & EEi( const int j )
	{
		assert( !IsD() );
		assert(j>=0 && j<2);
#if defined(__VCGLIB_EDGE_AE) 
		return zs[j];
#else
		assert(0);
		return *(char *)0; // tanto per farlo compilare...
#endif
	}

	inline const char & EEi( const int j ) const
	{
		assert( !IsD() );
		assert(j>=0 && j<2);
#if defined(__VCGLIB_EDGE_AE) 
		return zs[j];
#else
		assert(0);
		return *(char *)0;
#endif
	}

		/** Return the index that the face have in the j-th adjacent face.
	    @param j Index of the edge.
	 */
	inline char & UberZ( const int j )
	{
		assert(j>=0 && j<2);
#if defined(__VCGLIB_EDGE_AE) 
		return zs[j];
#elif defined(__VCGLIB_EDGE_SA) 
		return zs[j];
#else
		assert(0);
		static char dummy = 0;
		return dummy;
#endif
	}

	inline const char & UberZ( const int j ) const
	{
		assert(j>=0 & j<2);
#if defined(__VCGLIB_EDGE_AE) 
		return zs[j];
#elif defined(__VCGLIB_EDGE_SA) 
		return zs[j];
#else
		assert(0);
        static int dummy = 0;
        return dummy;
#endif
	}


	inline char & VEi( const int j )
	{
		assert( !IsD() );
		assert(j>=0 & j<2);
#ifdef __VCGLIB_EDGE_VA
		return zv[j];
#elif defined(__VCGLIB_EDGE_SA)
		return zs[j];
#else
		assert(0);
		static char dummy = 0;
        return dummy;
#endif
	}

	inline const char & VEi( const int j ) const
	{
		assert( !IsD() );
		assert(j>=0 & j<2);
#ifdef __VCGLIB_EDGE_VA
		return zv[j];
#elif defined(__VCGLIB_EDGE_SA)
		return zs[j];
#else
		assert(0);
		static char dummy = 0;
        return dummy;
#endif
	}

  //@}

/***********************************************/
/** @name Mark
    blah
    blah
**/
  //@{


#ifdef __VCGLIB_EDGE_EM
	/// Incremental mark (defines if FACE_I is defined)
	int imark;
#endif // Mark

	inline int & IMark()
	{
#ifdef __VCGLIB_EDGE_EM
		assert( !IsD() );
		assert( (_flags & NOTREAD) == 0 );
		assert( (_flags & NOTWRITE) == 0 );
		return imark;
#else
    return 0;
#endif // Mark
	}

	inline const int & IMark() const
	{
		assert( !IsD() );
#ifdef __VCGLIB_EDGE_EM		
		assert( (_flags & NOTREAD) == 0 );
		return imark;
#else
        static int dummy = 0;
        return dummy;		
#endif				
	}

	/// Initialize the imark system of the face
	inline void InitIMark()
	{
#ifdef __VCGLIB_EDGE_EM
		imark = 0;
#endif
	}

  
  //@}
/***********************************************/
/** @name Flags
    blah
    blah
**/
  //@{

	/// This are the _flags of face, the default value is 0
#ifdef __VCGLIB_EDGE_EF
	int  _flags;		
#endif
	enum {
		// This bit indicate that the face is deleted from the mesh
		DELETED     = 0x00000001,		// cancellato
		// This bit indicate that the face of the mesh is not readable
		NOTREAD     = 0x00000002,		// non leggibile (ma forse modificabile)
		// This bit indicate that the face is not modifiable
		NOTWRITE    = 0x00000004,		// non modificabile (ma forse leggibile) 
		// This bit indicate that the face is modified
		SELECTED    = 0x00000020,		// Selection _flags
		// Border _flags, it is assumed that BORDERi = BORDER0<<i 
		BORDER0     = 0x00000040,
		BORDER1     = 0x00000080,
		// First user bit
		USER0       = 0x00040000
		};
public:
	static int &LastBitFlag()
		{
			static int b =USER0;
			return b;
		}
	static inline int NewBitFlag()
		{
			LastBitFlag()=LastBitFlag()<<1;
			return LastBitFlag();
		}
	static inline bool DeleteBitFlag(int bitval)
		{	
			if(LastBitFlag()==bitval) {
					LastBitFlag()= LastBitFlag()>>1;
					return true;
			}
			assert(0);
			return false;
		}

  void ClearFlags() {
#ifdef __VCGLIB_EDGE_EF
  _flags=0;
#endif
  }

	/// Return the _flags.
	inline int & Flags ()
	{
#ifdef __VCGLIB_EDGE_EF
		assert( !IsD() );
		return _flags;
#else
    return *(int *)0;
#endif
	}

	inline const int & Flags () const
	{
#ifdef __VCGLIB_EDGE_EF
		assert( !IsD() );
		return _flags;
#else
    return 0;
#endif
	}
	/// Ritorna il _flags senza effettuare alcun controllo sui relativi bit
	inline int & UberFlags()
	{
#ifdef __VCGLIB_EDGE_EF
		return _flags;
#else
    assert(0);
    return *(int *)0;
#endif
	}

	inline const int UberFlags() const
	{
#ifdef __VCGLIB_EDGE_EF
		return _flags;
#else
    return 0;
#endif
	}

	/// This function checks if the face is deleted 
	bool IsD() const {
#ifdef __VCGLIB_EDGE_EF
    return (_flags & DELETED) != 0;
#else
    return false;
#endif
  }
  /// This function mark the face as deleted
  void SetD()		{
#ifdef __VCGLIB_EDGE_EF
    _flags |=DELETED;
#endif
  }
	/// This function mark the face as not deleted
	void ClearD()	{
#ifdef __VCGLIB_EDGE_EF
   _flags &= (~DELETED);
#endif
  }
	
	
	/// This function checks if the face is selected
	bool IsS() const {
#ifdef __VCGLIB_EDGE_EF
    return (_flags & SELECTED) != 0;
#else
    return false;
#endif
  }
	/// This function select the face
	void SetS()		{
#ifdef __VCGLIB_EDGE_EF
    _flags |=SELECTED;
#endif
  }
	/// This funcion execute the inverse operation of SetS()
	void ClearS()	{
#ifdef __VCGLIB_EDGE_EF
    _flags &= (~SELECTED);
#endif
  }

	/// This function checks if the edge is Border on a given side
	bool IsB(int i) const {
#ifdef __VCGLIB_EDGE_EF
    return (_flags & (BORDER0<<i)) != 0;
#else
    return false;
#endif
  }
	/// This function set edge as Border on a given side
	void SetB(int i)		{
#ifdef __VCGLIB_EDGE_EF
    _flags |=(BORDER0<<i);
#endif
  }
	/// This function clear edge as Border on a given side
	void ClearB(int i)	{
#ifdef __VCGLIB_EDGE_EF
    _flags &= (~(BORDER0<<i));
#endif
  }
	
	/// This function checks if the given user bit is true
	bool IsUserBit(int userBit){
#ifdef __VCGLIB_EDGE_EF
    return (_flags & userBit) != 0;
#else
    return false;
#endif
    }
	/// This function set  the given user bit 
	void SetUserBit(int userBit){
#ifdef __VCGLIB_EDGE_EF
    _flags |=userBit;
#endif
  }
	/// This function clear the given user bit 
	void ClearUserBit(int userBit){
#ifdef __VCGLIB_EDGE_EF
    _flags &= (~userBit);
#endif
  }


//@}
/*#*******************	
*  Bounding box *
**********************/

void GetBBox( BoxType & bb )
{
	bb.Set( v[0]->P() );
	bb.Add( v[1]->P() );
}

 /***********************************************/
 /** @name Reflection Functions 
 Static functions  that give information about the current vertex type.
Reflection is a mechanism making it possible to investigate yourself. Reflection is used to investigate format of objects at runtime, invoke methods and access fields of these objects. Here we provide static const functions that are resolved at compile time and they give information about the data (normal, color etc.) supported by the current vertex type.
 **/
 //@{

static bool HasEdgeNormal()  { 
#ifdef __VCGLIB_EDGE_FN 
  return true;
#else
  return false;
#endif
}
static bool HasEdgeQuality()  { 
#ifdef __VCGLIB_EDGE_FQ
  return true;
#else
  return false;
#endif
}
static bool HasEdgeColor()  { 
#ifdef __VCGLIB_EDGE_FC 
  return true;
#else
  return false;
#endif
}
static bool HasEEAdjacency()  { 
#if (defined(__VCGLIB_EDGE_AE) )
  return true;
#else
  return false;
#endif
}
static bool HasVEAdjacency()  { 
#if (defined(__VCGLIB_EDGE_AV) )
  return true;
#else
  return false;
#endif
}
static bool HasEdgeMark()  { 
#ifdef __VCGLIB_EDGE_FC 
  return true;
#else
  return false;
#endif
}

//@}

  /// operator to compare two edges
	inline bool operator == ( const EDGENAME & f ) const {
		if( (V(0) != f.V(0)) && (V(0) != f.V(1)) ) return false;
		if( (V(1) != f.V(0)) && (V(1) != f.V(1)) ) return false;
		return true;
	}

/** Calcola i coefficienti della combinazione convessa.
	@param bq Punto appartenente alla faccia
	@param a Valore di ritorno per il vertice V(0)
	@param b Valore di ritorno per il vertice V(1)
	@param _c Valore di ritorno per il vertice V(2)
	@return true se bq appartiene alla faccia, false altrimenti
*/
  bool InterpolationParameters(const CoordType & bq, typename VertexType::ScalarType &a, ScalarType &_b) const
{	
  typedef typename VertexType::ScalarType ScalarType;
const ScalarType EPSILON = ScalarType(0.000001);
ScalarType l;

#define x1 (cV(0)->P().x())
#define y1 (cV(0)->P().y())
#define z1 (cV(0)->P().z())
#define x2 (cV(1)->P().x())
#define y2 (cV(1)->P().y())
#define z2 (cV(1)->P().z())
#define px (bq.x())
#define py (bq.y())
#define pz (bq.z())
				a = (px-x1)/(x2-x1);
				l = (py-y1)/(y2-y1);
				if(	 (  l < a -EPSILON) ||  (  l > a +EPSILON))
					return false;

				l = (pz-z1)/(z2-z1);
				if(	 (  l < a -EPSILON) ||  (  l > a +EPSILON))
					return false;

				_b = 1-a;
				return true;

#undef x1
#undef y1
#undef z1
#undef x2
#undef y2
#undef z2
#undef px
#undef py
#undef pz
}



/// Return the DOUBLE of the area of the face
ScalarType Length() const
{
	return Norm( (V(1)->P() - V(0)->P()).Norm());
}

CoordType Barycenter() const
{
	return (V(0)->P()+V(1)->P())/ScalarType(2.0);
}

}; //end Class




}	 // end namespace

#endif



