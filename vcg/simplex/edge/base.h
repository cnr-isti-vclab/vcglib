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

#include <vcg/space/box3.h>
#include <vcg/space/tcoord2.h>

namespace vcg {

/**
\ingroup segment
    @name segment
		Class Edge.
    This is the base class for definition of a face of the mesh.
		@param SVTYPE (Templete Parameter) Specifies the vertex class type.
 */
template <class SVTYPE, class TCTYPE = TCoord2<float,1> > class EDGE_TYPE
{
public:
	///	The base type of the segment
	typedef EDGE_TYPE BaseEdgeType;
	/// The vertex type
	typedef SVTYPE VertexType;
	/// The type of the scalar field of the vertex coordinate
  typedef typename VertexType::ScalarType ScalarType;
	/// The type of the the vertex coordinate
	typedef Point3< ScalarType > CoordType;
	typedef Point3< ScalarType > NormalType;
	
  typedef typename SVTYPE::EdgeType EdgeFromVertType;
	/// The bounding box type
	typedef Box3<ScalarType> BoxType;
	
  /// Default Empty Costructor
  inline EDGE_TYPE(){}

	/// This are the _flags of face, the default value is 0
	int  _flags;		

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
		assert( (_flags & DELETED) == 0 );
		assert( (_flags & NOTREAD) == 0 ); 
		assert( (_flags & NOTWRITE) == 0 );
		assert(j >= 0);
		assert(j <  2);
		return v[j];
	}

	inline const SVTYPE * const & V( const int j ) const
	{
		assert( (_flags & DELETED) == 0 );
		assert( (_flags & NOTREAD) == 0 );
		assert(j>=0);
		assert(j<2);
		return v[j];
	}
	inline const SVTYPE * const & cV( const int j ) const
	{
		assert( (_flags & DELETED) == 0 );
		assert( (_flags & NOTREAD) == 0 );
		assert(j>=0);
		assert(j<2);
		return v[j];
	}

	// Shortcut per accedere ai punti delle facce
	inline CoordType & P( const int j )
	{	
		assert( (_flags & DELETED) == 0 );
		assert( (_flags & NOTREAD) == 0 ); 
		assert( (_flags & NOTWRITE) == 0 );
		assert(j>=0);
		assert(j<2);
		return v[j]->P();
	}

	inline const CoordType & P( const int j ) const
	{
		assert( (_flags & DELETED) == 0 );
		assert( (_flags & NOTREAD) == 0 );
		assert(j>=0);
		assert(j<2);
		return v[j]->cP();
	}
	inline const CoordType & cP( const int j ) const
	{
		assert( (_flags & DELETED) == 0 );
		assert( (_flags & NOTREAD) == 0 );
		assert(j>=0);
		assert(j<2);
		return v[j]->cP();
	}

	/** Return the pointer to the ((j+1)%3)-th vertex of the face.
		@param j Index of the face vertex.
	 */
	inline SVTYPE * & V0( const int j ) { return V(j);}
	inline SVTYPE * & V1( const int j ) { return V((j+1)%2);}
	inline const SVTYPE * const &  V0( const int j ) const { return V(j);}
	inline const SVTYPE * const &  V1( const int j ) const { return V((j+1)%3);}
	inline const SVTYPE * const & cV0( const int j ) const { return cV(j);}
	inline const SVTYPE * const & cV1( const int j ) const { return cV((j+1)%3);}

	/// Shortcut per accedere ai punti delle facce
	inline CoordType & P0( const int j ) { return V(j)->P();}
	inline CoordType & P1( const int j ) { return V((j+1)%3)->P();}
	inline const CoordType &  P0( const int j ) const { return V(j)->P();}
	inline const CoordType &  P1( const int j ) const { return V((j+1)%3)->P();}
	inline const CoordType & cP0( const int j ) const { return cV(j)->P();}
	inline const CoordType & cP1( const int j ) const { return cV((j+1)%3)->P();}

	inline SVTYPE * & UberV( const int j )
	{	
		assert(j>=0);
		assert(j<2);
		return v[j];
	}

	inline const SVTYPE * const & UberV( const int j ) const
	{
		assert(j>=0);
		assert(j<2);
		return v[j];
	}


  //@}

/***********************************************/
/** @name Normal
    blah
    blah
**/
  //@{

#ifdef __VCGLIB_EDGE_FN
	/// This vector indicates the normal of the face (defines if FACE_N is defined)
protected:
	CoordType _n;
public:
#endif

  /// Return the reference of the normal to the face (if __VCGLIB_EDGE_FN is defined).
	inline CoordType & N()
	{
#ifdef __VCGLIB_EDGE_FN
	return _n;
#else
	assert(0);
	return *(CoordType *)0;
#endif
	}
		/// Return the reference of the normal to the face (if __VCGLIB_EDGE_FN is defined).
	inline const CoordType & N() const
	{
#ifdef __VCGLIB_EDGE_FN
		return _n;
#else
	return *(CoordType *)0;
#endif
	}
	/// Return the reference of the normal to the face (if __VCGLIB_EDGE_FN is defined).
	inline const CoordType cN() const
	{
#ifdef __VCGLIB_EDGE_FN
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

#ifdef __VCGLIB_EDGE_FQ
protected:
	float _q;
#endif
public:
	float & Q()
	{
#ifdef __VCGLIB_EDGE_FQ
		return _q;
#else
		assert(0);
		return *(float*)(&_flags);
#endif
	}

const float & Q() const
	{
#ifdef __VCGLIB_EDGE_FQ
		return _q;
#else
		assert(0);
		return *(float*)(&_flags);
#endif
	}

  //@}

/***********************************************/
/** @name Texture
    blah
    blah
**/
  //@{

// Per Wedge Texture Coords
protected:
#ifdef __VCGLIB_EDGE_WT
	TCTYPE _wt[3];
#endif
public:
	TCTYPE & WT(const int i)
	{
#ifdef __VCGLIB_EDGE_WT
		return _wt[i];
#else
		assert(0);
		return *(TCTYPE*)(&_flags);
#endif
	}

	const TCTYPE & WT(const int i) const
	{
#ifdef __VCGLIB_EDGE_WT
		return _wt[i];
#else
		assert(0);
		return *(TCTYPE*)(&_flags);
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
#ifdef __VCGLIB_EDGE_FC
	Color4b _c;
#endif

public:
	Color4b & C()
	{
#ifdef __VCGLIB_EDGE_FC
		return _c;
#else
		assert(0);
		return *(Color4b*)(&_flags);
#endif
	}

	const Color4b C() const
	{
#ifdef __VCGLIB_EDGE_FC
		return _c;
#else
		return Color4b(Color4b::White);
#endif
	}

protected:
#ifdef __VCGLIB_EDGE_WC
	Color4b _wc[3];
#endif
public:
	Color4b & WC(const int i)
	{
#ifdef __VCGLIB_EDGE_WC
		return _wc[i];
#else
		assert(0);
		return *(Color4b*)(&_flags);
#endif
	}

const Color4b WC(const int i) const
	{
#ifdef __VCGLIB_EDGE_WC
		return _wc[i];
#else
		assert(0);
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

#if (defined(__VCGLIB_EDGE_EA) && defined(__VCGLIB_EDGE_EA))
	#error Error: You cannot specify face-to-face and shared topology together
#endif

#if (defined(__VCGLIB_EDGE_VA) && defined(__VCGLIB_EDGE_EA))
	#error Error: You cannot specify vertex-face and shared topology together
#endif

protected:
#if defined(__VCGLIB_EDGE_EA)
  /// Vector of face pointer, it's used to indicate the adjacency relations (defines if FACE_A is defined)
	EDGE_TYPE   *ss[3];				// Facce adiacenti
	/// Index of the face in the arrival face 
	char zs[4];									
#endif

#ifdef __VCGLIB_EDGE_VA
	///Vettore di puntatori a faccia, utilizzato per indicare le adiacenze vertice faccia
	EDGE_TYPE *sv[3];
	char zv[3];
#endif

#ifdef __VCGLIB_EDGE_EA
	///Vettore di puntatori a faccia, utilizzato per indicare le adiacenze vertice faccia
	EDGE_TYPE *ses[3];
	char zs[3];
#endif
public:




	/** Return the pointer to the j-th adjacent face.
	    @param j Index of the edge.
	 */
	inline EDGE_TYPE * & S( const int j )
	{
		assert( (_flags & DELETED) == 0 );
		assert( (_flags & NOTREAD) == 0 );
		assert( (_flags & NOTWRITE) == 0 );
		assert(j>=0);
	    assert(j<2);
#if defined(__VCGLIB_EDGE_EA)
		  return ss[j];
#elif defined(__VCGLIB_EDGE_EA)
			return ses[j];
#else 
		assert(0);
        static EDGE_TYPE *dum=0;
		return dum;
#endif
	}

	inline const EDGE_TYPE * const & S( const int j ) const
	{
		assert( (_flags & DELETED) == 0 );
		assert( (_flags & NOTREAD) == 0 );
		assert(j>=0);
	    assert(j<2);
#if defined(__VCGLIB_EDGE_EA)
		  return ss[j];
#elif defined(__VCGLIB_EDGE_EA)
			return ses[j];
#else
		  assert(0);
		  return (EDGE_TYPE *)this;
#endif
	}
	inline EDGE_TYPE * & S1( const int j ) { return F((j+1)%2);}
	inline const EDGE_TYPE * const&  S1( const int j ) const { return F((j+1)%2);}

/** Return the pointer to the j-th adjacent face.
	    @param j Index of the edge.
	 */
	inline EDGE_TYPE * & UberF( const int j )
	{
		assert(j>=0);
	  assert(j<2);
#if defined(__VCGLIB_EDGE_EA)
		  return ss[j];
#elif defined(__VCGLIB_EDGE_EA)
			return ses[j];
#else 
		assert(0); // if you stop here you are probably trying to use FF topology in a face without it
		return *((EDGE_TYPE **)(_flags));
#endif
	}

	inline const EDGE_TYPE * const & UberF( const int j ) const
	{
		assert(j>=0);
	  assert(j<2);
#if defined(__VCGLIB_EDGE_EA)
		  return ss[j];
#elif defined(__VCGLIB_EDGE_EA)
			return ses[j];
#else
		assert(0); // if you stop here you are probably trying to use FF topology in a face without it
		return *((EDGE_TYPE **)(_flags));
#endif
	}
	

	inline EDGE_TYPE * & Fv( const int j )
	{
		assert( (_flags & DELETED) == 0 );
		assert( (_flags & NOTREAD) == 0 );
		assert( (_flags & NOTWRITE) == 0 );
		assert(j>=0);
		assert(j<2);
#ifdef __VCGLIB_EDGE_VA
		return sv[j];
#elif defined(__VCGLIB_EDGE_EA)
		return ses[j];
#else
		assert(0); // you are probably trying to use VF topology in a vertex without it
		return *((EDGE_TYPE **)(_flags));
#endif
	}

	inline const EDGE_TYPE * const & Fv( const int j ) const
	{
		assert( (_flags & DELETED) == 0 );
		assert( (_flags & NOTREAD) == 0 );
		assert(j>=0);
		assert(j<2);
#ifdef __VCGLIB_EDGE_VA
		return sv[j];
#elif defined(__VCGLIB_EDGE_EA)
		return ses[j];
#else
		assert(0);
		return (EDGE_TYPE *)this;
#endif
	}


	/** Return the index that the face have in the j-th adjacent face.
	    @param j Index of the edge.
	 */
	inline char & Z( const int j )
	{
		assert( (_flags & DELETED) == 0 );
		assert( (_flags & NOTREAD) == 0 );
		assert( (_flags & NOTWRITE) == 0 );
		assert(j>=0);
		assert(j<2);
#if defined(__VCGLIB_EDGE_EA) 
		return zs[j];
#elif defined(__VCGLIB_EDGE_EA) 
		return zs[j];
#else
		assert(0);
		return *(char *)&_flags; // tanto per farlo compilare...
#endif
	}

	inline const char & Z( const int j ) const
	{
		assert( (_flags & DELETED) == 0 );
		assert( (_flags & NOTREAD) == 0 );
		assert(j>=0);
		assert(j<2);
#if defined(__VCGLIB_EDGE_EA) 
		return zs[j];
#elif defined(__VCGLIB_EDGE_EA) 
		return zs[j];
#else
		assert(0);
		return *(char *)&_flags;
#endif
	}

		/** Return the index that the face have in the j-th adjacent face.
	    @param j Index of the edge.
	 */
	inline char & UberZ( const int j )
	{
		assert(j>=0);
		assert(j<2);
#if defined(__VCGLIB_EDGE_EA) 
		return zs[j];
#elif defined(__VCGLIB_EDGE_EA) 
		return zs[j];
#else
		assert(0);
		return *(char *)&_flags;
#endif
	}

	inline const char & UberZ( const int j ) const
	{
		assert(j>=0);
		assert(j<2);
#if defined(__VCGLIB_EDGE_EA) 
		return zs[j];
#elif defined(__VCGLIB_EDGE_EA) 
		return zs[j];
#else
		assert(0);
		return *(char *)&_flags;
#endif
	}


	inline char & Zv( const int j )
	{
		assert( (_flags & DELETED) == 0 );
		assert( (_flags & NOTREAD) == 0 );
		assert( (_flags & NOTWRITE) == 0 );
		assert(j>=0);
		assert(j<2);
#ifdef __VCGLIB_EDGE_VA
		return zv[j];
#elif defined(__VCGLIB_EDGE_EA)
		return zs[j];
#else
		assert(0);
		return *(char *)&_flags;
#endif
	}

	inline const char & Zv( const int j ) const
	{
		assert( (_flags & DELETED) == 0 );
		assert( (_flags & NOTREAD) == 0 );
		assert(j>=0);
		assert(j<2);
#ifdef __VCGLIB_EDGE_VA
		return zv[j];
#elif defined(__VCGLIB_EDGE_EA)
		return zs[j];
#else
		assert(0);
		return *(char *)&_flags;
#endif
	}

  //@}

/***********************************************/
/** @name Mark
    blah
    blah
**/
  //@{


#ifdef __VCGLIB_EDGE_FM
	/// Incremental mark (defines if FACE_I is defined)
	int imark;
#endif // Mark
#ifdef __VCGLIB_EDGE_M
	inline int & IMark()
	{
		assert( (_flags & DELETED) == 0 );
		assert( (_flags & NOTREAD) == 0 );
		assert( (_flags & NOTWRITE) == 0 );
		return imark;
	}

	inline const int & IMark() const
	{
		assert( (_flags & DELETED) == 0 );
		assert( (_flags & NOTREAD) == 0 );
		return imark;
	}
#endif // Mark

	/// Initialize the imark system of the face
	inline void InitIMark()
	{
#ifdef __VCGLIB_EDGE_M
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

    void ClearFlags() {_flags=0;}

	/// Return the _flags.
	inline int & Flags ()
	{
		assert( (_flags & DELETED) == 0 );
		assert( (_flags & NOTREAD) == 0 );
		return _flags;
	}

	inline const int & Flags () const
	{
		assert( (_flags & DELETED) == 0 );
		assert( (_flags & NOTREAD) == 0 );
		return _flags;
	}
	/// Ritorna il _flags senza effettuare alcun controllo sui relativi bit
	inline int & UberFlags()
	{
		return _flags;
	}

	inline const int UberFlags() const
	{
		return _flags;
	}

	/// This function checks if the face is deleted 
	bool IsD() const {return (_flags & DELETED) != 0;}
	/// This function mark the face as deleted
	void SetD()		{_flags |=DELETED;}
	/// This function mark the face as not deleted
	void ClearD()	{_flags &= (~DELETED);}
	/// This function checks if the face is deleted 
	bool IsDeleted() const {return IsD();}
	
	/// This function checks if the face is readable 
	bool IsR() const {return (_flags & NOTREAD) == 0;}
	/// This function marks the face as readable
	void SetR()		{_flags &= (~NOTREAD);}
	/// This function marks the face as not readable
	void ClearR() {_flags |=NOTREAD;}
	
	/// This function checks if the face is readable 
	bool IsW() const {return (_flags & NOTWRITE)== 0;}
	/// This function marks the vertex as not writable
	void SetW() {_flags &=(~NOTWRITE);}
	/// This function marks the face as not writable
	void ClearW() {_flags |=NOTWRITE;}

	/// This funcion checks whether the face is both readable and modifiable
	bool IsRW() const {return (_flags & (NOTREAD | NOTWRITE)) == 0;}

	
	/// This function checks if the face is selected
	bool IsS() const {return (_flags & SELECTED) != 0;}
	/// This function select the face
	void SetS()		{_flags |=SELECTED;}
	/// This funcion execute the inverse operation of SetS()
	void ClearS()	{_flags &= (~SELECTED);}

	/// This function checks if the face is selected
	bool IsB(int i) const {return (_flags & (BORDER0<<i)) != 0;}
	/// This function select the face
	void SetB(int i)		{_flags |=(BORDER0<<i);}
	/// This funcion execute the inverse operation of SetS()
	void ClearB(int i)	{_flags &= (~(BORDER0<<i));}

	/// This function checks if the face is Crease  on side i
	bool IsEE(int i) const {return (_flags & (FEATURE0<<i)) != 0;}
	/// This function select the face flag
	void SetEE(int i)		{_flags |=(FEATURE0<<i);}
	/// This funcion execute the inverse operation of Set()
	void ClearEE(int i)	{_flags &= (~(FEATURE0<<i));}

	/// This function checks if the given user bit is true
	bool IsUserBit(int userBit){return (_flags & userBit) != 0;}
	/// This function set  the given user bit 
	void SetUserBit(int userBit){_flags |=userBit;}
	/// This function clear the given user bit 
	void ClearUserBit(int userBit){_flags &= (~userBit);}


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
#if (defined(__VCGLIB_EDGE_EA) || defined(__VCGLIB_EDGE_EA))
  return true;
#else
  return false;
#endif
}
static bool HasVSAdjacency()  { 
#if (defined(__VCGLIB_EDGE_VA) || defined(__VCGLIB_EDGE_EA))
  return true;
#else
  return false;
#endif
}
static bool HasSharedAdjacency()  { 
#if defined(__VCGLIB_EDGE_EA)
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

  /// operator to compare two faces
	inline bool operator == ( const EDGE_TYPE & f ) const {
		for(int i=0; i<3; ++i)
			if( (V(i) != f.V(0)) && (V(i) != f.V(1)) )
				return false;
		return true;
	}

/** Calcola i coefficienti della combinazione convessa.
	@param bq Punto appartenente alla faccia
	@param a Valore di ritorno per il vertice V(0)
	@param b Valore di ritorno per il vertice V(1)
	@param _c Valore di ritorno per il vertice V(2)
	@return true se bq appartiene alla faccia, false altrimenti
*/
bool InterpolationParameters(const CoordType & bq, ScalarType &a, ScalarType &_b) const
{	
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




