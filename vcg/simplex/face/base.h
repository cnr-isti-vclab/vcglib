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
Revision 1.1  2004/02/10 01:11:28  cignoni
Edited Comments and GPL license

****************************************************************************/




#ifndef FACE_TYPE 
#pragma message("\nYou should never directly include this file\n")
#else
/* People should subclass his vertex class from this one... */
#include <assert.h>
#include <vcg/utility.h>
#include <vcg/Mesh/MeshPos.h>

#include <vcg/Point3.h>
#include <vcg/Line3.h>
#include <vcg/Plane3.h>
#include <vcg/Box3.h>
#include <vcg/tools/ColorUB.h>

namespace vcg {

/** @name Face
		Class Face.
    This is the base class for definition of a face of the mesh.
		@param FVTYPE (Templete Parameter) Specifies the vertex class type.
 */
template <class FVTYPE, class TCTYPE = TCoord<float,1> > class FACE_TYPE
{
public:
	typedef typename FVTYPE::face_type face_from_vert_type;
	typedef typename FVTYPE::scalar_type FCTYPE;
	/// The type of the scalar field of the vertex coordinate
	typedef FCTYPE scalar_type;
	/// The type of the the vertex coordinate
	typedef Point3< FCTYPE > vectorial_type;
	///	The base type of the face
	typedef FACE_TYPE face_base;
	/// The vertex type
	typedef FVTYPE vertex_type;
	/// The bounding box type
	typedef Box3<scalar_type> bbox_type;
	/// The texture coordiante type
	typedef TCTYPE texcoord_type;


protected:
	/// Vector of vertex pointer incident in the face
	FVTYPE *v[3];
	/// This are the flags of face, the default value is 0
	int  flags;		


/*#*******************	
*  Quality           *
**********************/

#ifdef __VCGLIB_FACE_Q
protected:
	float quality;
#endif
public:
	float & Q()
	{
#ifdef __VCGLIB_FACE_Q
		return quality;
#else
		assert(0);
		return *(float*)(&flags);
#endif
	}

const float & Q() const
	{
#ifdef __VCGLIB_FACE_Q
		return quality;
#else
		assert(0);
		return *(float*)(&flags);
#endif
	}

/*#*******************	
*  Texture           *
**********************/

#ifdef __VCGLIB_FACE_T
	TCTYPE t;
#endif
public:
	TCTYPE & T()
	{
#ifdef __VCGLIB_FACE_T
		return t;
#else
		assert(0);
		return *(TCTYPE*)(&flags);
#endif
	}

const TCTYPE & T() const
	{
#ifdef __VCGLIB_FACE_T
		return t;
#else
		assert(0);
		return *(TCTYPE*)(&flags);
#endif
	}

// Per Wedge Texture Coords
protected:
#ifdef __VCGLIB_FACE_WT
	TCTYPE wt[3];
#endif
public:
	TCTYPE & WT(const int i)
	{
#ifdef __VCGLIB_FACE_WT
		return wt[i];
#else
		assert(0);
		return *(TCTYPE*)(&flags);
#endif
	}

	const TCTYPE & WT(const int i) const
	{
#ifdef __VCGLIB_FACE_WT
		return wt[i];
#else
		assert(0);
		return *(TCTYPE*)(&flags);
#endif
	}


/*#*******************	
*  Colori            *
**********************/
protected:
#ifdef __VCGLIB_FACE_C
	ColorUB c;
#endif

public:
	ColorUB & C()
	{
#ifdef __VCGLIB_FACE_C
		return c;
#else
		assert(0);
		return *(ColorUB*)(&flags);
#endif
	}

	const ColorUB C() const
	{
#ifdef __VCGLIB_FACE_C
		return c;
#else
		return ColorUB(ColorUB::White);
#endif
	}

protected:
#ifdef __VCGLIB_FACE_WC
	ColorUB wc[3];
#endif
public:
	ColorUB & WC(const int i)
	{
#ifdef __VCGLIB_FACE_WC
		return wc[i];
#else
		assert(0);
		return *(ColorUB*)(&flags);
#endif
	}

const ColorUB WC(const int i) const
	{
#ifdef __VCGLIB_FACE_WC
		return wc[i];
#else
		assert(0);
		return ColorUB(ColorUB::White);
#endif
	}


/*#*******************	
*  Normals           *
**********************/
protected:
#ifdef __VCGLIB_FACE_N
	/// This vector indicates the normal of the face (defines if FACE_N is defined)
	vectorial_type n;
#endif

#ifdef __VCGLIB_FACE_WN
	/// This vector indicates per wedge normal 
	vectorial_type wn[3];
#endif

public:
	vectorial_type & WN(const int i)
	{
#ifdef __VCGLIB_FACE_WN
		return wn[i];
#else
		assert(0);
		return *(vectorial_type *)(&flags);
#endif
	}

const vectorial_type & WN(const int i) const
	{
#ifdef __VCGLIB_FACE_WN
		return wn[i];
#else
		return vectorial_type();
#endif
	}
protected:

/*#*******************	
*  Adjacency data    *
**********************/
#if (defined(__VCGLIB_FACE_A) && defined(__VCGLIB_FACE_S))
	#error Error: You cannot specify face-to-face and shared topology together
#endif

#if (defined(__VCGLIB_FACE_V) && defined(__VCGLIB_FACE_S))
	#error Error: You cannot specify vertex-face and shared topology together
#endif


#if defined(__VCGLIB_FACE_A)
	/// Vector of face pointer, it's used to indicate the adjacency relations (defines if FACE_A is defined)
	FACE_TYPE   *ff[3];				// Facce adiacenti
	/// Index of the face in the arrival face 
	char zf[4];									
#endif

#ifdef __VCGLIB_FACE_V
	///Vettore di puntatori a faccia, utilizzato per indicare le adiacenze vertice faccia
	FACE_TYPE *fv[3];
	char zv[3];
#endif

#ifdef __VCGLIB_FACE_S
	///Vettore di puntatori a faccia, utilizzato per indicare le adiacenze vertice faccia
	FACE_TYPE *fs[3];
	char zs[3];
#endif

#ifdef __VCGLIB_FACE_M
	/// Incremental mark (defines if FACE_I is defined)
	int imark;
#endif // Mark


public:
	enum {
		OBJ_TYPE_A  = 0x0001,
		OBJ_TYPE_N  = 0x0002,
		OBJ_TYPE_M  = 0x0004,
		OBJ_TYPE_V  = 0x0008,
		OBJ_TYPE_S  = 0x0010,
		OBJ_TYPE_E  = 0x0020,
		OBJ_TYPE_C  = 0x0040,
		OBJ_TYPE_WN = 0x0080,
		OBJ_TYPE_WC = 0x0100,
		OBJ_TYPE_WT = 0x0200,
		OBJ_TYPE_Q  = 0x0400,
	};

	enum {
		OBJ_TYPE = 
#ifdef __VCGLIB_FACE_A
		OBJ_TYPE_A |
#endif
#ifdef __VCGLIB_FACE_N
		OBJ_TYPE_N |
#endif
#ifdef __VCGLIB_FACE_M
		OBJ_TYPE_M |
#endif
#ifdef __VCGLIB_FACE_V
		OBJ_TYPE_V |
#endif
#ifdef __VCGLIB_FACE_S
		OBJ_TYPE_S |
#endif
#ifdef __VCGLIB_FACE_E
		OBJ_TYPE_E |
#endif
#ifdef __VCGLIB_FACE_C
		OBJ_TYPE_C |
#endif
#ifdef __VCGLIB_FACE_WN
		OBJ_TYPE_WN |
#endif
#ifdef __VCGLIB_FACE_WC
		OBJ_TYPE_WC |
#endif
#ifdef __VCGLIB_FACE_WT
		OBJ_TYPE_WT |
#endif
#ifdef __VCGLIB_FACE_Q
		OBJ_TYPE_Q  |
#endif
		0
	};

	enum {
		// This bit indicate that the face is deleted from the mesh
		DELETED     = 0x00000001,		// cancellato
		// This bit indicate that the face of the mesh is not readable
		NOTREAD     = 0x00000002,		// non leggibile (ma forse modificabile)
		// This bit indicate that the face is not modifiable
		NOTWRITE    = 0x00000004,		// non modificabile (ma forse leggibile) 
		// This bit indicate that the face is modified
		MODIFIED    = 0x00000008,		// modificato 
		// This bit can be used to mark the visited face
		VISITED     = 0x00000010,		// Visited  
		// This bit can be used to select 
		SELECTED    = 0x00000020,		// Selection flags
		// Border flags, it is assumed that BORDERi = BORDER0<<i 
		BORDER0     = 0x00000040,
		BORDER1     = 0x00000080,
		BORDER2     = 0x00000100,
		// Face Orientation Flags, used efficiently compute point face distance  
		NORMX		= 0x00000200,
		NORMY		= 0x00000400,
		NORMZ		= 0x00000800,
		// Complex flags, it is assumed that BORDERi = BORDER0<<i 
		COMPLEX0    = 0x00001000,
		COMPLEX1    = 0x00002000,
		COMPLEX2    = 0x00004000,
		// Crease flags,  it is assumed that FEATUREi = FEATURE0<<i 
		FEATURE0    = 0x00008000,
		FEATURE1    = 0x00010000,
		FEATURE2    = 0x00020000,
		// Flags per Marcatura degli halfedge
		MHEDGE0     = 0x00040000,
		MHEDGE1     = 0x00080000,
		MHEDGE2     = 0x00100000,
		// First user bit
		USER0       = 0x00040000
		};
	
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

	inline FACE_TYPE() {
/*
#ifdef _DEBUG 
		flags=0;
#ifdef __VCGLIB_FACE_A
		f[0]=f[1]=f[2]=0;
#endif

#endif
*/
	};

/*#*******************	
*  Bounding box *
**********************/

void GetBBox( bbox_type & bb )
{
	bb.Set( v[0]->P() );
	bb.Add( v[1]->P() );
	bb.Add( v[2]->P() );
}

	/// Return the number of vertices of the face
	int size() const {return 3;}

	/** Return the pointer to the j-th vertex of the face.
		@param j Index of the face vertex.
	 */
	inline FVTYPE * & V( const int j )
	{	
		assert( (flags & DELETED) == 0 );
		assert( (flags & NOTREAD) == 0 ); 
		assert( (flags & NOTWRITE) == 0 );
		assert(j>=0);
		assert(j<size());
		return v[j];
	}

	inline const FVTYPE * const & V( const int j ) const
	{
		assert( (flags & DELETED) == 0 );
		assert( (flags & NOTREAD) == 0 );
		assert(j>=0);
		assert(j<size());
		return v[j];
	}
	inline const FVTYPE * const & cV( const int j ) const
	{
		assert( (flags & DELETED) == 0 );
		assert( (flags & NOTREAD) == 0 );
		assert(j>=0);
		assert(j<size());
		return v[j];
	}

	// Shortcut per accedere ai punti delle facce
	inline vectorial_type & P( const int j )
	{	
		assert( (flags & DELETED) == 0 );
		assert( (flags & NOTREAD) == 0 ); 
		assert( (flags & NOTWRITE) == 0 );
		assert(j>=0);
		assert(j<size());
		return v[j]->P();
	}

	inline const vectorial_type & P( const int j ) const
	{
		assert( (flags & DELETED) == 0 );
		assert( (flags & NOTREAD) == 0 );
		assert(j>=0);
		assert(j<size());
		return v[j]->cP();
	}
	inline const vectorial_type & cP( const int j ) const
	{
		assert( (flags & DELETED) == 0 );
		assert( (flags & NOTREAD) == 0 );
		assert(j>=0);
		assert(j<size());
		return v[j]->cP();
	}

	/** Return the pointer to the ((j+1)%3)-th vertex of the face.
		@param j Index of the face vertex.
	 */
	inline FVTYPE * & V0( const int j ) { return V(j);}
	inline FVTYPE * & V1( const int j ) { return V((j+1)%3);}
	inline FVTYPE * & V2( const int j ) { return V((j+2)%3);}
	inline const FVTYPE * const &  V0( const int j ) const { return V(j);}
	inline const FVTYPE * const &  V1( const int j ) const { return V((j+1)%3);}
	inline const FVTYPE * const &  V2( const int j ) const { return V((j+2)%3);}
	inline const FVTYPE * const & cV0( const int j ) const { return cV(j);}
	inline const FVTYPE * const & cV1( const int j ) const { return cV((j+1)%3);}
	inline const FVTYPE * const & cV2( const int j ) const { return cV((j+2)%3);}

	// Shortcut per accedere ai punti delle facce
	inline vectorial_type & P0( const int j ) { return V(j)->P();}
	inline vectorial_type & P1( const int j ) { return V((j+1)%3)->P();}
	inline vectorial_type & P2( const int j ) { return V((j+2)%3)->P();}
	inline const vectorial_type &  P0( const int j ) const { return V(j)->P();}
	inline const vectorial_type &  P1( const int j ) const { return V((j+1)%3)->P();}
	inline const vectorial_type &  P2( const int j ) const { return V((j+2)%3)->P();}
	inline const vectorial_type & cP0( const int j ) const { return cV(j)->P();}
	inline const vectorial_type & cP1( const int j ) const { return cV((j+1)%3)->P();}
	inline const vectorial_type & cP2( const int j ) const { return cV((j+2)%3)->P();}

	inline FVTYPE * & Supervisor_V( const int j )
	{	
		assert(j>=0);
		assert(j<size());
		return v[j];
	}

	inline const FVTYPE * const & Supervisor_V( const int j ) const
	{
		assert(j>=0);
		assert(j<size());
		return v[j];
	}

	/// Return the flags.
	inline int & Flags ()
	{
		assert( (flags & DELETED) == 0 );
		assert( (flags & NOTREAD) == 0 );
		return flags;
	}

	inline const int & Flags () const
	{
		assert( (flags & DELETED) == 0 );
		assert( (flags & NOTREAD) == 0 );
		return flags;
	}
	/// Ritorna il flags senza effettuare alcun controllo sui relativi bit
	inline int & Supervisor_Flags()
	{
		return flags;
	}

	inline const int Supervisor_Flags() const
	{
		return flags;
	}

	/// Return the reference of the normal to the face (if __VCGLIB_FACE_N is defined).
	inline vectorial_type & N()
	{
#ifdef __VCGLIB_FACE_N
	return n;
#else
	assert(0);
	return *(vectorial_type *)0;
#endif
	}
		/// Return the reference of the normal to the face (if __VCGLIB_FACE_N is defined).
	inline const vectorial_type & N() const
	{
#ifdef __VCGLIB_FACE_N
		return n;
#else
		return vcg::Normal(V(0)->P(), V(1)->P(), V(2)->P());
#endif
	}
	/// Return the reference of the normal to the face (if __VCGLIB_FACE_N is defined).
	inline const vectorial_type cN() const
	{
		return Normal();
	}

	/** Return the pointer to the j-th adjacent face.
	    @param j Index of the edge.
	 */
	inline FACE_TYPE * & F( const int j )
	{
		assert( (flags & DELETED) == 0 );
		assert( (flags & NOTREAD) == 0 );
		assert( (flags & NOTWRITE) == 0 );
		assert(j>=0);
	    assert(j<size());
#if defined(__VCGLIB_FACE_A)
		  return ff[j];
#elif defined(__VCGLIB_FACE_S)
			return fs[j];
#else 
		assert(0);
        static FACE_TYPE *dum=0;
		return dum;
#endif
	}

	inline const FACE_TYPE * const & F( const int j ) const
	{
		assert( (flags & DELETED) == 0 );
		assert( (flags & NOTREAD) == 0 );
		assert(j>=0);
	    assert(j<size());
#if defined(__VCGLIB_FACE_A)
		  return ff[j];
#elif defined(__VCGLIB_FACE_S)
			return fs[j];
#else
		  assert(0);
		  return (FACE_TYPE *)this;
#endif
	}
	inline FACE_TYPE * & F1( const int j ) { return F((j+1)%3);}
	inline FACE_TYPE * & F2( const int j ) { return F((j+2)%3);}
	inline const FACE_TYPE * const&  F1( const int j ) const { return F((j+1)%3);}
	inline const FACE_TYPE * const&  F2( const int j ) const { return F((j+2)%3);}


/** Return the pointer to the j-th adjacent face.
	    @param j Index of the edge.
	 */
	inline FACE_TYPE * & Supervisor_F( const int j )
	{
		assert(j>=0);
	  assert(j<size());
#if defined(__VCGLIB_FACE_A)
		  return ff[j];
#elif defined(__VCGLIB_FACE_S)
			return fs[j];
#else 
		assert(0); // if you stop here you are probably trying to use FF topology in a face without it
		return *((FACE_TYPE **)(flags));
#endif
	}

	inline const FACE_TYPE * const & Supervisor_F( const int j ) const
	{
		assert(j>=0);
	  assert(j<size());
#if defined(__VCGLIB_FACE_A)
		  return ff[j];
#elif defined(__VCGLIB_FACE_S)
			return fs[j];
#else
		assert(0); // if you stop here you are probably trying to use FF topology in a face without it
		return *((FACE_TYPE **)(flags));
#endif
	}
	

	inline FACE_TYPE * & Fv( const int j )
	{
		assert( (flags & DELETED) == 0 );
		assert( (flags & NOTREAD) == 0 );
		assert( (flags & NOTWRITE) == 0 );
		assert(j>=0);
		assert(j<size());
#ifdef __VCGLIB_FACE_V
		return fv[j];
#elif defined(__VCGLIB_FACE_S)
		return fs[j];
#else
		assert(0); // you are probably trying to use VF topology in a vertex without it
		return *((FACE_TYPE **)(flags));
#endif
	}

	inline const FACE_TYPE * const & Fv( const int j ) const
	{
		assert( (flags & DELETED) == 0 );
		assert( (flags & NOTREAD) == 0 );
		assert(j>=0);
		assert(j<size());
#ifdef __VCGLIB_FACE_V
		return fv[j];
#elif defined(__VCGLIB_FACE_S)
		return fs[j];
#else
		assert(0);
		return (FACE_TYPE *)this;
#endif
	}


	/** Return the index that the face have in the j-th adjacent face.
	    @param j Index of the edge.
	 */
	inline char & Z( const int j )
	{
		assert( (flags & DELETED) == 0 );
		assert( (flags & NOTREAD) == 0 );
		assert( (flags & NOTWRITE) == 0 );
		assert(j>=0);
		assert(j<size());
#if defined(__VCGLIB_FACE_A) 
		return zf[j];
#elif defined(__VCGLIB_FACE_S) 
		return zs[j];
#else
		assert(0);
		return *(char *)&flags; // tanto per farlo compilare...
#endif
	}

	inline const char & Z( const int j ) const
	{
		assert( (flags & DELETED) == 0 );
		assert( (flags & NOTREAD) == 0 );
		assert(j>=0);
		assert(j<size());
#if defined(__VCGLIB_FACE_A) 
		return zf[j];
#elif defined(__VCGLIB_FACE_S) 
		return zs[j];
#else
		assert(0);
		return *(char *)&flags;
#endif
	}

		/** Return the index that the face have in the j-th adjacent face.
	    @param j Index of the edge.
	 */
	inline char & Supervisor_Z( const int j )
	{
		assert(j>=0);
		assert(j<size());
#if defined(__VCGLIB_FACE_A) 
		return zf[j];
#elif defined(__VCGLIB_FACE_S) 
		return zs[j];
#else
		assert(0);
		return *(char *)&flags;
#endif
	}

	inline const char & Supervisor_Z( const int j ) const
	{
		assert(j>=0);
		assert(j<size());
#if defined(__VCGLIB_FACE_A) 
		return zf[j];
#elif defined(__VCGLIB_FACE_S) 
		return zs[j];
#else
		assert(0);
		return *(char *)&flags;
#endif
	}


	inline char & Zv( const int j )
	{
		assert( (flags & DELETED) == 0 );
		assert( (flags & NOTREAD) == 0 );
		assert( (flags & NOTWRITE) == 0 );
		assert(j>=0);
		assert(j<size());
#ifdef __VCGLIB_FACE_V
		return zv[j];
#elif defined(__VCGLIB_FACE_S)
		return zs[j];
#else
		assert(0);
		return *(char *)&flags;
#endif
	}

	inline const char & Zv( const int j ) const
	{
		assert( (flags & DELETED) == 0 );
		assert( (flags & NOTREAD) == 0 );
		assert(j>=0);
		assert(j<size());
#ifdef __VCGLIB_FACE_V
		return zv[j];
#elif defined(__VCGLIB_FACE_S)
		return zs[j];
#else
		assert(0);
		return *(char *)&flags;
#endif
	}

//#endif

	/* non posso usare vettori di facce
	inline FVTYPE * & operator [] ( const int i ){
			assert(i>=0 && i<3);
			return v[i];
	}
	inline const FVTYPE * & operator [] ( const int i ) const {
			assert(i>=0 && i<3);
			return v[i];
	}
	*/

	/// operator to compare two faces
	inline bool operator == ( const FACE_TYPE & f ) const {
		for(int i=0; i<size(); ++i)
			if( (V(i) != f.V(0)) && (V(i) != f.V(1)) && (V(i) != f.V(2)) )
				return false;
		return true;
	}


	/*inline bool operator < ( const FACE_TYPE & f ) const {
	}*/



	/// This function checks if the face is deleted 
	bool IsD() const {return (flags & DELETED) != 0;}
	/// This function mark the face as deleted
	void SetD()		{flags |=DELETED;}
	/// This function mark the face as not deleted
	void ClearD()	{flags &= (~DELETED);}
	/// This function checks if the face is deleted 
	bool IsDeleted() const {return IsD();}
	
	/// This function checks if the face is readable 
	bool IsR() const {return (flags & NOTREAD) == 0;}
	/// This function marks the face as readable
	void SetR()		{flags &= (~NOTREAD);}
	/// This function marks the face as not readable
	void ClearR() {flags |=NOTREAD;}
	/// This function checks if the face is readable
	bool IsReadable() const {return IsR();}
	
	/// This function checks if the face is readable 
	bool IsW() const {return (flags & NOTWRITE)== 0;}
	/// This function marks the vertex as not writable
	void SetW() {flags &=(~NOTWRITE);}
	/// This function marks the face as not writable
	void ClearW() {flags |=NOTWRITE;}

	/// This function checks if the face is modifiable
	bool IsWritable() const {return IsW();}
	/// This funcion checks whether the face is both readable and modifiable
	bool IsRW() const {return (flags & (NOTREAD | NOTWRITE)) == 0;}

	/// This function checks if the face is Modified
	bool IsM() const {return (flags & MODIFIED)!= 0;}
	/// This function marks the face as modified. It's necessary to mark all modified faces to have a consistent mesh
	void SetM()		{flags |=MODIFIED;}
	/// This function marks the face as not visited. This flag, initially, is setted to random value, therefore, to the beginnig of every function it is necessary to clean up the flag
	void ClearM() {flags &= (~MODIFIED);}

	/// This function checks if the face is marked as visited
	bool IsV() const {return (flags & VISITED)!= 0;}
	/// This funcion marks the face as visited
	void SetV()		{flags |=VISITED;}
	/// This function marks the face as not visited. This flag, initially, is setted to random value, therefore, to the beginnig of every function it is necessary to clean up the flag
	void ClearV() {flags &= (~VISITED);}
	
	/// This function checks if the face is selected
	bool IsS() const {return (flags & SELECTED) != 0;}
	/// This function select the face
	void SetS()		{flags |=SELECTED;}
	/// This funcion execute the inverse operation of SetS()
	void ClearS()	{flags &= (~SELECTED);}

	/// This function checks if the face is selected
	bool IsB(int i) const {return (flags & (BORDER0<<i)) != 0;}
	/// This function select the face
	void SetB(int i)		{flags |=(BORDER0<<i);}
	/// This funcion execute the inverse operation of SetS()
	void ClearB(int i)	{flags &= (~(BORDER0<<i));}

	/// This function checks if the face is Complex  on side i
	bool IsCF(int i) const {return (flags & (COMPLEX0<<i)) != 0;}
	/// This function select the face
	void SetCF(int i)		{flags |=(COMPLEX0<<i);}
	/// This funcion execute the inverse operation of SetS()
	void ClearCF(int i)	{flags &= (~(COMPLEX0<<i));}

	/// This function checks if the face is Crease  on side i
	bool IsFF(int i) const {return (flags & (FEATURE0<<i)) != 0;}
	/// This function select the face flag
	void SetFF(int i)		{flags |=(FEATURE0<<i);}
	/// This funcion execute the inverse operation of Set()
	void ClearFF(int i)	{flags &= (~(FEATURE0<<i));}

	/// This function checks if the half edge i is marked
	bool IsHEdgeM(int i) const {return (flags & (MHEDGE0<<i)) != 0;}
	/// This function set the half edge mark
	void SetHEdgeM(int i)		{flags |=(MHEDGE0<<i);}
	/// This funcion execute the inverse operation of the previuse one
	void ClearHEdgeM(int i)	{flags &= (~(MHEDGE0<<i));}

	/// This function checks if the given user bit is true
	bool IsUserBit(int userBit){return (flags & userBit) != 0;}
	/// This function set  the given user bit 
	void SetUserBit(int userBit){flags |=userBit;}
	/// This function clear the given user bit 
	void ClearUserBit(int userBit){flags &= (~userBit);}

/*#*******************	
*  Normals           *
**********************/
/// Calculate the normal to the face, the value is store in the field n of the face
void ComputeNormal() 
{
#ifdef __VCGLIB_FACE_N
	n = vcg::Normal(V(0)->cP(), V(1)->cP(), V(2)->cP());
#else
	assert(0);
#endif
}
void ComputeNormalizedNormal() 
{
#ifdef __VCGLIB_FACE_N
	n = vcg::NormalizedNormal(V(0)->cP(), V(1)->cP(), V(2)->cP());
#else
	assert(0);
#endif
}

/// Return the value of the face normal; warning: if __VCGLIB_FACE_N is not defined the value is computed each time
vectorial_type Normal() const
{
#ifdef __VCGLIB_FACE_N
	return n;
#else
	return vcg::Normal(V(0)->P(), V(1)->P(), V(2)->P());
#endif
}

/** Calcola i coefficienti della combinazione convessa.
	@param bq Punto appartenente alla faccia
	@param a Valore di ritorno per il vertice V(0)
	@param b Valore di ritorno per il vertice V(1)
	@param c Valore di ritorno per il vertice V(2)
	@return true se bq appartiene alla faccia, false altrimenti
*/
bool InterpolationParameters(const vectorial_type & bq, scalar_type &a, scalar_type &b, scalar_type &c ) const
{
	/********** VECCHIA VERSIONE *********************
	//Calcolo degli assi dominanti
	int axis; // asse piu' perpendicolare alla faccia
	scalar_type max;
	vectorial_type fnorm = Normal();
	max = Abs(fnorm[0]);
	axis = 0;
	if( max < Abs(fnorm[1]) )
	{
		max = Abs(fnorm[1]);
		axis = 1;
	}
	if( max < Abs(fnorm[2]) )
	{
		max = Abs(fnorm[2]);
		axis = 2;
	}

	scalar_type x[3];
	scalar_type C[3][3+1];
	
	C[0][0]	= cV(0)->P()[(axis+1)%3]; 
	C[0][1]	= cV(1)->P()[(axis+1)%3];
	C[0][2]	= cV(2)->P()[(axis+1)%3];
	C[0][3]	=       bq  [(axis+1)%3];
	
	C[1][0]	= cV(0)->P()[(axis+2)%3]; 
	C[1][1]	= cV(1)->P()[(axis+2)%3];
	C[1][2]	= cV(2)->P()[(axis+2)%3];
	C[1][3]	=       bq  [(axis+2)%3];
	
	C[2][0]	= 1; 
	C[2][1]	= 1;
	C[2][2]	= 1;
	C[2][3]	= 1;

	if(Gauss33(x,C))
		{
			a=x[0];
			b=x[1];
			c=x[2];
			return true;
		}
	else
		{
			a=b=c=1.0/3.0;
			return false;
		}
********** FINE VECCHIA VERSIONE ***************/
	
	const scalar_type EPSILON = scalar_type(0.000001);


#define x1 (cV(0)->P().x())
#define y1 (cV(0)->P().y())
#define z1 (cV(0)->P().z())
#define x2 (cV(1)->P().x())
#define y2 (cV(1)->P().y())
#define z2 (cV(1)->P().z())
#define x3 (cV(2)->P().x())
#define y3 (cV(2)->P().y())
#define z3 (cV(2)->P().z())
#define px (bq.x())
#define py (bq.y())
#define pz (bq.z())

     scalar_type t1  = px*y2;
     scalar_type t2  = px*y3;
     scalar_type t3  = py*x2;
     scalar_type t4  = py*x3;
     scalar_type t5  = x2*y3;
     scalar_type t6  = x3*y2;
     scalar_type t8  = x1*y2;
     scalar_type t9  = x1*y3;
     scalar_type t10 = y1*x2;
     scalar_type t11 = y1*x3;
     scalar_type t13 = t8-t9-t10+t11+t5-t6;
     if(fabs(t13)>=EPSILON)
	 {
         scalar_type t15 = px*y1;
         scalar_type t16 = py*x1;
         a =  (t1 -t2-t3 +t4+t5-t6 )/t13;
         b = -(t15-t2-t16+t4+t9-t11)/t13;
         c =  (t15-t1-t16+t3+t8-t10)/t13;
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
		scalar_type t15 = px*z1;
		scalar_type t16 = pz*x1;
		a =  (t1 -t2-t3 +t4+t5-t6 )/t13;
		b = -(t15-t2-t16+t4+t9-t11)/t13;
		c =  (t15-t1-t16+t3+t8-t10)/t13;
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
         scalar_type t15 = pz*y1;
         scalar_type t16 = py*z1;
         a =  (t1 -t2-t3 +t4+t5-t6 )/t13;
         b = -(t15-t2-t16+t4+t9-t11)/t13;
         c =  (t15-t1-t16+t3+t8-t10)/t13;
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


/*#*******************	
*  Adjacency Members *
**********************/


/** Return a boolean that indicate if the face is complex.
    @param j Index of the edge
	@return true se la faccia e' manifold, false altrimenti
*/
inline bool IsManifold( const int j ) const
{
#if (defined(__VCGLIB_FACE_A) || defined(__VCGLIB_FACE_S))
	return ( F(j)==this || this == F(j)->F(Z(j)) );
#endif
	return true;
	assert(0);
}

/** Return a boolean that indicate if the j-th edge of the face is a border.
	@param j Index of the edge
	@return true if j is an edge of border, false otherwise
*/
inline bool IsBorder( const int j ) const
{
#if (defined(__VCGLIB_FACE_A) || defined(__VCGLIB_FACE_S))
	return F(j)==this;
#endif
	return true;
	assert(0);
}


/// This function counts the boreders of the face
inline int BorderCount() const
{
#if (defined(__VCGLIB_FACE_A) || defined(__VCGLIB_FACE_S))
	int t = 0;
	if( IsBorder(0) ) ++t;
	if( IsBorder(1) ) ++t;
	if( IsBorder(2) ) ++t;
	return t;
#endif
	assert(0);
	return 3;
}


/// This function counts the number of incident faces in a complex edge
inline int ComplexSize(const int e) const
{
#if (defined(__VCGLIB_FACE_A) || defined(__VCGLIB_FACE_S))
	int cnt=0;
	FACE_TYPE *fi=(FACE_TYPE *)this;
	int nzi,zi=e;
	do
	{
		nzi=fi->Z(zi);
		fi=fi->F(zi);
		zi=nzi;
		++cnt;
	}
	while(fi!=this);
	return cnt;
#endif
	assert(0);
	return 2;
}

/*Funzione di detach che scollega una faccia da un ciclo 
(eventualmente costituito da due soli elementi) incidente su un edge*/
/** This function detach the face from the adjacent face via the edge e. It's possible to use it also in non-two manifold situation.
		The function cannot be applicated if the adjacencies among faces aren't define.
		@param e Index of the edge
*/
void Detach(const int e)
{
	typedef FEdgePosB< FACE_TYPE > ETYPE;

	assert(!IsBorder(e));
	ETYPE EPB(this,e);  // la faccia dall'altra parte
	EPB.NextF();
	int cnt=0;
	while ( EPB.f->F(EPB.z) != this)
	{ 
		assert(!IsManifold(e));   // Si entra in questo loop solo se siamo in una situazione non manifold.
		assert(!EPB.f->IsBorder(EPB.z));
		EPB.NextF();
		cnt++;
	}
	assert(EPB.f->F(EPB.z)==this);

	EPB.f->F(EPB.z) = F(e);
	EPB.f->Z(EPB.z) = Z(e);
	
	F(e) = this;
	Z(e) = e;

	EPB.f->SetM();
	this->SetM();
}

void OldDetach(const int e)
{
	typedef EdgePosB< FACE_TYPE > ETYPE;

	assert(!IsBorder(e));
	ETYPE EPB(this,e);
	ETYPE TEPB(0,-1);
	EPB.NextF();
	while ( EPB.f != this)
	{
		TEPB = EPB;
		assert(!EPB.f->IsBorder(EPB.z));
		EPB.NextF();
	}
	assert(TEPB.f->F(TEPB.z)==this);
	TEPB.f->F(TEPB.z) = F(e);
	TEPB.f->Z(TEPB.z) = Z(e);
	F(e) = this;
	Z(e) = e;
	TEPB.f->SetM();
	this->SetM();
}


/** This function attach the face (via the edge z1) to another face (via the edge z2). It's possible to use it also in non-two manifold situation.
		The function cannot be applicated if the adjacencies among faces aren't define.
		@param z1 Index of the edge
		@param f2 Pointer to the face
		@param z2 The edge of the face f2 
*/
void Attach(int z1, face_base *&f2, int z2)
{
	typedef FEdgePosB< FACE_TYPE > ETYPE;
	ETYPE EPB(f2,z2);
	ETYPE TEPB;
	TEPB = EPB;
	EPB.NextF();
	while( EPB.f != f2)  //Alla fine del ciclo TEPB contiene la faccia che precede f2
	{
		TEPB = EPB;
		EPB.NextF();
	}
	//Salvo i dati di f1 prima di sovrascrivere
	face_base *f1prec = this->F(z1);  
	int z1prec = this->Z(z1);
	//Aggiorno f1
	this->F(z1) = TEPB.f->F(TEPB.z);  
	this->Z(z1) = TEPB.f->Z(TEPB.z);
	//Aggiorno la faccia che precede f2
	TEPB.f->F(TEPB.z) = f1prec;
	TEPB.f->Z(TEPB.z) = z1prec;
}


void AssertAdj()
{
	assert(F(0)->F(Z(0))==this);
	assert(F(1)->F(Z(1))==this);
	assert(F(2)->F(Z(2))==this);

	assert(F(0)->Z(Z(0))==0);
	assert(F(1)->Z(Z(1))==1);
	assert(F(2)->Z(Z(2))==2); 
}
//#endif // Adjacency

/*#**************
*  Mark Members *
*****************/
/// Return the incremental mark of the face
#ifdef __VCGLIB_FACE_M
	inline int & IMark()
	{
		assert( (flags & DELETED) == 0 );
		assert( (flags & NOTREAD) == 0 );
		assert( (flags & NOTWRITE) == 0 );
		return imark;
	}

	inline const int & IMark() const
	{
		assert( (flags & DELETED) == 0 );
		assert( (flags & NOTREAD) == 0 );
		return imark;
	}
#endif // Mark

	/// Initialize the imark system of the face
	inline void InitIMark()
	{
#ifdef __VCGLIB_FACE_M
		imark = 0;
#endif
	}

/// Return the DOUBLE of the area of the face
FCTYPE Area() const
{
	return Norm( (V(1)->P() - V(0)->P()) ^ (V(2)->P() - V(0)->P()) );
}

vectorial_type Barycenter() const
{
	return (V(0)->P()+V(1)->P()+V(2)->P())/FCTYPE(3.0);
}

FCTYPE Perimeter() const
{
	return Distance(V(0)->P(),V(1)->P())+
		     Distance(V(1)->P(),V(2)->P())+
				 Distance(V(2)->P(),V(0)->P());
}

/// Return the quality of the face, the return value is in [0,sqrt(3)/2] = [0 - 0.866.. ]
FCTYPE QualityFace( ) const
{
	
	return Quality(V(0)->P(), V(1)->P(), V(2)->P());
	/*
	vectorial_type d10 = V(1)->P() - V(0)->P();
	vectorial_type d20 = V(2)->P() - V(0)->P();
	vectorial_type d12 = V(1)->P() - V(2)->P();

	vectorial_type x = d10^d20;

	FCTYPE a = Norm( x );		// doppio dell' Area
	FCTYPE b;
	
	b = Norm2( d10 );
	FCTYPE t = b; 
	t = Norm2( d20 ); if( b<t ) b = t;
	t = Norm2( d12 ); if( b<t ) b = t;

	assert(b!=0.0);

	return a/b;*/

}

// Funzione di supporto
inline void Nexts( face_base *&f,int &z )
{
    int t;
    t = z;
    z = (*f).Z(z);
    f = (*f).F(t);
}

/** This function change the orientation of the face. Inverting the index of two vertex 
@param z Index of the edge
*/
void Swap ( const int z )
{

  int i;
  face_base *tmp, *prec;
  int t, precz;

  swap ( V((z  )%3),V((z+1)%3));

  if( OBJ_TYPE & (OBJ_TYPE_A|OBJ_TYPE_S ) )
  {
	swap ( F((z+1)%3),F((z+2)%3));
	swap ( Z((z+1)%3),Z((z+2)%3));

	for(i = 1; i < 3; i++)
	{

      tmp = this;
      t = (z+i)%3;
      do {
					prec = tmp;
					precz = t;
					Nexts(tmp,t);
      }
      while (tmp != this);
  
      (*prec).Z(precz) = (z+i)%3;
    }
  }
}

/*void Swap( const int e )
{
	swap( V(e), V((e+1)%3) );

#ifdef __VCGLIB_FACE_A
	
	swap( F((e+1)%3), F((e+2)%3) );
	swap( Z((e+1)%3), Z((e+2)%3) );

	F((e+1)%3)->Z(Z((e+1)%3)) = (e+1)%3;
	F((e+2)%3)->Z(Z((e+2)%3)) = (e+2)%3;
#endif
}*/

// Stacca la faccia corrente dalla catena di facce incidenti sul vertice z, 
// NOTA funziona SOLO per la topologia VF!!!
// usata nelle classi di collapse
void VFDetach(int z)
{
	
	if(V(z)->Fp()==this )
	{
		int fz = V(z)->Zp();
		V(z)->Fp() = (face_from_vert_type *) F(fz);
		V(z)->Zp() = Z(fz);
	}
	else
	{
			VEdgePosB<FACE_TYPE> x,y;

		x.f = V(z)->Fp();
		x.z = V(z)->Zp();

		for(;;)
		{
			y = x;
			x.NextF();
			assert(x.f!=0);
			if(x.f==this)
			{
				y.f->F(y.z) = F(z);
				y.f->Z(y.z) = Z(z);
				break;
			}
		}
	}
}

	// Sezione dist e ray
#ifdef __VCGLIB_FACE_E
	vectorial_type edge[3];
	Plane3<scalar_type> plane;

	void ComputeE()
	{	
			// Primo calcolo degli edges
		edge[0] = V(1)->P(); edge[0] -= V(0)->P();
		edge[1] = V(2)->P(); edge[1] -= V(1)->P();
		edge[2] = V(0)->P(); edge[2] -= V(2)->P();
			// Calcolo di plane
		plane.n = edge[0]^edge[1];
		plane.d = plane.n * V(0)->P();
		plane.Normalize();
			// Calcolo migliore proiezione
		scalar_type nx = Abs(plane.n[0]);
		scalar_type ny = Abs(plane.n[1]);
		scalar_type nz = Abs(plane.n[2]);
		scalar_type d;
		if(nx>ny && nx>nz) { flags |= NORMX; d = 1/plane.n[0]; }
		else if(ny>nz)     { flags |= NORMY; d = 1/plane.n[1]; }
		else               { flags |= NORMZ; d = 1/plane.n[2]; }

			// Scalatura spigoli
		edge[0] *= d;
		edge[1] *= d;
		edge[2] *= d;
	}

/*
   Point face distance
   trova il punto <p> sulla faccia piu' vicino a <q>, con possibilità di 
   rejection veloce su se la distanza trovata è maggiore di <rejdist>

 Commenti del 12/11/02
 Funziona solo se la faccia e di quelle di tipo E (con edge e piano per faccia gia' calcolati)
 algoritmo:
	1) si calcola la proiezione <p> di q sul piano della faccia
	2) se la distanza punto piano e' > rejdist ritorna
	3) si lavora sul piano migliore e si cerca di capire se il punto sta dentro il triangolo:
	   a) prodotto vettore tra edge triangolo (v[i+1]-v[i]) e (p-v[i])
		 b) se il risultato e' negativo (gira in senso orario) allora il punto
		    sta fuori da quella parte e si fa la distanza punto segmento.
     c) se il risultato sempre positivo allora sta dentro il triangolo
	4) e si restituisce la distanza punto /piano gia` calcolata 

	Note sulla robustezza:
	il calcolo del prodotto vettore e` la cosa piu` delicata:
	possibili fallimenti quando a^b ~= 0
	1) doveva essere <= 0 e viene positivo (q era fuori o sulla linea dell'edge)
	   allora capita che si faccia la distanza punto piano anziche` la distanza punto seg
  2) doveva essere > 0 e viene <=0 (q era dentro il triangolo)

*/
	bool Dist( const vectorial_type & q, scalar_type & dist, vectorial_type & p )
	{
		//const scalar_type EPSILON = scalar_type( 0.000001);
		const scalar_type EPSILON = 0.00000001;
		scalar_type b,b0,b1,b2;
			// Calcolo distanza punto piano
		scalar_type d = Distance( plane, q );
		if( d>dist || d<-dist )			// Risultato peggiore: niente di fatto
			return false;

			// Calcolo del punto sul piano
		// NOTA: aggiunto un '-d' in fondo Paolo C.
		vectorial_type t = plane.n;
		t[0] *= -d;
		t[1] *= -d;
		t[2] *= -d;
		p = q; p += t;
		
	#define PP(i)	(v[i]->P())
	#define E(i)    (edge[i])
    
		switch( flags & (NORMX|NORMY|NORMZ) )
		{
		case NORMX:
			b0 = E(1)[1]*(p[2] - PP(1)[2]) - E(1)[2]*(p[1] - PP(1)[1]);
			if(b0<=0)
			{
				b0 = PSDist(q,V(1)->P(),V(2)->P(),p);
				if(dist>b0) { dist = b0; return true; }
				else return false;
			}
			b1 = E(2)[1]*(p[2] - PP(2)[2]) - E(2)[2]*(p[1] - PP(2)[1]);
			if(b1<=0)
			{
				b1 = PSDist(q,V(2)->P(),V(0)->P(),p);
				if(dist>b1) { dist = b1; return true; }
				else return false;
			}
			b2 = E(0)[1]*(p[2] - PP(0)[2]) - E(0)[2]*(p[1] - PP(0)[1]);
			if(b2<=0)
			{
				b2 = PSDist(q,V(0)->P(),V(1)->P(),p);
				if(dist>b2) { dist = b2; return true; }
				else return false;
			}
			// sono tutti e tre > 0 quindi dovrebbe essere dentro;
			// per sicurezza se il piu' piccolo dei tre e' < epsilon (scalato rispetto all'area della faccia
			// per renderlo dimension independent.) allora si usa ancora la distanza punto 
			// segmento che e' piu robusta della punto piano, e si fa dalla parte a cui siamo piu' 
			// vicini (come prodotto vettore)
			// Nota: si potrebbe rendere un pochino piu' veloce sostituendo Area()
			// con il prodotto vettore dei due edge in 2d lungo il piano migliore.
			if( (b=min(b0,min(b1,b2))) < EPSILON*Area()) 
      {
				scalar_type bt;
				if(b==b0) 	    bt = PSDist(q,V(1)->P(),V(2)->P(),p);
				else if(b==b1) 	bt = PSDist(q,V(2)->P(),V(0)->P(),p);
				else if(b==b2) 	bt = PSDist(q,V(0)->P(),V(1)->P(),p);
				//printf("Warning area:%g %g %g %g thr:%g bt:%g\n",Area(), b0,b1,b2,EPSILON*Area(),bt);
				if(dist>bt) { dist = bt; return true; }
				else return false;
			}
			break;

		  case NORMY:
			b0 = E(1)[2]*(p[0] - PP(1)[0]) - E(1)[0]*(p[2] - PP(1)[2]);
			if(b0<=0)
			{
				b0 = PSDist(q,V(1)->P(),V(2)->P(),p);
				if(dist>b0) { dist = b0; return true; }
				else return false;
			}
			b1 = E(2)[2]*(p[0] - PP(2)[0]) - E(2)[0]*(p[2] - PP(2)[2]);
			if(b1<=0)
			{
				b1 = PSDist(q,V(2)->P(),V(0)->P(),p);
				if(dist>b1) { dist = b1; return true; }
				else return false;
			}
			b2 = E(0)[2]*(p[0] - PP(0)[0]) - E(0)[0]*(p[2] - PP(0)[2]);
			if(b2<=0)
			{
				b2 = PSDist(q,V(0)->P(),V(1)->P(),p);
				if(dist>b2) { dist = b2; return true; }
				else return false;
			}
			if( (b=min(b0,min(b1,b2))) < EPSILON*Area()) 
      {
				scalar_type bt;
				if(b==b0) 	    bt = PSDist(q,V(1)->P(),V(2)->P(),p);
				else if(b==b1) 	bt = PSDist(q,V(2)->P(),V(0)->P(),p);
				else if(b==b2) 	bt = PSDist(q,V(0)->P(),V(1)->P(),p);
				//printf("Warning area:%g %g %g %g thr:%g bt:%g\n",Area(), b0,b1,b2,EPSILON*Area(),bt);
				if(dist>bt) { dist = bt; return true; }
				else return false;
			}
			break;

		  case NORMZ:
			b0 = E(1)[0]*(p[1] - PP(1)[1]) - E(1)[1]*(p[0] - PP(1)[0]);
			if(b0<=0)
			{
				b0 = PSDist(q,V(1)->P(),V(2)->P(),p);
				if(dist>b0) { dist = b0; return true; }
				else return false;
			}
			b1 = E(2)[0]*(p[1] - PP(2)[1]) - E(2)[1]*(p[0] - PP(2)[0]);
			if(b1<=0)
			{
				b1 = PSDist(q,V(2)->P(),V(0)->P(),p);
				if(dist>b1) { dist = b1; return true; }
				else return false;
			}
			b2 = E(0)[0]*(p[1] - PP(0)[1]) - E(0)[1]*(p[0] - PP(0)[0]);
			if(b2<=0)
			{
				b2 = PSDist(q,V(0)->P(),V(1)->P(),p);
				if(dist>b2) { dist = b2; return true; }
				else return false;
			}
			if( (b=min(b0,min(b1,b2))) < EPSILON*Area()) 
      {
				scalar_type bt;
				if(b==b0) 	    bt = PSDist(q,V(1)->P(),V(2)->P(),p);
				else if(b==b1) 	bt = PSDist(q,V(2)->P(),V(0)->P(),p);
				else if(b==b2) 	bt = PSDist(q,V(0)->P(),V(1)->P(),p);
				//printf("Warning area:%g %g %g %g thr:%g bt:%g\n",Area(), b0,b1,b2,EPSILON*Area(),bt);
				
				if(dist>bt) { dist = bt; return true; }
				else return false;
			}
			break;
		}

	#undef E
	#undef PP

		dist = scalar_type(fabs(d));
		//dist = Distance(p,q);
		return true;
	}


	//Intersect ray with triangle.  This is an optimized version of the
    //intersection routine from Snyder and Barr's '87 SIGGRAPH paper.
    //Si distingue fra distanza negativa e positiva (N e P)
    //Returns 1 or 0
	static inline bool CheckVal( const scalar_type b ){
		return b<-1e-20 || b>1.0+1e-20 ;
	}

	bool Intersect( Line3<scalar_type> const & r, scalar_type &d )
	{
		const scalar_type EPSILON = 1e-20;

		vectorial_type p;

		scalar_type k = plane.n * r.dire;           // Plane intersection.
		if (k< EPSILON && k>-EPSILON)
			return false;

		d = (plane.d - plane.n * r.orig) / k;
		p = r.orig + r.dire*d;

		scalar_type b;

		vectorial_type p0 = v[0]->P();
		vectorial_type p1 = v[1]->P();
		vectorial_type p2 = v[2]->P();

		switch(flags & (NORMX|NORMY|NORMZ)){
		  case NORMX:
			b= edge[1][1]*(p[2]-p1[2]) - edge[1][2]*(p[1]-p1[1]); if(CheckVal(b)) return false;
			b= edge[2][1]*(p[2]-p2[2]) - edge[2][2]*(p[1]-p2[1]); if(CheckVal(b)) return false;
			b= edge[0][1]*(p[2]-p0[2]) - edge[0][2]*(p[1]-p0[1]); if(CheckVal(b)) return false;
			break;
		  case NORMY:
			b= edge[1][2]*(p[0]-p1[0]) - edge[1][0]*(p[2]-p1[2]); if(CheckVal(b)) return false;
			b= edge[2][2]*(p[0]-p2[0]) - edge[2][0]*(p[2]-p2[2]); if(CheckVal(b)) return false;
			b= edge[0][2]*(p[0]-p0[0]) - edge[0][0]*(p[2]-p0[2]); if(CheckVal(b)) return false;
			break;
		  case NORMZ:
			b= edge[1][0]*(p[1]-p1[1]) - edge[1][1]*(p[0]-p1[0]); if(CheckVal(b)) return false;
			b= edge[2][0]*(p[1]-p2[1]) - edge[2][1]*(p[0]-p2[0]); if(CheckVal(b)) return false;
			b= edge[0][0]*(p[1]-p0[1]) - edge[0][1]*(p[0]-p0[0]); if(CheckVal(b)) return false;
			break;
		}

		return true;
	}

#endif

	// Intersect ray with triangle.
	// returns barycientric coordinate and dist
	bool Intersect( Line3<scalar_type> const & ray, scalar_type & dist, scalar_type &a, scalar_type &b) {
		//static double a,b;
	  return Intersection( ray, v[0]->P(), v[1]->P(), v[2]->P(), a,b, dist );
	};

	/// return the index [0..2] of a vertex in a face
	inline int VertexIndex( const FVTYPE * w ) const
	{
			 if( v[0]==w ) return  0;
		else if( v[1]==w ) return  1;
		else if( v[2]==w ) return  2;
		else               return -1;
	}

		/// Return the texture distorsion
	scalar_type texture_distorsion( int nt=0 ) const
	{
#ifndef __VCGLIB_FACE_WT
		assert(0);
#endif
		scalar_type e = 0;
		for(int nz=0;nz<3;++nz)
		{
			scalar_type l0 = Distance( V1(nz)->P(),V(nz)->P() );
			scalar_type l1 = Distance( WT((nz+1)%3).t(nt), WT(nz).t(nt) );
			scalar_type dl = l1-l0;
			e += dl*dl;
		}
		e = (e/3)/ Area();
		return e;
	}

		/// Return the texture distorsion
	scalar_type hoppe_distorsion( int nt=0 ) const
	{
#ifndef __VCGLIB_FACE_WT
		assert(0);
#endif
		scalar_type A = Area();

		scalar_type s1 = WT(0).u(nt);
		scalar_type t1 = WT(0).v(nt);
		scalar_type s2 = WT(1).u(nt);
		scalar_type t2 = WT(1).v(nt);
		scalar_type s3 = WT(2).u(nt);
		scalar_type t3 = WT(2).v(nt);

		vectorial_type Ss = (V(0)->P()*(t2-t3)+V(1)->P()*(t3-t1)+V(2)->P()*(t1-t2))/(2*A);
		vectorial_type St = (V(0)->P()*(s3-s2)+V(1)->P()*(s1-s3)+V(2)->P()*(s2-s1))/(2*A);

		scalar_type a = Ss*Ss;
		scalar_type b = Ss*St;
		scalar_type c = St*St;

		return sqrt((a+c)/2+sqrt((a-c)*(a-c)+4*b*b));
	}

}; //end Class




}	 // end namespace


#endif

