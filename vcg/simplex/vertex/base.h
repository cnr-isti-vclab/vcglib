/*#***************************************************************************
 * VertexBase.h                                                         o o    *
 *                                                                  o     o  *
 * Visual Computing Group                                           _  O  _  *
 * IEI Institute, CNUCE Institute, CNR Pisa                          \/)\/   *
 *                                                                  /\/|     *
 * Copyright(C) 1999 by Paolo Cignoni, Paolo Pingi, Claudio Rocchini   |     *
 * All rights reserved.                                                \     *
 *                                                                           *
 * Permission  to use, copy, modify, distribute  and sell this  software and *
 * its documentation for any purpose is hereby granted without fee, provided *
 * that  the above copyright notice appear  in all copies and that both that *
 * copyright   notice  and  this  permission  notice  appear  in  supporting *
 * documentation. the author makes  no representations about the suitability *
 * of this software for any purpose. It is provided  "as is" without express *
 * or implied warranty.                                                      *
 *                                                                           *
 * NOTE THAT THIS FILE SHOULD NOT DIRECTL BE INCLUDED                        *
 * It is automatically included by Mesh.h                                    *
 *                                                                           *
 ***************************************************************************#*/
/*#**************************************************************************
  History

 2000	Jan 31 First Working release
 2000 Feb 04 USER0 added to the flag's enum
					11 Aggiunta funzione InitIMark()
			Jun 13 Aggiunta adiacenze vertice faccia				 
			Jun 26 Aggiunto cP() per forzare l'accesso costante alle coord.
					27 Vertex e' stato templatato anche sul tipo della faccia per
						 far resitutire a Fp il tipo giusto. Tale tipo ha un valore di 
						 default = a DUMMYFACETYPE, per permettere l'uso della classe
						 secondo lo stile precedente. 
					28 Aggiunti Flag NOTBORDER e NOTMANIFOLD 
      Sep 27 Aggiunto cN() per forzare l'accesso costante alla Normale.
			Oct 31 Tolti i flag del bordo per vertice e le funzioni collegate 
						 che erano inutili e scorretti (pc)
      Nov 01 Aggiunte assert(0) e commenti se si tenta di usare vf topology 
						 senza averla
      Nov 30 Cambiato il tipo flags da int in unsigned int (??);
			Dec 18 Aggiunto NewBitFlag() e DeleteBitFlag()
 2001 Jan 03 Aggiunto Supervisor_Normal e assert(0) in normal
      Jan 27 Aggiunta flags BORDER (C.R.)
		  Feb 16 Aggiunto colore
	           Aggiunte coordinate texture
				     Corretti public e protected per le facce;
			Mar 08 Aggiunto assert(0) se si cerca di accedere a C() in 
						 modo non costante e il vertice non ha il colore (pc)
          20 Corretto VF per i casi sbagliati (	return *(VFTYPE **)flags;  invece di 	return (void *)this; 	) 
  	  May 16 Aggiunta gestione qualita' (CR)
			       Aggiunta gestione OBJ per qualita' (CR)		
		  Jun 12 Aggiunte assert(0) ai lettori di dati inesistenti
				  13 Cambiato scalare coordinata texture default
					19 Commentate funzioni normal
						 Modificate funzione N, cN , in modo da rispettare lo standard
 			Jul 27 Aggiunto supervisor_flags const (pc)
			Sep 28 Aggiunto Supervisor_N() (pc)
			Oct 16 Tolte un paio di parentesi a DeleteBitFlag (facevano un warning nel compilatore intel) (pc)
 2002 Jan Modificato in operator [] v[i] con V(i) (PP)
			Dic Tolto IsMF() (gano)
 2003 Mag Aggiunto dati per TensorMass(particle)
			July 10: Add 2 properties to the particle (in case of an explicit FEM).
                        Damping and fixed status.   (cesar)
			Oct 7, 2003  Damping and fixed status of the particle also for TensorMass  
			Oct 21  Aggiunte IsUserBit(USERBIT),ClearUserBit(..) e SetUserBit(..) (gano)
****************************************************************************/

/*
People should subclass his vertex class from these one...
*/

#ifndef VERTEX_TYPE 
#pragma message("\nYou should never directly include this file\n")
#else


class DUMMYFACETYPE;

namespace vcg {

/** @name Vertex
    Class Vertex.
    This is the base class for definition of a vertex of the mesh.
	@param FLTYPE (Template Parameter) Specifies the scalar field of the vertex coordinate type.
	@param VFTYPE (Template Parameter) Specifies the type for the face, needed only for VF adjacency.
 */
template <class FLTYPE, class VFTYPE = DUMMYFACETYPE, class TCTYPE = TCoord<float,1> > class VERTEX_TYPE
{
public:

	/// The scalar type
	typedef FLTYPE         scalar_type;
	/// The coordinate type
	typedef Point3<FLTYPE> coord_type;
	/// The type base of the vertex
	typedef VERTEX_TYPE    vertex_base;
	typedef VFTYPE    face_type;
protected:
	/// Spatial coordinates of the vertex
	coord_type p;
	/// This are the flags of vertex, the default value is 0
	int flags;		

		// Definizione texture
#ifdef __VCGLIB_VERTEX_T
	TCTYPE t;
#endif

public:
	TCTYPE & T()
	{
#ifdef __VCGLIB_VERTEX_T
		return t;
#else
		assert(0);
		return *(TCTYPE*)(&flags);
#endif
	}

	const TCTYPE & T() const
	{
#ifdef __VCGLIB_VERTEX_T
		return t;
#else
		assert(0);
		return *(TCTYPE*)(&flags);
#endif
	}

		// Definizione del colore

#ifdef __VCGLIB_VERTEX_C
protected:
	ColorUB c;
#endif

public:
	ColorUB & C()
	{
#ifdef __VCGLIB_VERTEX_C
		return c;
#else
		assert(0);
		return *(ColorUB*)(&flags);
#endif
	}

	const ColorUB & C() const
	{
#ifdef __VCGLIB_VERTEX_C
		return c;
#else
		return ColorUB(ColorUB::White);
#endif
	}

		// Definizione Qualita'
#ifdef __VCGLIB_VERTEX_Q
protected:
	float quality;
#endif

public:
	float & Q()
	{
#ifdef __VCGLIB_VERTEX_Q
		return quality;
#else
		assert(0);
		return *(float*)(&flags);
#endif
	}

	const float & Q() const
	{
#ifdef __VCGLIB_VERTEX_Q
		return quality;
#else
		return 1;
#endif
	}

// Field to contains the index of the object in the CONTAINER
protected:
/*#*********************************
* Puntatore ad una faccia di v star*
***********************************/
#if ((defined __VCGLIB_VERTEX_A) || (defined __VCGLIB_VERTEX_AS)) 
	/// Puntatore ad una faccia appartenente alla stella del vertice, implementa l'adiacenza vertice-faccia
	VFTYPE *fp;
	int zp;
#endif

/*#**************
*  Mark Members *	
*****************/
#ifdef __VCGLIB_VERTEX_M
	/// The incremental vertex mark
	int imark;
#endif // Mark

/*#**************
*  Normal Members *
*****************/
#ifdef __VCGLIB_VERTEX_N
	/// The normal to the vertex
	coord_type n;
#endif // Normal

public:
	/// Return the spatial coordinate of the vertex
	inline coord_type & P()
	{
	  assert( (flags & DELETED) == 0 );
		assert( (flags & NOTREAD) == 0 );
		assert( (flags & NOTWRITE) == 0 );
		return p;
	}

	/// Return the constant spatial coordinate of the vertex
	inline const coord_type & P() const
	{
		assert( (flags & DELETED) == 0 );
		assert( (flags & NOTREAD) == 0 );
		return p;
	}

	/// Return the constant spatial coordinate of the vertex
	inline const coord_type & cP() const
	{
		assert( (flags & DELETED) == 0 );
		assert( (flags & NOTREAD) == 0 );
		return p;
	}

	/// Return the spatial coordinate of the vertex, senza effettuare controlli sul flag
	inline coord_type & Supervisor_P()
	{
		return p;
	}

	/// Return the constant spatial coordinate of the vertex, senza effettuare controlli sul flag
	inline const coord_type & Supervisor_P() const
	{
		return p;
	}



	/// Return the Normal of the vertex
	inline coord_type & Normal()
	{
	  assert( (flags & DELETED) == 0 );
		assert( (flags & NOTREAD) == 0 );
		assert( (flags & NOTWRITE) == 0 );
#ifdef __VCGLIB_VERTEX_N
		return n;
#else
		assert(0);
		return *(coord_type *)this;
#endif
	}


#if 0 // Inizio commentatura vecchio stile normali

	/// Return the constant normal of the vertex
	inline const coord_type & Normal() const
	{
		assert( (flags & DELETED) == 0 );
		assert( (flags & NOTREAD) == 0 );
#ifdef __VCGLIB_VERTEX_N
		return n;
#else
		assert(0);
		return *(coord_type *)this;
#endif
	}

#endif // Fine commentatura vecchio stile normali

/// Return the Normal of the vertex
	inline coord_type & Supervisor_N()
	{
#ifdef __VCGLIB_VERTEX_N
		return n;
#else
		assert(0);
		return *(coord_type *)this;
#endif
	}

	/// Return the constant normal of the vertex
	inline const coord_type & Supervisor_N() const
	{
#ifdef __VCGLIB_VERTEX_N
		return n;
#else
		assert(0);
		return *(coord_type *)this;
#endif
	}



	/// Return the vector of flags
	inline int & Flags ()
	{
			assert( (flags & DELETED) == 0 );
			assert( (flags & NOTREAD) == 0 );
			return flags;
	}

	/// Return the vector of flags, senza effettuare controlli sui bit
	inline int & Supervisor_Flags ()
	{
			return flags;
	}
	inline const int Supervisor_Flags() const
	{
		return flags;
	}



	/// Return the vertex normal
	inline coord_type & N()
	{
		assert( (flags & DELETED) == 0 );
		assert( (flags & NOTREAD) == 0 );
		assert( (flags & NOTWRITE) == 0 );
#ifdef __VCGLIB_VERTEX_N
		return n;
#else
		assert(0);
		return *(coord_type *)this;
#endif
	}

	/// Return the constant vertex normal
	inline const coord_type & N() const
	{
		assert( (flags & DELETED) == 0 );
		assert( (flags & NOTREAD) == 0 );
#ifdef __VCGLIB_VERTEX_N
		return n;
#else
		assert(0);
		return *(coord_type *)this;
#endif
	}

	inline const coord_type  cN() const
	{
		assert( (flags & DELETED) == 0 );
		assert( (flags & NOTREAD) == 0 );
#ifdef __VCGLIB_VERTEX_N
		return n;
#else
		return coord_type(0,0,0);
#endif
	}


#ifdef __VCGLIB_VERTEX_M
	/// This function return the vertex incremental mark
	inline int & IMark()
	{
		assert( (flags & DELETED) == 0 );
		assert( (flags & NOTREAD) == 0 );
		assert( (flags & NOTWRITE) == 0 );
		return imark;
	}

	/// This function return the constant vertex incremental mark
	inline const int & IMark() const
	{
		assert( (flags & DELETED) == 0 );
		assert( (flags & NOTREAD) == 0 );
		return imark;
	}
#endif

	/// Initialize the imark system of the vertex
	inline void InitIMark()
	{
#ifdef __VCGLIB_VERTEX_M
		imark = 0;
#endif
	}



inline VFTYPE * & Fp()
	{
#if ((defined __VCGLIB_VERTEX_A) || (defined __VCGLIB_VERTEX_AS))
		  return fp;
#else
    assert(0);// you are probably trying to use VF topology in a vertex without it
		return *((VFTYPE **)(flags));  
#endif
	}

inline const VFTYPE * & Fp() const
	{
#if ((defined __VCGLIB_VERTEX_A) || (defined __VCGLIB_VERTEX_AS))
		  return fp;
#else
		assert(0);// you are probably trying to use VF topology in a vertex without it
		return (VFTYPE *)this;
#endif
	}

inline int & Zp()
	{
#if ((defined __VCGLIB_VERTEX_A) || (defined __VCGLIB_VERTEX_AS))

		  return zp;
#else
    assert(0);// you are probably trying to use VF topology in a vertex without it
		return flags;
#endif
	}

inline const int & Zp() const
	{
#if ((defined __VCGLIB_VERTEX_A) || (defined __VCGLIB_VERTEX_AS))
		  return zp;
#else
		assert(0);// you are probably trying to use VF topology in a vertex without it
		return (void *)this;
#endif
	}

#ifdef __PARTICLE
        // variable declaration
        /** external force acting on the particle */
		coord_type extForce;
		/** internal force acting on the particle */
		coord_type intForce;
        /** mass of the particle */
         double mas;
        /** velocity of the particle */
		coord_type vel;
		/** accelleration of the particle */
		coord_type acc;
        /** current position of the particle */
		coord_type pos;

	    /** damping of the particle */
	    coord_type _damping;

       /** Fixed particle.  */
        bool _pointFixed;

		void computeAccelleration()
		{   
			 // acc=( (  extForce  + intForce  )/mas);
			  acc=( (  extForce  + intForce + _damping )/mas);
		};

        void resImpFor(){extForce = coord_type(0.0,0.0,0.0);}

		/** ComputeExternal forces */ 
		void computeExternalForces( coord_type value)
			{
				extForce = value;
			}

		bool fixedParticle( void ) { return _pointFixed; }
		void fixParticle( bool value ) { _pointFixed = value;}
#endif

	enum {
		OBJ_TYPE_N =  0x0001,
		OBJ_TYPE_M =  0x0002,
		OBJ_TYPE_A =  0x0004,
		OBJ_TYPE_AS = 0x0008,
		OBJ_TYPE_C  = 0x0010,
		OBJ_TYPE_T  = 0x0020,
		OBJ_TYPE_Q  = 0x0040,
	};

	enum {
		OBJ_TYPE = 
#ifdef __VCGLIB_VERTEX_N
		OBJ_TYPE_N |
#endif
#ifdef __VCGLIB_VERTEX_M
		OBJ_TYPE_M |
#endif
#ifdef __VCGLIB_VERTEX_A
		OBJ_TYPE_A |
#endif
#ifdef __VCGLIB_VERTEX_AS
		OBJ_TYPE_AS |
#endif
#ifdef __VCGLIB_VERTEX_C
		OBJ_TYPE_C |
#endif
#ifdef __VCGLIB_VERTEX_T
		OBJ_TYPE_T |
#endif
#ifdef __VCGLIB_VERTEX_Q
		OBJ_TYPE_Q |
#endif
		0
	};
	

	enum { 
		// This bit indicate that the vertex is deleted from the mesh
		DELETED    = 0x0001,		// cancellato
		// This bit indicate that the vertex of the mesh is not readable
		NOTREAD    = 0x0002,		// non leggibile (ma forse modificabile) 
		// This bit indicate that the vertex is not modifiable
		NOTWRITE   = 0x0004,		// non modificabile (ma forse leggibile) 
		// This bit indicate that the vertex is modified
		MODIFIED   = 0x0008,		// modificato 
		// This bit can be used to mark the visited vertex
		VISITED    = 0x0010,		// Visited  
		// This bit can be used to select 
		SELECTED   = 0x0020,		// Selection flag
		// Border Flag
		BORDER     = 0x0100,
		// First user bit
		USER0      = 0x0200			// Fisrt user bit
			};

/*
Queste funzioni servono per ottenere a runtime un bit per i flag

*/	
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

	/** Return the i-th spatial value of the vertex coordinate.
	    @param i Index of the spatial vertex coordinate (x=0 y=1 z=2).
	 */
	inline FLTYPE & operator [] ( const int i ){
			assert(i>=0 && i<3);
			return P().V(i);
	}
	/** Return the i-th spatial value of the const vertex coordinate.
	    @param i Index of the spatial vertex coordinate (x=0 y=1 z=2).
	 */
	inline const FLTYPE & operator [] ( const int i ) const {
			assert(i>=0 && i<3);
			return P().V(i);
	}
	/// Operator to compare two vertices using lexicographic order
	inline bool operator < ( const VERTEX_TYPE & ve) const {
		return p < ve.p;
		}
	inline VERTEX_TYPE() {
#ifdef _DEBUG 
		flags=0;
#endif
	};

	/// This function checks if the vertex is deleted
	bool IsD() const {return (flags & DELETED) != 0;}
	/// This function checks if the vertex is readable
	bool IsR() const {return (flags & NOTREAD) == 0;}
	/// This function checks if the vertex is modifiable
	bool IsW() const {return (flags & NOTWRITE)== 0;}
	/// This funcion checks whether the vertex is both readable and modifiable
	bool IsRW() const {return (flags & (NOTREAD | NOTWRITE)) == 0;}
	/// This function checks if the vertex is Modified
	bool IsM() const {return (flags & MODIFIED)!= 0;}
	/// This function checks if the vertex is marked as visited
	bool IsV() const {return (flags & VISITED) != 0;}
	/// This function checks if the vertex is selected
	bool IsS() const {return (flags & SELECTED) != 0;}
	/// This function checks if the vertex is readable
	bool IsB() const {return (flags & BORDER) != 0;}
//	bool IsMF() const {return (flags & NOTMANIFOLD) == 0;}

	/// This function checks if the vertex is deleted from the mesh
	bool IsDeleted() const {return IsD();}
	/// This function checks if the vertex is readable
	bool IsReadable() const {return IsR();}
	/** Set the flag value
		@param flagp Valore da inserire nel flag
	*/
	void SetFlags(int flagp) {flags=flagp;}

	/// This function deletes the vertex from the mesh
	void SetD() {flags |=DELETED;}
	/// This funcion execute the inverse operation of SetD()
	void ClearD() {flags &=(~DELETED);}
	/// This function marks the vertex as modified. It's necessary to mark all modified vertex to have a consistent mesh
	void SetM() {flags |=MODIFIED;}
	/// This function marks the vertex as not modified
	void ClearM() {flags &=(~MODIFIED);}
	/// This function marks the vertex as readable
	void SetR() {flags &=(~NOTREAD);}
	/// This function marks the vertex as not readable
	void ClearR() {flags |=NOTREAD;}
	/// This function marks the vertex as writable
	void ClearW() {flags |=NOTWRITE;}
	/// This function marks the vertex as not writable
	void SetW() {flags &=(~NOTWRITE);}
	/// This funcion marks the vertex as visited
	void SetV() {flags |=VISITED;}
	/// This function marks the vertex as not visited. This flag, initially, is setted to random value, therefore, to the beginnig of every function it is necessary to clean up the flag
	void ClearV() {flags &=(~VISITED);}
	/// This function select the vertex
	void SetS()		{flags |=SELECTED;}
	/// This funcion execute the inverse operation of SetS()
	void ClearS()	{flags &= ~SELECTED;}
	void SetB()		{flags |=BORDER;}
	void ClearB()	{flags &=~BORDER;}
	
	/// This function checks if the given user bit is true
	bool IsUserBit(int userBit){return (flags & userBit) != 0;}
	/// This function set  the given user bit 
	void SetUserBit(int userBit){flags |=userBit;}
	/// This function clear the given user bit 
	void ClearUserBit(int userBit){flags &= (~userBit);}
};


}	 // end namespace
#endif

/*
 * mode: c++
 * tab-width: 3
 * c-basic-offset: 3
 */
