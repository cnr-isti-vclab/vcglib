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

#pragma warning( disable : 4804 )

/*
People should subclass his vertex class from these one...
*/

#ifndef __VCG_MESH
#define __VCG_MESH

#include <assert.h>
#include <list>
#include <vector>
#include <set>
#include <map>
#include <stack>
#include <algorithm>
#include <iterator>
//#include <vcg/Mesh/Selection.h>

#include <vcg/TriTriIntersection.h>
#include <vcg/tools/plylib.h>

namespace vcg {

/** Class Mesh.
    This is class for definition of a mesh.
		@param STL_VERT_CONT (Template Parameter) Specifies the type of the vertices container any the vertex type.
		@param STL_FACE_CONT (Template Parameter) Specifies the type of the faces container any the face type.
 */
template < class STL_VERT_CONT, class STL_FACE_CONT >
class Mesh{
	public:
	/// The face container
	typedef STL_FACE_CONT face_container;
	/// The face container
	typedef STL_VERT_CONT vertex_container;
	/// The vertex type 
	typedef typename STL_VERT_CONT::value_type MVTYPE;
	/// The face type 
	typedef typename STL_FACE_CONT::value_type MFTYPE;
	/// The scalar type
	typedef typename MVTYPE::scalar_type MCTYPE;
	/// The type of the vectors
	typedef typename MFTYPE::vectorial_type vectorial_type;
	/// The type of the scalars
	typedef MCTYPE scalar_type;
	/// The vertex type
	typedef MVTYPE vertex_type;
	/// Tipo vertice originario
	typedef typename MVTYPE::vertex_base vertex_base;
	/// The face type
	typedef MFTYPE face_type;
	/// The type of vertex iterator
	typedef typename STL_VERT_CONT::iterator vertex_iterator;
	/// The type of face iterator
	typedef typename STL_FACE_CONT::iterator face_iterator;
	/// The type of constant vertex iterator
	typedef typename STL_VERT_CONT::const_iterator const_vertex_iterator;
	/// The type of constant face iterator
	typedef typename STL_FACE_CONT::const_iterator const_face_iterator;
	/// The vertex pointer type
	typedef MVTYPE * vertex_pointer;
	/// The face pointer type
	typedef MFTYPE * face_pointer;
	/// The type of the constant vertex pointer
	typedef const MVTYPE * const_vertex_pointer;
	/// The type of the constant face pointer
	typedef const MFTYPE * const_face_pointer;
	/// The vertex base pointer type
	typedef typename MVTYPE::vertex_base * vertex_base_pointer;
	/// The type of the constant vertex base pointer
	typedef const typename MVTYPE::vertex_base * const_vertex_base_pointer;
	/// The face base pointer type
	typedef typename MFTYPE::face_base face_base;
	typedef typename MFTYPE::face_base * face_base_pointer;
	/// The type of the constant face base pointer
	typedef const typename MFTYPE::face_base * const_face_base_pointer;
	/// The mesh type
	typedef Mesh<STL_VERT_CONT,STL_FACE_CONT> MMTYPE;
	/// The edge type for FF topology
	typedef FEdgePosB<MFTYPE> fedgepos_type;
	/// The edge type for VF topology
	typedef VEdgePosB<MFTYPE> vedgepos_type;
	/// The half edge type
	typedef HEdgePosB<MFTYPE> hedgepos_type;
	/// The half edge type with the normal
	typedef HEdgePosBN<MFTYPE> hedgeposn_type;
	/// The ear type
	typedef Ear<MFTYPE> ear_type;
	/// The Box3 type
	typedef Box3<MCTYPE> BOX_TYPE;

	/// Set of vertices 
	STL_VERT_CONT vert;
	/// Real number of vertices
	int vn;
	/// Set of faces
	STL_FACE_CONT face;
	/// Real number of faces
	int fn;
	/// Bounding box of the mesh
	Box3<MCTYPE> bbox;
	/// Internal status
	int status;
	
  /// Nomi di textures
	vector<string> textures;
	vector<string> normalmaps;

		/// La camera
	Camera<scalar_type> camera;

		/// Il colore della mesh
private:
	ColorUB c;
public:

	inline const ColorUB & C() const
	{
		return c;
	}

	inline ColorUB & C()
	{
		return c;
	}


	/// Default constructor
	Mesh()
	{
		fn = vn = 0;
		imark = 0;
	}

	inline int MemUsed() const
	{
		return sizeof(MMTYPE)+sizeof(MVTYPE)*vert.size()+sizeof(MFTYPE)*face.size();
	}

	inline int MemNeeded() const
	{
		return sizeof(MMTYPE)+sizeof(MVTYPE)*vn+sizeof(MFTYPE)*fn;
	}



/// Function to destroy the mesh
MMTYPE & Clear()
{
	vert.clear();
	face.clear();
	textures.clear();
	normalmaps.clear();
	vn = 0;
	fn = 0;
	return *this;
}

/* Funzioni di info sulle caratteristiche della mesh */ 

static bool HasPerVertexNormal()  { return bool(vertex_type::OBJ_TYPE & (vertex_type::OBJ_TYPE_N)); }
static bool HasPerVertexColor()   { return bool(vertex_type::OBJ_TYPE & (vertex_type::OBJ_TYPE_C)); }
static bool HasPerVertexMark()    { return bool(vertex_type::OBJ_TYPE & (vertex_type::OBJ_TYPE_M)); }
static bool HasPerVertexQuality() { return bool(vertex_type::OBJ_TYPE & (vertex_type::OBJ_TYPE_Q)); }
static bool HasPerVertexTexture() { return bool(vertex_type::OBJ_TYPE & (vertex_type::OBJ_TYPE_T)); }

static bool HasPerFaceColor()     { return bool(face_type::OBJ_TYPE & (face_type::OBJ_TYPE_C)); }
static bool HasPerFaceNormal()    { return bool(face_type::OBJ_TYPE & (face_type::OBJ_TYPE_N)); }
static bool HasPerFaceMark()      { return bool(face_type::OBJ_TYPE & (face_type::OBJ_TYPE_M)); }
static bool HasPerFaceQuality()   { return bool(face_type::OBJ_TYPE & (face_type::OBJ_TYPE_Q)); }

static bool HasPerWedgeColor()    { return bool(face_type::OBJ_TYPE & (face_type::OBJ_TYPE_WC)); }
static bool HasPerWedgeNormal()   { return bool(face_type::OBJ_TYPE & (face_type::OBJ_TYPE_WN)); }
static bool HasPerWedgeTexture()  { return bool(face_type::OBJ_TYPE & (face_type::OBJ_TYPE_WT)); }

static bool HasFFTopology()       { return bool(face_type::OBJ_TYPE & (face_type::OBJ_TYPE_A)) || HasSTopology();  }
static bool HasVFTopology()       { return bool(face_type::OBJ_TYPE & (face_type::OBJ_TYPE_V)) || HasSTopology(); }
static bool HasSTopology()        { return bool(face_type::OBJ_TYPE & (face_type::OBJ_TYPE_S)); }
static bool HasTopology()         { return HasFFTopology() || HasVFTopology(); }




/// Initialize the imark-system of the faces
void InitFaceIMark()
{
	face_iterator f;
	
	for(f=face.begin();f!=face.end();++f)
		if( !(*f).IsDeleted() && (*f).IsR() && (*f).IsW() )
			(*f).InitIMark();
}

/// Initialize the imark-system of the vertices
void InitVertexIMark()
{
	vertex_iterator vi;

	for(vi=vert.begin();vi!=vert.end();++vi)
		if( !(*vi).IsDeleted() && (*vi).IsRW() )
			(*vi).InitIMark();
}

// Warning assignament should take a const mesh in input
/** Assignment operator for mesh. The mesh content is losed.
*/
inline MMTYPE & operator = (MMTYPE & m )
{
	Clear();
	SelectedMerge(m,true);
	return *this;
}
/// The incremental mark
int imark;

/** Check if the vertex incremental mark matches the one of the mesh. 
	@param v Vertex pointer
*/
inline bool IsMarked( MVTYPE * const v ) const { return v->IMark() == imark; }
/** Check if the face incremental mark matches the one of the mesh. 
	@param v Face pointer
*/
inline bool IsMarked( MFTYPE * const f ) const { return f->IMark() == imark; }
/** Set the vertex incremental mark of the vertex to the one of the mesh.
	@param v Vertex pointer
*/
inline void Mark( MVTYPE * const v ) const { v->IMark() = imark; }
/** Set the face incremental mark of the vertex to the one of the mesh.
	@param v Vertex pointer
*/
inline void Mark( MFTYPE * const f ) const { f->IMark() = imark; }
/// Unmark the mesh
inline void UnMarkAll() { ++imark; }


/** Function to add n vertices to the mesh. The second parameter hold a vector of 
	pointers to pointer to elements of the mesh that should be updated after a 
	possible vector realloc. 
	@param n Il numero di vertici che si vuole aggiungere alla mesh.
	@param local_var Vettore di variabili locali che rappresentano puntatori a vertici. 
	restituisce l'iteratore al primo elemento aggiunto.
*/
vertex_iterator AddVertices(int n, vector<vertex_base **> &local_var)
{
	vertex_iterator oldbegin, newbegin;
	oldbegin = vert.begin();
  vertex_iterator last=vert.end();
	if(vert.empty()) last=0;  // if the vector is empty we cannot find the last valid element
	else --last;
	unsigned int siz=0;
#ifdef __STL_CONFIG_H	
if(last!=0) distance(vert.begin(),last,siz);
#else
if(last!=0) siz=distance(vert.begin(),last);
#endif
	for(int i=0; i<n; ++i)
	{
		vert.push_back(MVTYPE());
		vert.back().Supervisor_Flags() = 0;
	}
	vn+=n;
	newbegin = vert.begin();
	if(newbegin != oldbegin)
		{
			face_iterator f;
			for (f=face.begin(); f!=face.end(); ++f)
				if(!(*f).IsD())
				for(int k=0; k<(*f).size(); ++k)
					(*f).V(k) = (*f).V(k)-&*oldbegin+&*newbegin;
			for(int j=0; j<local_var.size(); ++j)
				if((*local_var[j]) !=0 ) *local_var[j] = *local_var[j]-&*oldbegin+&*newbegin;

			// deve restituire l'iteratore alla prima faccia aggiunta;
			// e poiche' lo spazio e' cambiato si ricalcola last da zero  
			if(last!=0) 
			{ 
				last = vert.begin(); 
				advance(last,siz+1);
			}
			else last=vert.begin(); 
		}
	else 
	{ 
		// se non e'cambiato lo spazio (vector abbastanza grande o lista)
		if(last==0) last = vert.begin(); // se il vettore era vuoto si restituisce begin
		           else advance(last,1); // altrimenti il primo dopo quello che era in precedenza l'ultimo valido.
	}
	return last;
}

vertex_iterator AddVertices(int n)
{
	vector<vertex_base **> local_var;
	return AddVertices(n,local_var);
}

/** Function to add n faces to the mesh.
	@param n Il numero di facce che si vuole aggiungere alla mesh
*/
face_iterator AddFaces(int n)
{
	vector<face_base **> local_var;
	return AddFaces(n,local_var);
}
/** Function to add n faces to the mesh. 
  NOTA: Aggiorna fn;
	The second parameter hold a vector of 
	pointers to pointer to elements of the mesh that should be updated after a 
	possible vector realloc.
	@param n Facce da aggiungere
	@param local_var Vettore di variabili locali che rappresentano puntatori a facce, occorre, 
	perche' questi valori siano consistenti, aggiornarli ogni qual volta venga eseguito un resize
	del contenitore delle facce.
*/
face_iterator AddFaces(int n, vector<face_base **> &local_var)
{
	face_iterator oldbegin, newbegin;
	oldbegin = face.begin();
	face_iterator last=face.end();
	if(face.empty()) last=0;
	            else last--;

	unsigned int siz=0;
#ifdef __STL_CONFIG_H	
	if(last!=0) distance(face.begin(),last,siz);
#else
	if(last!=0) siz=distance(face.begin(),last);
#endif
	MFTYPE dum;
	dum.Supervisor_Flags()=0;
	for(int i=0; i<n; ++i)
		face.push_back(dum);
	
	fn+=n;
	newbegin = face.begin();
	if(newbegin != oldbegin)// se e' cambiato lo spazio (vector abbastanza grande o lista)
	{
		if(MFTYPE::OBJ_TYPE & MFTYPE::OBJ_TYPE_A) 
		{
			face_iterator f;
			for (f=face.begin(); f!=face.end(); ++f)
				for(int k=0; k<(*f).size(); ++k)if(!(*f).IsD())
					(*f).F(k) = (*f).F(k)-&*oldbegin+&*newbegin;
		}
		vector<face_base **>::iterator jit;
		for(jit=local_var.begin(); jit!=local_var.end(); ++jit)
			if((**jit) !=0 ) **jit = **jit-&*oldbegin+&*newbegin;
		
		// deve restituire l'iteratore alla prima faccia aggiunta;
		if(last!=0) 
		{ 
			last = face.begin(); 
			advance(last,siz+1);
		}
		else last=face.begin(); 
	}
	else // 
	{ assert(newbegin == oldbegin);
		// se non e'cambiato lo spazio (vector abbastanza grande o lista)
		if(last==0) last = face.begin(); // se il vettore era vuoto si restituisce begin
		           else advance(last,1); // altrimenti il primo dopo quello che era in precedenza l'ultimo valido.
	}

	return last;

}


/// Calcolo del volume di una mesh chiusa
scalar_type Volume()
{
 
  face_iterator f;
  int j,k;
  scalar_type V = 0;
  vectorial_type T,N,B;
 
  for(f = face.begin(); f!=face.end(); ++f)
  {
	for(j = 0; j < 3; ++j)
	{
	  /*calcolo tangente, normale e binormale (6 volte)*/
	  k = (j+1)%3;
	  T = (*f).V(k)->P() - (*f).V(j)->P();
	  T.Normalize();
	  T = ( (*f).V( k     )->P() - (*f).V(j)->P() ) ^
	   ( (*f).V((k+1)%3)->P() - (*f).V(j)->P() ) ;
	  B.Normalize();
	  N = T ^ B;
   
	  vectorial_type pj = (*f).V(j)->P();
	  vectorial_type pk = (*f).V(k)->P();
 

	  V +=  (pj*  T )*(pj*N)*(pj*B);
	  V +=  (pk*(-T))*(pk*N)*(pk*B);
    }
  }
	return V/6;
}


};	// end class Mesh


}	 // end namespace


#endif

