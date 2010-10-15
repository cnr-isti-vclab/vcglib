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

#ifndef __VCG_POLYGON_COMPONENT
#define __VCG_POLYGON_COMPONENT

//#include <vector>
//#include <vcg/space/triangle3.h>
//#include <vcg/space/texcoord2.h>
//#include <vcg/space/color4.h>

namespace vcg {
  namespace face {
/*
Some naming Rules
All the Components that can be added to a vertex should be defined in the namespace vert:

*/

/*-------------------------- PolInfo -------------------------------------------*/
template <class T> class EmptyPolyInfo: public T {
protected:
	inline void SetVN(const int & n) {assert(0);}
public:
	typedef typename T::HEdgeType HEdgeType;
	typedef typename T::HEdgePointer HEdgePointer;

	/* Note: the destructor will not be called in general because there are no virtual destructors.
		Instead, the job of deallocating the memory will be done bu the edge allocator.
		This destructor is only done for those who istance a face alone (outside a mesh)
	*/
		static bool HasPolyInfo()   { return false; }
		inline void Alloc(const int & ns){T::Alloc(ns);};// it should be useless
		inline void Dealloc(){T::Dealloc();};// it should be useless

		// EmptyPFHAdj
		HEdgePointer       &FHp( )       { static typename T::HEdgePointer fp=0;  assert(0); return fp; }
		HEdgePointer const cFHp( ) const { static typename T::HEdgePointer const fp=0; return fp; }
		static bool HasFHAdjacency()   {   return false; }
};




template <class T> class PolyInfo: public T {
protected:
	inline void SetVN(const int & n) {_ns = n;}
public:
	PolyInfo(){ _ns = -1; }
	/* Note: the destructor will not be called in general because there are no virtual destructors.
		Instead, the job of deallocating the memory will be done bu the face allocator.
		This destructor is only done for those who istance a face alone (outside a mesh)
	*/
		static bool HasPolyInfo()   { return true; }
		inline const int &  VN() const { return _ns;}
		inline int Prev(const int & i){ return (i+(VN()-1))%VN();}
		inline int Next(const int & i){ return (i+1)%VN();}
		inline void Alloc(const int & ns){};
		inline void Dealloc(){};

		// EmptyPFHAdj
		typename T::HEdgePointer       &FHp(const int & )       { static typename T::HEdgePointer fp=0;  assert(0); return fp; }
		typename T::HEdgePointer const cFHp(const int & ) const { static typename T::HEdgePointer const fp=0; return fp; }
		static bool HasFHAdjacency()						{   return false; }
private:
	int _ns;
};

/*-------------------------- VERTEX ----------------------------------------*/ 
template <class T> class PEmptyFVAdj: public T {
public:
	typedef typename T::VertexType::CoordType CoordType;
	typedef typename T::VertexType::ScalarType ScalarType;
// typedef typename T::CoordType CoordType;
  inline typename T::VertexType *       & V( const int j ) 	    {	assert(0);		static typename T::VertexType *vp=0; return vp; }
  inline typename T::VertexType * const & V( const int j ) const {	assert(0);		static typename T::VertexType *vp=0; return vp; }
        inline typename T::VertexType * cV( const int j ) const {	assert(0);		static typename T::VertexType *vp=0; return vp;	}
	inline       typename T::CoordType & P( const int j ) 	    {	assert(0);		static typename T::CoordType coord(0, 0, 0); return coord;	}
	inline const typename T::CoordType & P( const int j ) const {	assert(0);		static typename T::CoordType coord(0, 0, 0); return coord;	}
	inline const typename T::CoordType &cP( const int j ) const	{	assert(0);		static typename T::CoordType coord(0, 0, 0); return coord;	}
	template <class LeftF>
	void ImportData(const LeftF & leftF) {T::ImportData(leftF);}
  	static bool HasFVAdjN()   { return false; }
	static void Name(std::vector<std::string> & name){T::Name(name);}
	inline void Alloc(const int & ns){T::Alloc();};
	inline void Dealloc(){T::Dealloc();}

};
template <class T> class PFVAdj: public T {
public:
	typedef typename  T::VertexType::CoordType CoordType;
 	typedef typename  T::VertexType::ScalarType ScalarType;
 	typedef typename  T::VertexType VertexType;

	PFVAdj(){_vpoly = NULL;}
 	inline typename T::VertexType *       & V( const int j ) 	     { assert(j>=0 && j<this->VN()); return  _vpoly[j]; }
	inline typename T::VertexType * const & V( const int j ) const { assert(j>=0 && j<this->VN()); return  _vpoly[j]; }
    inline typename T::VertexType *        cV( const int j ) const { assert(j>=0 && j<this->VN());	return _vpoly[j]; }


	/** Return the pointer to the ((j+1)%3)-th vertex of the face.
		@param j Index of the face vertex.
	 */
	inline        VertexType *       &  V0( const int j )       { return V(j);}
	inline        VertexType *       &  V1( const int j )       { return V((j+1)%this->VN());}
	inline        VertexType *       &  V2( const int j )       { return V((j+2)%this->VN());}
	inline const  VertexType * const &  V0( const int j ) const { return V(j);}
	inline const  VertexType * const &  V1( const int j ) const { return V((j+1)%this->VN());}
	inline const  VertexType * const &  V2( const int j ) const { return V((j+2)%this->VN());}
	inline const  VertexType * const & cV0( const int j ) const { return cV(j);}
	inline const  VertexType * const & cV1( const int j ) const { return cV((j+1)%this->VN());}
	inline const  VertexType * const & cV2( const int j ) const { return cV((j+2)%this->VN());}

	// Shortcut per accedere ai punti delle facce
	inline        CoordType & P( const int j ) 	    {	assert(j>=0 && j<this->VN());		return _vpoly[j]->P();	}
	inline const  CoordType & P( const int j ) const	{	assert(j>=0 && j<this->VN());		return _vpoly[j]->cP(); }
	inline const  CoordType &cP( const int j ) const	{	assert(j>=0 && j<this->VN());		return _vpoly[j]->cP(); }

	/// Shortcut per accedere ai punti delle facce
	inline        CoordType &  P0( const int j )       { return V(j)->P();}
	inline        CoordType &  P1( const int j )       { return V((j+1)%this->VN())->P();}
	inline        CoordType &  P2( const int j )       { return V((j+2)%this->VN())->P();}
	inline const  CoordType &  P0( const int j ) const { return V(j)->P();}
	inline const  CoordType &  P1( const int j ) const { return V((j+1)%this->VN())->P();}
	inline const  CoordType &  P2( const int j ) const { return V((j+2)%this->VN())->P();}
	inline const  CoordType & cP0( const int j ) const { return cV(j)->P();}
	inline const  CoordType & cP1( const int j ) const { return cV((j+1)%this->VN())->P();}
	inline const  CoordType & cP2( const int j ) const { return cV((j+2)%this->VN())->P();}

	template <class LeftF>
	void ImportData(const LeftF & leftF){  for(int i  =0; i < this->VN(); ++i) V(i) = NULL; T::ImportData(leftF);}
	inline void Alloc(const int & ns) {		
		if(_vpoly == NULL){this->SetVN(ns);
			_vpoly = new   typename T::VertexType*[this->VN()]; 
				for(int i = 0; i < this->VN(); ++i) _vpoly[i] = 0;
		}
		T::Alloc(ns);
	}
	inline void Dealloc() {	if(_vpoly!=NULL){
														delete [] _vpoly; 
														_vpoly = NULL;
														}
													T::Dealloc();
	}

  static bool HasFVAdjacency()   { return true; }
	static void Name(std::vector<std::string> & name){name.push_back(std::string("PFVAdj"));T::Name(name);}

  private:
	typename T::VertexPointer *_vpoly;
};

/*----------------------------- PVFADJ ------------------------------*/ 
template <class T> class EmptyPVFAdj: public T {
public:
 typedef typename T::VertexType VertexType;
	typedef int VFAdjType;
  typename T::FacePointer       &VFp(const int)       { static typename T::FacePointer fp=0;  assert(0); return fp; }
  typename T::FacePointer const cVFp(const int) const { static typename T::FacePointer const fp=0; return fp; }
  typename T::FacePointer       &FFp(const int)       { static typename T::FacePointer fp=0;  assert(0); return fp; }
  typename T::FacePointer const cFFp(const int) const { static typename T::FacePointer const fp=0; return fp; }
  char &VFi(const int j){static char z=0;  assert(0); return z;};
  char &FFi(const int j){static char z=0;  assert(0); return z;};
  const char &cVFi(const int j){static char z=0; return z;};
  const char &cFFi(const int j) const {static char z=0; return z;};
	template <class LeftF>
	void ImportData(const LeftF & leftF){ T::ImportData(leftF);}
	void Alloc(const int & ns){T::Alloc(ns);}
	void Dealloc(){T::Dealloc();}
  static bool HasVFAdjacency()   {   return false; }
  static bool HasFFAdjacency()   {   return false; }
  static bool HasFFAdjacencyOcc()   {   return false; }
  static bool HasVFAdjacencyOcc()   {   return false; }
  static void Name(std::vector<std::string> & name){T::Name(name);}

};

template <class T> class PVFAdj: public T {
public:

	PVFAdj(){_vfiP = NULL; _vfiP = NULL;}
	/* Note: the destructor will not be called in general because there are no virtual destructors.
		Instead, the job of deallocating the memory will be done bu the edge allocator.
		This destructor is only done for those who istance a face alone (outside a mesh)
	*/
 typedef typename T::VertexType VertexType;
 typedef typename T::FaceType FaceType;
  typename T::FacePointer       &VFp(const int j)        { assert(j>=0 && j<this->VN());  return _vfpP[j]; }
  typename T::FacePointer const  VFp(const int j) const  { assert(j>=0 && j<this->VN());  return _vfpP[j]; }
  typename T::FacePointer const cVFp(const int j) const  { assert(j>=0 && j<this->VN());  return _vfpP[j]; }
  char &VFi(const int j) {return _vfiP[j]; }
	template <class LeftF>
	void ImportData(const LeftF & leftF){T::ImportData(leftF);}
	inline void Alloc(const int & ns) { 
		if(_vfpP == NULL){
			this->SetVN(ns);
			_vfpP = new  FaceType*[this->VN()]; 
			_vfiP = new  char[this->VN()];
			for(int i = 0; i < this->VN(); ++i) {_vfpP[i] = 0;_vfiP = -1;}
		}
		T::Alloc(ns);
													
	}
	inline void Dealloc() {	if(_vfpP!=NULL){
														delete [] _vfpP; _vfpP = NULL;
														delete [] _vfiP; _vfiP = NULL;
													}
													T::Dealloc();
	}

  static bool HasVFAdjacency()      {   return true; }
  static bool HasVFAdjacencyOcc()   {   return false; }
  static void Name(std::vector<std::string> & name){name.push_back(std::string("PVFAdj"));T::Name(name);}

private:
  typename T::FacePointer *_vfpP ;    
  char *_vfiP ;    
};

/*----------------------------- FFADJ ------------------------------*/ 

template <class T> class PFFAdj: public T {
public:
 typedef typename T::FaceType FaceType;
	PFFAdj(){_ffpP = NULL; _ffiP = NULL; }
	typename T::FacePointer       &FFp(const int j)        { assert(j>=0 && j<this->VN());  return _ffpP[j]; }
  typename T::FacePointer const  FFp(const int j) const  { assert(j>=0 && j<this->VN());  return _ffpP[j]; }
  typename T::FacePointer const cFFp(const int j) const  { assert(j>=0 && j<this->VN());  return _ffpP[j]; }
  char        &FFi(const int j)       { return _ffiP[j]; }
  const char &cFFi(const int j) const { return _ffiP[j]; }
	
	template <class LeftF>
	void ImportData(const LeftF & leftF){T::ImportData(leftF);}
	inline void Alloc(const int & ns) {
		if( _ffpP == NULL){	
			this->SetVN(ns);
			_ffpP = new  FaceType*[this->VN()]; 
			_ffiP = new  char[this->VN()];
			for(int i = 0; i < this->VN(); ++i) {_ffpP[i] = 0;_ffiP[i] = 0;}
		}
		T::Alloc(ns);									
	}
	inline void Dealloc() {	if(_ffpP!=NULL){
														delete [] _ffpP; _ffpP = NULL;
														delete [] _ffiP; _ffiP = NULL;
													}
													T::Dealloc();
	}

  static bool HasFFAdjacency()      {   return true; }
  static bool HasFFAdjacencyOcc()   {   return false; }
  static void Name(std::vector<std::string> & name){name.push_back(std::string("PFFAdj"));T::Name(name);}

//private:
  typename T::FacePointer *_ffpP ;    
  char *_ffiP ;    
};

/*----------------------------- PFEADJ ------------------------------*/ 

template <class T> class PFEAdj: public T {
public:
 typedef typename T::EdgeType EdgeType;
	PFEAdj(){_fepP = NULL;  }
	typename T::EdgePointer       &FEp(const int j)        { assert(j>=0 && j<this->VN());  return _fepP[j]; }
  typename T::EdgePointer const  FEp(const int j) const  { assert(j>=0 && j<this->VN());  return _fepP[j]; }
  typename T::EdgePointer const cFEp(const int j) const  { assert(j>=0 && j<this->VN());  return _fepP[j]; }
	
	template <class LeftF>
	void ImportData(const LeftF & leftF){T::ImportData(leftF);}
	inline void Alloc(const int & ns) {
	if( _fepP == NULL){	
		this->SetVN(ns);
		_fepP = new  EdgeType *[this->VN()]; 
		for(int i = 0; i < this->VN(); ++i) {_fepP[i] = 0;}
	}
		T::Alloc(ns);									
	}
	inline void Dealloc() {	if(_fepP!=NULL) {delete [] _fepP; _fepP = NULL;} T::Dealloc();}

	static bool HasFEAdjacency()      {   return true; }
	static bool HasFEAdjacencyOcc()   {   return false; }
	static void Name(std::vector<std::string> & name){name.push_back(std::string("PFEAdj"));T::Name(name);}

//private:
	typename T::EdgePointer *_fepP ;
};


/*----------------------------- PFHADJ ------------------------------*/

template <class T> class PFHAdj: public T {
public:
 typedef typename T::HEdgeType HEdgeType;
 typedef typename T::HEdgePointer HEdgePointer;

	PFHAdj(){_fhP = NULL;  }
	typename T::HEdgePointer       &FHp()        {  return _fhP; }
	typename T::HEdgePointer const cFHp() const  {  return _fhP; }

	template <class LeftF>
	void ImportData(const LeftF & leftF){T::ImportData(leftF);}
	inline void Alloc(const int & ns) {T::Alloc(ns);}
	inline void Dealloc() {	 T::Dealloc();}

	static bool HasFHAdjacency()      {   return true; }
	static bool HasFHAdjacencyOcc()   {   return false; }
	static void Name(std::vector<std::string> & name){name.push_back(std::string("PFHAdj"));T::Name(name);}

//private:
	typename T::HEdgePointer  _fhP ;
};



  } // end namespace face
}// end namespace vcg
#endif
