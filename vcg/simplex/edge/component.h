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

#ifndef __VCG_EDGE_PLUS_COMPONENT
#define __VCG_EDGE_PLUS_COMPONENT
#include <vector>
#include <string>
//#include <vcg/space/point3.h>
//#include <vcg/space/texcoord2.h>
#include <vcg/space/color4.h>

namespace vcg {
  namespace edge {
/*
Some naming Rules
All the Components that can be added to a vertex should be defined in the namespace edge:

*/

/*-------------------------- VERTEX ----------------------------------------*/ 
template <class T> class EmptyVertexRef: public T {
public:
 // typedef typename T::VertexType VertexType;
 // typedef typename T::CoordType CoordType;
  inline typename T::VertexType *       & V( const int j ) 	    {	assert(0);		static typename T::VertexType *vp=0; return vp; }
  inline typename T::VertexType * const & V( const int j ) const {	assert(0);		static typename T::VertexType *vp=0; return vp; }
        inline typename T::VertexType *  cV( const int j ) const {	assert(0);		static typename T::VertexType *vp=0; return vp;	}
	inline       typename T::CoordType & P( const int j ) 	    {	assert(0);		static typename T::CoordType coord(0, 0, 0); return coord;	}
	inline const typename T::CoordType & P( const int j ) const {	assert(0);		static typename T::CoordType coord(0, 0, 0); return coord;	}
	inline const typename T::CoordType &cP( const int j ) const	{	assert(0);		static typename T::CoordType coord(0, 0, 0); return coord;	}
	template <class LeftF>
	void ImportData(const LeftF & leftF) {T::ImportData(leftF);}
  static bool HasVertexRef()   { return false; }
	static void Name(std::vector<std::string> & name){T::Name(name);}

};
template <class T> class VertexRef: public T {
public:
	VertexRef(){
		v[0]=0;
		v[1]=0;
	}

  inline typename T::VertexType *       & V( const int j ) 	     { assert(j>=0 && j<2); return v[j]; }
  inline typename T::VertexType * const & V( const int j ) const { assert(j>=0 && j<2); return v[j]; }
        inline typename T::VertexType *  cV( const int j ) const { assert(j>=0 && j<2);	return v[j]; }

	// Shortcut per accedere ai punti delle facce
	inline       typename T::CoordType & P( const int j ) 	    {	assert(j>=0 && j<2);		return v[j]->P();	}
	inline const typename T::CoordType &cP( const int j ) const	{	assert(j>=0 && j<2);		return v[j]->cP(); }

	/** Return the pointer to the ((j+1)%3)-th vertex of the face.
		@param j Index of the face vertex.
	 */
	inline       typename T::VertexType *       &  V0( const int j )       { return V(j);}
	inline       typename T::VertexType *       &  V1( const int j )       { return V((j+1)%2);}
	inline const typename T::VertexType * const &  V0( const int j ) const { return V(j);}
	inline const typename T::VertexType * const &  V1( const int j ) const { return V((j+1)%2);}
	inline const typename T::VertexType * const & cV0( const int j ) const { return cV(j);}
	inline const typename T::VertexType * const & cV1( const int j ) const { return cV((j+1)%2);}
	
	/// Shortcut per accedere ai punti delle facce
	inline       typename T::CoordType &  P0( const int j )       { return V(j)->P();}
	inline       typename T::CoordType &  P1( const int j )       { return V((j+1)%2)->P();}
	inline const typename T::CoordType &  P0( const int j ) const { return V(j)->P();}
	inline const typename T::CoordType &  P1( const int j ) const { return V((j+1)%2)->P();}
	inline const typename T::CoordType & cP0( const int j ) const { return cV(j)->P();}
	inline const typename T::CoordType & cP1( const int j ) const { return cV((j+1)%2)->P();}

	inline       typename T::VertexType *       & UberV( const int j )	      { assert(j>=0 && j<2); return v[j]; }
	inline const typename T::VertexType * const & UberV( const int j ) const	{ assert(j>=0 && j<2);	return v[j];	}

	template <class LeftF>
	void ImportData(const LeftF & leftF){ T::ImportData(leftF);}

  static bool HasVertexRef()   { return true; }
	static void Name(std::vector<std::string> & name){name.push_back(std::string("VertexRef"));T::Name(name);}


  private:
  typename T::VertexType *v[2];
};



/*-------------------------- INCREMENTAL MARK  ----------------------------------------*/ 

template <class T> class EmptyMark: public T {
public:
  static bool HasMark()   { return false; }
  static bool HasMarkOcc()   { return false; }
  inline void InitIMark()    {  }
  inline int & IMark()       { assert(0); static int tmp=-1; return tmp;}
  inline const int & IMark() const {return 0;}
	template < class LeftV>
	void ImportData(const LeftV  & left ) { T::ImportData( left); }
	static void Name(std::vector<std::string> & name){T::Name(name);}

};
template <class T> class Mark: public T {
public:
  static bool HasMark()      { return true; }
  static bool HasMarkOcc()   { return true; }
  inline void InitIMark()    { _imark = 0; }
  inline int & IMark()       { return _imark;}
  inline const int & IMark() const {return _imark;}
	template < class LeftV>
	void ImportData(const LeftV  & left ) { IMark() = left.IMark(); T::ImportData( left); }
	static void Name(std::vector<std::string> & name){name.push_back(std::string("Mark"));T::Name(name);}
    
 private:
	int _imark;
};

/*------------------------- FLAGS -----------------------------------------*/ 
template <class T> class EmptyBitFlags: public T {
public:
	typedef int FlagType;
  /// Return the vector of Flags(), senza effettuare controlli sui bit
  int &Flags() { static int dummyflags(0);  assert(0); return dummyflags; }
  int Flags() const { return 0; }
	template < class LeftV>
	void ImportData(const LeftV  & left ) { T::ImportData( left); }
  static bool HasFlags()   { return false; }
	static void Name(std::vector<std::string> & name){T::Name(name);}

};

template <class T> class BitFlags:  public T {
public:
	BitFlags(){_flags=0;}
  typedef int FlagType;
  int &Flags() {return _flags; }
  int Flags() const {return _flags; }
	template < class LeftV>
	void ImportData(const LeftV  & left ) { Flags() = left.Flags(); T::ImportData( left); }
  static bool HasFlags()   { return true; }
	static void Name(std::vector<std::string> & name){name.push_back(std::string("BitFlags"));T::Name(name);}

private:
  int  _flags;    
};

/*-------------------------- EMPTY COLOR & QUALITY ----------------------------------*/ 

template <class T> class EmptyColorQuality: public T {
public:
  typedef float QualityType;
  QualityType &Q() { static QualityType dummyQuality(0);  assert(0); return dummyQuality; }
	static bool HasQuality()   { return false; }

  typedef vcg::Color4b ColorType;
  ColorType &C() { static ColorType dumcolor(vcg::Color4b::White); assert(0); return dumcolor; }
	template < class LeftV>
	void ImportData(const LeftV  & left ) { T::ImportData( left); }
  static bool HasColor()   { return false; }
	static void Name(std::vector<std::string> & name){T::Name(name);}
};

/*-------------------------- Color  ----------------------------------*/ 

template <class A, class T> class Color: public T {
public:
  Color():_color(vcg::Color4b::White) {}
  typedef A ColorType;
  ColorType &C() { return _color; }
  const ColorType &C() const { return _color; }
  const ColorType &cC() const { return _color; }
	template < class LeftV>
	void ImportData(const LeftV  & left ) { C() = left.cC(); T::ImportData( left); }
  static bool HasColor()   { return true; }
	static void Name(std::vector<std::string> & name){name.push_back(std::string("Color"));T::Name(name);}

private:
  ColorType _color;    
};

template <class TT> class Color4b: public edge::Color<vcg::Color4b, TT> {
	public: static void Name(std::vector<std::string> & name){name.push_back(std::string("Color4b"));TT::Name(name);}
};

/*-------------------------- Quality  ----------------------------------*/ 

template <class A, class TT> class Quality: public TT {
public:
  typedef A QualityType;
  QualityType &Q() { return _quality; }
  const QualityType & cQ() const {return _quality; }
	template < class LeftV>
	void ImportData(const LeftV  & left ) { Q() = left.cQ(); TT::ImportData( left); }
  static bool HasQuality()   { return true; }
	static void Name(std::vector<std::string> & name){name.push_back(std::string("Quality"));TT::Name(name);}

private:
  QualityType _quality;    
};

template <class TT> class Qualitys: public Quality<short, TT> {
public: static void Name(std::vector<std::string> & name){name.push_back(std::string("Qualitys"));TT::Name(name);}
};
template <class TT> class Qualityf: public Quality<float, TT> {
public: static void Name(std::vector<std::string> & name){name.push_back(std::string("Qualityf"));TT::Name(name);}
};
template <class TT> class Qualityd: public Quality<double, TT> {
public: static void Name(std::vector<std::string> & name){name.push_back(std::string("Qualityd"));TT::Name(name);}
};

/*----------------------------- EVADJ ------------------------------*/ 
template <class T> class EmptyEVAdj: public T {
public:
	typename T::VertexPointer &EVp(const int &) { static typename T::VertexPointer ep=0;  assert(0); return ep; }
	 typename T::VertexPointer const cEVp(const int &) const  { static typename T::VertexPointer ep=0;  assert(0); return ep; }
	 typename T::VertexPointer &V(const int &i) { return EVp(i);}
		typename T::VertexPointer const cV(const int &i) const  { return cEVp(i); }
	template < class LeftV>
		void ImportData(const LeftV  & left ) { T::ImportData( left); }
  static bool HasEVAdjacency()   {   return false; }
  static bool HasEVAdjacencyOcc()   {   return false; }
	static void Name(std::vector<std::string> & name){ T::Name(name);}
};

template <class T> class EVAdj: public T {
public:
  EVAdj(){_vp[0]= _vp[1] =0;}
	typename T::VertexPointer			&	EVp(const int & i) {return _vp[i]; }
	typename T::VertexPointer const cEVp(const int & i) const {return _vp[i]; }
	typename T::VertexPointer &V(const int &i) { return EVp(i);}
	typename T::VertexPointer const cV(const int &i) const  { return cEVp(i); }

	template < class LeftV>
	void ImportData(const LeftV  & left ) {  T::ImportData( left); }
  static bool HasEVAdjacency()   {   return true; }
  static bool HasEVAdjacencyOcc()   {   return true; }
	static void Name(std::vector<std::string> & name){name.push_back(std::string("EVAdj"));T::Name(name);}
	
private:
   typename T::VertexPointer	 _vp[2] ;    
};


/*----------------------------- EEADJ ------------------------------*/ 
template <class T> class EmptyEEAdj: public T {
public:
  typename T::EdgePointer &EEp(const int &  ) { static typename T::EdgePointer ep=0;  assert(0); return ep; }
	const typename T::EdgePointer cEEp(const int & ) const { static typename T::EdgePointer ep=0;  assert(0); return ep; }
	int &EEi(const int &){static int z=0; assert(0); return z;};
	const int &cEEi(const int &) const {static int z=0; assert(0); return z;};
	template < class LeftV>
		void ImportData(const LeftV  & left ) { T::ImportData( left); }
  static bool HasEEAdjacency()   {   return false; }
  static bool HasEEAdjacencyOcc()   {   return false; }
	static void Name(std::vector<std::string> & name){ T::Name(name);}
};

template <class T> class EEAdj: public T {
public:
  EEAdj(){_ep=0;}
  typename T::EdgePointer &EEp(const int & i) {return _ep[i]; }
  typename T::EdgePointer cEEp(const int & i) {return _ep[i]; }
	int &EEi(const int & i){ return _zp[i];};
	const int &cEEi(const int &i ){return _zp[i];};

	template < class LeftV>
	void ImportData(const LeftV  & left ) {  T::ImportData( left); }
  static bool HasEEAdjacency()   {   return true; }
  static bool HasEEAdjacencyOcc()   {   return true; }
	static void Name(std::vector<std::string> & name){name.push_back(std::string("EEAdj"));T::Name(name);}

private:
  typename T::EdgePointer _ep[2] ;    
  int _zp[2] ;    
};


/*----------------------------- EHADJ ------------------------------*/ 
template <class T> class EmptyEHAdj: public T {
public:
  typename T::HEdgePointer &EHp(  ) { static typename T::HEdgePointer hp=0;  assert(0); return hp; }
	const typename T::HEdgePointer cEHp(  ) const { static typename T::HEdgePointer hp=0;  assert(0); return hp; }
  
	template < class LeftV>
		void ImportData(const LeftV  & left ) { T::ImportData( left); }
  static bool HasEHAdjacency()   {   return false; }
  static bool HasEHAdjacencyOcc()   {   return false; }
	static void Name(std::vector<std::string> & name){ T::Name(name);}
};

template <class T> class EHAdj: public T {
public:
  EHAdj(){_hp=0;}
  typename T::HEdgePointer &EHp( ) {return _hp ; }
	const typename T::HEdgePointer cEHp( ) const {return _hp ; }
 
	template < class LeftV>
	void ImportData(const LeftV  & left ) { T::ImportData( left); }
  static bool HasEHAdjacency()   {   return true; }
  static bool HasEHAdjacencyOcc()   {   return true; }
	static void Name(std::vector<std::string> & name){name.push_back(std::string("EHAdj"));T::Name(name);}

private:
  typename T::HEdgePointer _hp ;
};

/*----------------------------- ETADJ ------------------------------*/ 


template <class T> class EmptyETAdj: public T {
public:
	typename T::TetraPointer &ETp() { static typename T::TetraPointer tp = 0;  assert(0); return tp; }
	typename T::TetraPointer cETp() { static typename T::TetraPointer tp = 0;  assert(0); return tp; }
	int &VTi() { static int z = 0; return z; };
	static bool HasETAdjacency() { return false; }
	static bool HasETAdjacencyOcc() { return false; }
	static void Name( std::vector< std::string > & name ) { T::Name(name); }
};

template <class T> class ETAdj: public T {
public:
	ETAdj() { _tp = 0; }
	typename T::TetraPointer &ETp() { return _tp; }
	typename T::TetraPointer cETp() { return _tp; }
	int &ETi() {return _zp; }
	static bool HasETAdjacency() { return true; }
	static bool HasETAdjacencyOcc()   {   return true; }
	static void Name( std::vector< std::string > & name ) { name.push_back( std::string("ETAdj") ); T::Name(name); }

private:
	typename T::TetraPointer _tp ;    
	int _zp ;    
};



/*----------------------------- EFADJ ------------------------------*/ 

template <class T> class EmptyEFAdj: public T {
public:
  typename T::FacePointer &EFp() { static typename T::FacePointer fp=0;  assert(0); return fp; }
	const typename T::FacePointer cEFp() const  { static typename T::FacePointer fp=0;  assert(0); return fp; }
	int &EFi()   {static int z=0; return z;};
	const int &cEFi() const {static int z=0; return z;};
	template < class LeftV>
	void ImportData(const LeftV  & left ) { T::ImportData( left); }
  static bool HasEFAdjacency()   {   return false; }
  static bool HasEFAdjacencyOcc()   {   return false; }
	static void Name(std::vector<std::string> & name){ T::Name(name);}
};

template <class T> class EFAdj: public T {
public:
  EFAdj(){_fp=0;}
  typename T::FacePointer &EFp() {return _fp; }
	const typename T::FacePointer cEFp() const {return _fp; }
	int &EFi()   {static int z=0; return z;};
	const int &cEFi() const  {return _zp; }
	template < class LeftV>
	void ImportData(const LeftV  & left ) {  T::ImportData( left); }
  static bool HasEFAdjacency()   {   return true; }
  static bool HasEFAdjacencyOcc()   {   return true; }
	static void Name(std::vector<std::string> & name){name.push_back(std::string("EFAdj"));T::Name(name);}

private:
  typename T::FacePointer _fp ;    
  int _zp ;    
};


  } // end namespace edge
}// end namespace vcg
#endif
