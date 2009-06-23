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
//#include <vector>
//#include <string>
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
	void ImportLocal(const LeftF & leftF) {T::ImportLocal(leftF);}
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
	void ImportLocal(const LeftF & leftF){ V(0) = NULL; V(1) = NULL; V(2) = NULL; T::ImportLocal(leftF);}

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
	void ImportLocal(const LeftV  & left ) { T::ImportLocal( left); }
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
	void ImportLocal(const LeftV  & left ) { IMark() = left.IMark(); T::ImportLocal( left); }
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
	void ImportLocal(const LeftV  & left ) { T::ImportLocal( left); }
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
	void ImportLocal(const LeftV  & left ) { Flags() = left.Flags(); T::ImportLocal( left); }
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
	void ImportLocal(const LeftV  & left ) { T::ImportLocal( left); }
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
	void ImportLocal(const LeftV  & left ) { C() = left.cC(); T::ImportLocal( left); }
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
	void ImportLocal(const LeftV  & left ) { Q() = left.cQ(); TT::ImportLocal( left); }
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
  typename T::VertexPointer &V(const int &) { static typename T::VertexPointer ep=0;  assert(0); return ep; }
  typename T::VertexPointer cV(const int &) { static typename T::VertexPointer ep=0;  assert(0); return ep; }
  int &EVi(){static int z=0; return z;};
	template < class LeftV>
		void ImportLocal(const LeftV  & left ) { T::ImportLocal( left); }
  static bool HasEVAdjacency()   {   return false; }
  static bool HasEVAdjacencyOcc()   {   return false; }
	static void Name(std::vector<std::string> & name){ T::Name(name);}
};

template <class T> class EVAdj: public T {
public:
  EVAdj(){_vp[0]= _vp[1] =0;}
  typename T::VertexPointer			&	V(const int & i) {return _vp[i]; }
  const typename T::VertexPointer cV(const int & i) const {return _vp[i]; }
	template < class LeftV>
	void ImportLocal(const LeftV  & left ) { V() = NULL; T::ImportLocal( left); }
  static bool HasEVAdjacency()   {   return true; }
  static bool HasEVAdjacencyOcc()   {   return true; }
	static void Name(std::vector<std::string> & name){name.push_back(std::string("EVAdj"));T::Name(name);}
	
private:
   typename T::VertexPointer	 _vp[2] ;    
};

/*----------------------------- HEVADJ ------------------------------*/ 
template <class T> class EmptyHEVAdj: public T {
public:
  typename T::VertexPointer &HEVp() { static typename T::VertexPointer ep=0;  assert(0); return ep; }
  typename T::VertexPointer cHEVp() { static typename T::VertexPointer ep=0;  assert(0); return ep; }
  int &EVi(){static int z=0; return z;};
	template < class LeftV>
		void ImportLocal(const LeftV  & left ) { T::ImportLocal( left); }
  static bool HasHEVAdjacency()   {   return false; }
  static bool HasHEVAdjacencyOcc()   {   return false; }
	static void Name(std::vector<std::string> & name){ T::Name(name);}
};

template <class T> class HEVAdj: public T {
public:
  HEVAdj(){_vp =0;}
  typename T::VertexPointer			&	HEVp() {return _vp ; }
  const typename T::VertexPointer cHEVp() const {return _vp ; }
	template < class LeftV>
	void ImportLocal(const LeftV  & left ) { this->V() = NULL; T::ImportLocal( left); }
  static bool HasHEVAdjacency()   {   return true; }
  static bool HasHEVAdjacencyOcc()   {   return true; }
	static void Name(std::vector<std::string> & name){name.push_back(std::string("HEVAdj"));T::Name(name);}
	
private:
   typename T::VertexPointer	 _vp ;    
};

/*----------------------------- EEADJ ------------------------------*/ 
template <class T> class EmptyEEAdj: public T {
public:
  typename T::EdgePointer &EEp(const int &  ) { static typename T::EdgePointer ep=0;  assert(0); return ep; }
  typename T::EdgePointer cEEp(const int & ) { static typename T::EdgePointer ep=0;  assert(0); return ep; }
  int &EEi(){static int z=0; return z;};
	template < class LeftV>
		void ImportLocal(const LeftV  & left ) { T::ImportLocal( left); }
  static bool HasEEAdjacency()   {   return false; }
  static bool HasEEAdjacencyOcc()   {   return false; }
	static void Name(std::vector<std::string> & name){ T::Name(name);}
};

template <class T> class EEAdj: public T {
public:
  EEAdj(){_ep=0;}
  typename T::EdgePointer &EEp(const int & i) {return _ep[i]; }
  typename T::EdgePointer cEEp(const int & i) {return _ep[i]; }
  int &EEi(const int & i) {return _zp[i]; }
	template < class LeftV>
	void ImportLocal(const LeftV  & left ) { EEp() = NULL; T::ImportLocal( left); }
  static bool HasEEAdjacency()   {   return true; }
  static bool HasEEAdjacencyOcc()   {   return true; }
	static void Name(std::vector<std::string> & name){name.push_back(std::string("EEAdj"));T::Name(name);}

private:
  typename T::EdgePointer _ep[2] ;    
  int _zp[2] ;    
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



/*----------------------------- HENextADJ ------------------------------*/ 
template <class T> class EmptyHENextAdj: public T {
public:
  typename T::EdgePointer &HENp( ) { static typename T::EdgePointer ep=0;  assert(0); return ep; }
  typename T::EdgePointer cHEp( ) { static typename T::EdgePointer ep=0;  assert(0); return ep; }
	template < class LeftV>
		void ImportLocal(const LeftV  & left ) { T::ImportLocal( left); }
  static bool HasHENextAdjacency()   {   return false; }
  static bool HasHENextAdjacencyOcc()   {   return false; }
	static void Name(std::vector<std::string> & name){ T::Name(name);}
};

template <class T> class HENextAdj: public T {
public:
  HENextAdj(){_nep=0;}
  typename T::EdgePointer &HENp() {return _nep; }
  typename T::EdgePointer cHENp() {return _nep; }
	template < class LeftV>
	void ImportLocal(const LeftV  & left ) { this->EEp() = NULL; T::ImportLocal( left); }
  static bool HasHENextAdjacency()   {   return true; }
	static void Name(std::vector<std::string> & name){name.push_back(std::string("HENextAdj"));T::Name(name);}
	
private:
  typename T::EdgePointer _nep ;    
};

/*----------------------------- HEOppADJ ------------------------------*/ 
template <class T> class EmptyHEOppAdj: public T {
public:
  typename T::EdgePointer &HEOp(const int & i ) { static typename T::EdgePointer ep=0;  assert(0); return ep; }
  typename T::EdgePointer cHOp(const int & i) { static typename T::EdgePointer ep=0;  assert(0); return ep; }
  int &EEi(){static int z=0; return z;};
	template < class LeftV>
		void ImportLocal(const LeftV  & left ) { T::ImportLocal( left); }
  static bool HasHEOppAdjacency()   {   return false; }
  static bool HasHEOpptAdjacencyOcc()   {   return false; }
	static void Name(std::vector<std::string> & name){ T::Name(name);}
};

template <class T> class HEOppAdj: public T {
public:
  HEOppAdj(){_oep=0;}
  typename T::EdgePointer &HEOp() {return _oep; }
  typename T::EdgePointer cHEOp() {return _oep; }
	template < class LeftV>
	void ImportLocal(const LeftV  & left ) { this->EEp() = NULL; T::ImportLocal( left); }
  static bool HasHEOppAdjacency()   {   return true; }
	static void Name(std::vector<std::string> & name){name.push_back(std::string("HEOpptAdj"));T::Name(name);}
	
private:
  typename T::EdgePointer _oep ;    
 
};
/*----------------------------- HEPrevADJ ------------------------------*/ 
template <class T> class EmptyHEPrevAdj: public T {
public:
  typename T::EdgePointer &HEPp() { static typename T::EdgePointer ep=0;  assert(0); return ep; }
  typename T::EdgePointer cHPp() { static typename T::EdgePointer ep=0;  assert(0); return ep; }
  int &EEi(){static int z=0; return z;};
	template < class LeftV>
		void ImportLocal(const LeftV  & left ) { T::ImportLocal( left); }
  static bool HasHEPrevAdjacency()   {   return false; }
  static bool HasHEPrevAdjacencyOcc()   {   return false; }
  static void Name(std::vector<std::string> & name){ T::Name(name);}
};

template <class T> class HEPrevAdj: public T {
public:
  HEPrevAdj(){_pep=0;}
  typename T::EdgePointer &HEPp() {return _pep; }
  typename T::EdgePointer cHEPp() {return _pep; }
  int &EEi(const int & i) {return this->_nei[i]; }
	template < class LeftV>
	void ImportLocal(const LeftV  & left ) { this->EEp() = NULL; T::ImportLocal( left); }
  static bool HasHEPrevAdjacency()   {   return true; }
	static void Name(std::vector<std::string> & name){name.push_back(std::string("HEPrevAdj"));T::Name(name);}
	
private:
  typename T::EdgePointer _pep ;    
};
/*----------------------------- EFADJ ------------------------------*/ 

template <class T> class EmptyEFAdj: public T {
public:
  typename T::FacePointer &EFp() { static typename T::FacePointer fp=0;  assert(0); return fp; }
  typename T::FacePointer cEFp() { static typename T::FacePointer fp=0;  assert(0); return fp; }
  int &EFi(){static int z=0; return z;};
	template < class LeftV>
	void ImportLocal(const LeftV  & left ) { T::ImportLocal( left); }
  static bool HasEFAdjacency()   {   return false; }
  static bool HasEFAdjacencyOcc()   {   return false; }
	static void Name(std::vector<std::string> & name){ T::Name(name);}
};

template <class T> class EFAdj: public T {
public:
  EFAdj(){_fp=0;}
  typename T::FacePointer &EFp() {return _fp; }
  typename T::FacePointer cEFp() {return _fp; }
  int &EFi() {return _zp; }
	template < class LeftV>
	void ImportLocal(const LeftV  & left ) { this->EFp() = NULL; T::ImportLocal( left); }
  static bool HasEFAdjacency()   {   return true; }
  static bool HasEFAdjacencyOcc()   {   return true; }
	static void Name(std::vector<std::string> & name){name.push_back(std::string("EFAdj"));T::Name(name);}

private:
  typename T::FacePointer _fp ;    
  int _zp ;    
};


/*----------------------------- EFADJ ------------------------------*/ 
/**
 HEdgeData keep all the data for the half edge
*/
template <class T> 
class EmptyHEdgeData : public	EmptyEFAdj<		// pointer to the face
							EmptyHEOppAdj <		// pointer to the opposite half edge
							EmptyHENextAdj <	// pointer to the next half edge along the face
							EmptyHEVAdj <		// pointer to the vertex
							EmptyHEPrevAdj<
							T > > > > >{};


template <class T> 
class HEdgeData : public	EFAdj<			// pointer to the face
							HEOppAdj <		// pointer to the opposite half edge
							HENextAdj <		// pointer to the next half edge along the face
							HEVAdj <		// pointer to the vertex
							T > > > >{

	// functions to make the half edge user confortable
	typename T::VertexPointer & Vertex()					{ return this->HEVp();}
	const typename T::VertexPointer &  cVertex()	const	{ return this->cHEVp();}
	typename T::EdgePointer Opposite()						{ return &this->HEOp();}
	const typename T::EdgePointer & cOpposite()		const	{ return this->cHEOp();}
	typename T::EdgePointer & Next()						{ return this->HENp();}
	const typename T::EdgePointer &  cNext()			const 	{ return this->HENp();}

};

  } // end namespace edge
}// end namespace vcg
#endif
