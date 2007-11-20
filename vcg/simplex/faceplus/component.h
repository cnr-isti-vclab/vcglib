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
Revision 1.17  2007/05/04 16:16:04  ganovelli
added include to texcoor2

Revision 1.16  2007/03/12 15:42:11  tarini
Texture coord name change!  "TCoord" and "Texture" are BAD. "TexCoord" is GOOD.

Revision 1.15  2007/03/12 15:37:19  tarini
Texture coord name change!  "TCoord" and "Texture" are BAD. "TexCoord" is GOOD.

Revision 1.14  2007/02/27 09:32:00  cignoni
Added constructor to the VFadj component to comply to the allocator needs

Revision 1.13  2007/02/12 19:01:23  ganovelli
added Name(std:vector<std::string>& n) that fills n with the names of the attribute of the face type

Revision 1.12  2007/01/11 10:22:39  cignoni
Added intialization of vertexRef to 0.

Revision 1.11  2006/12/06 00:08:57  cignoni
Added FFp1 and FFp2 shortcuts

Revision 1.10  2006/12/04 11:00:02  ganovelli
Cambiate Has*Opt in Has*Occ e aggiunti typedef per la compilazione di Occ

Revision 1.9  2006/11/28 22:34:28  cignoni
Added default constructor with null initialization to adjacency members.
AddFaces and AddVertices NEED to know if the topology is correctly computed to update it.

Revision 1.8  2006/10/07 09:59:42  cignoni
Added missing const to EmptyFF

Revision 1.7  2006/01/09 13:58:55  cignoni
Added Initialization of Color in Vertex and Face Components

Revision 1.6  2005/11/22 15:49:39  cignoni
removed two spurious computenormal

Revision 1.5  2005/11/21 21:44:47  cignoni
Moved ComputeNormal and ComputeNormalizedNormal out of the face class (no more a member function!)

Revision 1.4  2005/11/18 15:44:49  cignoni
Access to constant normal changed from by val to by reference

Revision 1.3  2005/11/16 22:58:17  cignoni
Added IncrementalMark and WedgeTexCoord
Standardized name of flags. It is plural becouse each simplex has many flag.

Revision 1.2  2005/11/12 18:43:14  cignoni
added missing cFFi

Revision 1.1  2005/10/14 15:07:58  cignoni
First Really Working version


****************************************************************************/
#ifndef __VCG_FACE_PLUS_COMPONENT
#define __VCG_FACE_PLUS_COMPONENT

#include <vector>
#include <vcg/space/triangle3.h>
#include <vcg/space/texcoord2.h>
#include <vcg/space/color4.h>

namespace vcg {
  namespace face {
/*
Some naming Rules
All the Components that can be added to a vertex should be defined in the namespace vert:

*/

/*-------------------------- VERTEX ----------------------------------------*/ 
template <class T> class EmptyVertexRef: public T {
public:
 // typedef typename T::VertexType VertexType;
 // typedef typename T::CoordType CoordType;
  inline typename T::VertexType *       & V( const int j ) 	    {	assert(0);		static typename T::VertexType *vp=0; return vp; }
  inline typename T::VertexType * const & V( const int j ) const {	assert(0);		static typename T::VertexType *vp=0; return vp; }
	inline typename T::VertexType * const  cV( const int j ) const {	assert(0);		static typename T::VertexType *vp=0; return vp;	}
	inline       typename T::CoordType & P( const int j ) 	    {	assert(0);		static typename T::CoordType coord(0, 0, 0); return coord;	}
	inline const typename T::CoordType & P( const int j ) const {	assert(0);		static typename T::CoordType coord(0, 0, 0); return coord;	}
	inline const typename T::CoordType &cP( const int j ) const	{	assert(0);		static typename T::CoordType coord(0, 0, 0); return coord;	}
  static bool HasVertexRef()   { return false; }
	static void Name(std::vector<std::string> & name){T::Name(name);}

};
template <class T> class VertexRef: public T {
public:
	VertexRef(){
		v[0]=0;
		v[1]=0;
		v[2]=0;
	}

  inline typename T::VertexType *       & V( const int j ) 	     { assert(j>=0 && j<3); return v[j]; }
  inline typename T::VertexType * const & V( const int j ) const { assert(j>=0 && j<3); return v[j]; }
	inline typename T::VertexType * const  cV( const int j ) const { assert(j>=0 && j<3);	return v[j]; }

	// Shortcut per accedere ai punti delle facce
	inline       typename T::CoordType & P( const int j ) 	    {	assert(j>=0 && j<3);		return v[j]->P();	}
	inline const typename T::CoordType & P( const int j ) const	{	assert(j>=0 && j<3);		return v[j]->cP(); }
	inline const typename T::CoordType &cP( const int j ) const	{	assert(j>=0 && j<3);		return v[j]->cP(); }

	/** Return the pointer to the ((j+1)%3)-th vertex of the face.
		@param j Index of the face vertex.
	 */
	inline       typename T::VertexType *       &  V0( const int j )       { return V(j);}
	inline       typename T::VertexType *       &  V1( const int j )       { return V((j+1)%3);}
	inline       typename T::VertexType *       &  V2( const int j )       { return V((j+2)%3);}
	inline const typename T::VertexType * const &  V0( const int j ) const { return V(j);}
	inline const typename T::VertexType * const &  V1( const int j ) const { return V((j+1)%3);}
	inline const typename T::VertexType * const &  V2( const int j ) const { return V((j+2)%3);}
	inline const typename T::VertexType * const & cV0( const int j ) const { return cV(j);}
	inline const typename T::VertexType * const & cV1( const int j ) const { return cV((j+1)%3);}
	inline const typename T::VertexType * const & cV2( const int j ) const { return cV((j+2)%3);}

	/// Shortcut per accedere ai punti delle facce
	inline       typename T::CoordType &  P0( const int j )       { return V(j)->P();}
	inline       typename T::CoordType &  P1( const int j )       { return V((j+1)%3)->P();}
	inline       typename T::CoordType &  P2( const int j )       { return V((j+2)%3)->P();}
	inline const typename T::CoordType &  P0( const int j ) const { return V(j)->P();}
	inline const typename T::CoordType &  P1( const int j ) const { return V((j+1)%3)->P();}
	inline const typename T::CoordType &  P2( const int j ) const { return V((j+2)%3)->P();}
	inline const typename T::CoordType & cP0( const int j ) const { return cV(j)->P();}
	inline const typename T::CoordType & cP1( const int j ) const { return cV((j+1)%3)->P();}
	inline const typename T::CoordType & cP2( const int j ) const { return cV((j+2)%3)->P();}

	inline       typename T::VertexType *       & UberV( const int j )	      { assert(j>=0 && j<3); return v[j]; }
	inline const typename T::VertexType * const & UberV( const int j ) const	{ assert(j>=0 && j<3);	return v[j];	}
  static bool HasVertexRef()   { return true; }
	static void Name(std::vector<std::string> & name){name.push_back(std::string("VertexRef"));T::Name(name);}


  private:
  typename T::VertexType *v[3];
};



/*-------------------------- NORMAL ----------------------------------------*/ 

template <class T> class EmptyNormal: public T {
public:
  //typedef vcg::Point3s NormalType;
  typedef typename T::VertexType::NormalType NormalType;
  NormalType &N() { static NormalType dummy_normal(0, 0, 0); return dummy_normal; }
  const NormalType &cN() const { static NormalType dummy_normal(0, 0, 0); return dummy_normal; }
  NormalType &WN(int) { static NormalType dummy_normal(0, 0, 0); return dummy_normal; }
  const NormalType cWN(int) const { static NormalType dummy_normal(0, 0, 0); return dummy_normal; }

  static bool HasWedgeNormal()   { return false; }
  static bool HasFaceNormal()   { return false; }
  static bool HasWedgeNormalOcc()   { return false; }
  static bool HasFaceNormalOcc()   { return false; }
//  void ComputeNormal() {assert(0);}
//  void ComputeNormalizedNormal() {assert(0);}
	static void Name(std::vector<std::string> & name){ T::Name(name);}

};
template <class T> class NormalFromVert: public T {
public:
  typedef typename T::VertexType::NormalType NormalType;
  NormalType &N() { return _norm; }
  NormalType &cN() const { return _norm; }
  static bool HasFaceNormal()   { return true; }
//  void ComputeNormal() {	_norm = vcg::Normal<typename T::FaceType>(*(static_cast<typename T::FaceType *>(this))); }
//  void ComputeNormalizedNormal() {	_norm = vcg::NormalizedNormal(*this);}  
	static void Name(std::vector<std::string> & name){name.push_back(std::string("NormalFromVert"));T::Name(name);}

private:
  NormalType _norm;    
};


template <class T>
void ComputeNormal(T &f) {	f.N() = vcg::Normal<T>(f); }

template <class T>
void ComputeNormalizedNormal(T &f) {	f.N() = vcg::NormalizedNormal<T>(f); }

template <class A, class T> class NormalAbs: public T {
public:
  typedef A NormalType;
  NormalType &N() { return _norm; }
  NormalType cN() const { return _norm; }
  static bool HasFaceNormal()   { return true; }
  static void Name(std::vector<std::string> & name){name.push_back(std::string("NormalAbs"));T::Name(name);}

private:
  NormalType _norm;    
};

template <class T> class WedgeNormal: public T {
public:
  typedef typename T::VertexType::NormalType NormalType;
  NormalType &WN(const int j) { return _wnorm[j]; }
  const NormalType cWN(const int j) const { return _wnorm[j]; }
  static bool HasWedgeNormal()   { return true; }
  static void Name(std::vector<std::string> & name){name.push_back(std::string("WedgeNormal"));T::Name(name);}

private:
  NormalType _wnorm[3];    
};


template <class T> class Normal3s: public NormalAbs<vcg::Point3s, T> {
public:static void Name(std::vector<std::string> & name){name.push_back(std::string("Normal3s"));T::Name(name);}
};
template <class T> class Normal3f: public NormalAbs<vcg::Point3f, T> {
public:  static void Name(std::vector<std::string> & name){name.push_back(std::string("Normal3f"));T::Name(name);}
};
template <class T> class Normal3d: public NormalAbs<vcg::Point3d, T> {
public: static void Name(std::vector<std::string> & name){name.push_back(std::string("Normal3d"));T::Name(name);}
};


/*-------------------------- TexCoord ----------------------------------------*/ 

template <class TT> class EmptyWedgeTexCoord: public TT {
public:
	typedef int WedgeTexCoordType;
  typedef vcg::TexCoord2<float,1> TexCoordType;
  TexCoordType &WT(const int) { static TexCoordType dummy_texture; return dummy_texture;}
  TexCoordType const &cWT(const int) const { static TexCoordType dummy_texture; return dummy_texture;}
  static bool HasWedgeTexCoord()   { return false; }
  static bool HasWedgeTexCoordOcc()   { return false; }
  static void Name(std::vector<std::string> & name){TT::Name(name);}

};
template <class A, class TT> class WedgeTexCoord: public TT {
public:
	typedef int WedgeTexCoordType;
  typedef A TexCoordType;
  TexCoordType &WT(const int i) { return _wt[i]; }
  TexCoordType const &cWT(const int i) const { return _wt[i]; }
  static bool HasWedgeTexCoord()   { return true; }
  static void Name(std::vector<std::string> & name){name.push_back(std::string("WedgeTexCoord"));TT::Name(name);}

private:
  TexCoordType _wt[3];    
};

template <class TT> class WedgeTexCoord2s: public WedgeTexCoord<TexCoord2<short,1>, TT> {
public:  static void Name(std::vector<std::string> & name){name.push_back(std::string("WedgeTexCoord2s"));TT::Name(name);}
};
template <class TT> class WedgeTexCoord2f: public WedgeTexCoord<TexCoord2<float,1>, TT> {
public:  static void Name(std::vector<std::string> & name){name.push_back(std::string("WedgeTexCoord2f"));TT::Name(name);}
};
template <class TT> class WedgeTexCoord2d: public WedgeTexCoord<TexCoord2<double,1>, TT> {
public: static void Name(std::vector<std::string> & name){name.push_back(std::string("WedgeTexCoord2d"));TT::Name(name);}
};

/*------------------------- FLAGS -----------------------------------------*/ 
template <class T> class EmptyBitFlags: public T {
public:
	/// Return the vector of Flags(), senza effettuare controlli sui bit
  int &Flags() { static int dummyflags(0); return dummyflags; }
  const int Flags() const { return 0; }
  static bool HasFlags()   { return false; }
  static bool HasFlagsOcc()   { return false; }
	static void Name(std::vector<std::string> & name){T::Name(name);}

};

template <class T> class BitFlags:  public T {
public:
  BitFlags(){_flags=0;}
   int &Flags() {return _flags; }
   const int Flags() const {return _flags; }
  static bool HasFlags()   { return true; }
  static void Name(std::vector<std::string> & name){name.push_back(std::string("BitFlags"));T::Name(name);}


private:
  int  _flags;    
};

/*-------------------------- COLOR ----------------------------------*/ 

template <class T> class EmptyColorQuality: public T {
public:
    typedef float QualityType;
  typedef vcg::Color4b ColorType;
  ColorType &C() { static ColorType dumcolor(vcg::Color4b::White); return dumcolor; }
  ColorType &WC(const int) { static ColorType dumcolor(vcg::Color4b::White); return dumcolor; }
  QualityType &Q() { static QualityType dummyQuality(0); return dummyQuality; }
  static bool HasFaceColor()   { return false; }
  static bool HasWedgeColor()   { return false; }
  static bool HasFaceQuality()   { return false; }
	static bool HasFaceColorOcc() { return false;}
  static void Name(std::vector<std::string> & name){T::Name(name);}

};
template <class A, class T> class Color: public T {
public:
  typedef A ColorType;
  Color():_color(vcg::Color4b::White) {}
  ColorType &C() { return _color; }
  static bool HasFaceColor()   { return true; }
  static void Name(std::vector<std::string> & name){name.push_back(std::string("Color"));T::Name(name);}

private:
  ColorType _color;    
};

template <class A, class T> class WedgeColor: public T {
public:
  typedef A ColorType;
  ColorType &WC(const int i) { return _color[i]; }
  static bool HasFaceColor()   { return true; }
  static void Name(std::vector<std::string> & name){name.push_back(std::string("WedgeColor"));T::Name(name);}

private:
  ColorType _color[3];    
};

template <class T> class Color4b: public Color<vcg::Color4b, T> {
  static void Name(std::vector<std::string> & name){name.push_back(std::string("Color4b"));T::Name(name);}
};

/*-------------------------- Quality  ----------------------------------*/ 

template <class T> class EmptyQuality: public T {
public:
	static void Name(std::vector<std::string> & name){T::Name(name);}

};
template <class A, class T> class Quality: public T {
public:
  typedef A QualityType;
  QualityType &Q() { return _quality; }
  static bool HasFaceQuality()   { return true; }
  static bool HasFaceQualityOcc()   { return true; }
  static void Name(std::vector<std::string> & name){name.push_back(std::string("Quality"));T::Name(name);}
private:
  QualityType _quality;    
};

template <class T> class Qualitys: public Quality<short, T> {
public:  static void Name(std::vector<std::string> & name){name.push_back(std::string("Qualitys"));T::Name(name);}
};
template <class T> class Qualityf: public Quality<float, T> {
public:  static void Name(std::vector<std::string> & name){name.push_back(std::string("Qualityf"));T::Name(name);}
};
template <class T> class Qualityd: public Quality<double, T> {
public:  static void Name(std::vector<std::string> & name){name.push_back(std::string("Qualityd"));T::Name(name);}
};
/*-------------------------- INCREMENTAL MARK  ----------------------------------------*/ 

template <class T> class EmptyMark: public T {
public:
	typedef int MarkType;
  static bool HasMark()   { return false; }
  static bool HasMarkOcc()   { return false; }
  inline void InitIMark()    {  }
  inline int & IMark()       { assert(0); static int tmp=-1; return tmp;}
  inline const int IMark() const {return 0;}
  static void Name(std::vector<std::string> & name){T::Name(name);}

};
template <class T> class Mark: public T {
public:
  static bool HasMark()      { return true; }
  static bool HasMarkOcc()   { return true; }
  inline void InitIMark()    { _imark = 0; }
  inline int & IMark()       { return _imark;}
  inline const int & IMark() const {return _imark;}
  static void Name(std::vector<std::string> & name){name.push_back(std::string("Mark"));T::Name(name);}
    
 private:
	int _imark;
};


/*----------------------------- VFADJ ------------------------------*/ 


template <class T> class EmptyAdj: public T {
public:
	typedef int VFAdjType;
  typename T::FacePointer       &VFp(const int)       { static typename T::FacePointer fp=0; return fp; }
  typename T::FacePointer const cVFp(const int) const { static typename T::FacePointer const fp=0; return fp; }
  typename T::FacePointer       &FFp(const int)       { static typename T::FacePointer fp=0; return fp; }
  typename T::FacePointer const cFFp(const int) const { static typename T::FacePointer const fp=0; return fp; }
  char &VFi(const int j){static char z=0; return z;};
  char &FFi(const int j){static char z=0; return z;};
  static bool HasVFAdjacency()   {   return false; }
  static bool HasFFAdjacency()   {   return false; }
  static bool HasFFAdjacencyOcc()   {   return false; }
  static bool HasVFAdjacencyOcc()   {   return false; }
  static void Name(std::vector<std::string> & name){T::Name(name);}

};

template <class T> class VFAdj: public T {
public:
	VFAdj(){
		_vfp[0]=0;
		_vfp[1]=0;
		_vfp[2]=0;
	}
  typename T::FacePointer       &VFp(const int j)        { assert(j>=0 && j<3);  return _vfp[j]; }
  typename T::FacePointer const  VFp(const int j) const  { assert(j>=0 && j<3);  return _vfp[j]; }
  typename T::FacePointer const cVFp(const int j) const  { assert(j>=0 && j<3);  return _vfp[j]; }
  char &VFi(const int j) {return _vfi[j]; }
  static bool HasVFAdjacency()      {   return true; }
  static bool HasVFAdjacencyOcc()   {   return false; }
  static void Name(std::vector<std::string> & name){name.push_back(std::string("VFAdj"));T::Name(name);}

private:
  typename T::FacePointer _vfp[3] ;    
  char _vfi[3] ;    
};
/*----------------------------- FFADJ ------------------------------*/ 


template <class T> class FFAdj: public T {
public:
	FFAdj(){
		_ffp[0]=0;
		_ffp[1]=0;
		_ffp[2]=0;
	}
  typename T::FacePointer       &FFp(const int j)        { assert(j>=0 && j<3);  return _ffp[j]; }
  typename T::FacePointer const  FFp(const int j) const  { assert(j>=0 && j<3);  return _ffp[j]; }
  typename T::FacePointer const cFFp(const int j) const  { assert(j>=0 && j<3);  return _ffp[j]; }
  char        &FFi(const int j)       { return _ffi[j]; }
  const char &cFFi(const int j) const { return _ffi[j]; }

  typename T::FacePointer        &FFp1( const int j )       { return FFp((j+1)%3);}
	typename T::FacePointer        &FFp2( const int j )       { return FFp((j+2)%3);}
	typename T::FacePointer  const  FFp1( const int j ) const { return FFp((j+1)%3);}
	typename T::FacePointer  const  FFp2( const int j ) const { return FFp((j+2)%3);}

  static bool HasFFAdjacency()      {   return true; }
  static bool HasFFAdjacencyOcc()   {   return false; }
  static void Name(std::vector<std::string> & name){name.push_back(std::string("FFAdj"));T::Name(name);}

private:
  typename T::FacePointer _ffp[3] ;    
  char _ffi[3] ;    
};

  } // end namespace vert
}// end namespace vcg
#endif
