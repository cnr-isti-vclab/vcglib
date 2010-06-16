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
Revision 1.21  2008/02/04 21:26:45  ganovelli
added ImportData which imports all local attributes into vertexplus and faceplus.
A local attribute is everything (N(), C(), Q()....) except pointers to other simplices
(i.e. FFAdj, VFAdj, VertexRef) which are set to NULL.
Added some function for const attributes

Revision 1.20  2008/01/28 08:42:51  cignoni
added assert when writing on empty data members

Revision 1.19  2008/01/19 17:49:05  ganovelli
missing const  cVF added

Revision 1.18  2007/11/20 09:43:53  ganovelli
added missing include to color4

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
  inline typename T::VertexType *       & V( const int ) 	    {	assert(0);		static typename T::VertexType *vp=0; return vp; }
  inline typename T::VertexType * const & V( const int ) const {	assert(0);		static typename T::VertexType *vp=0; return vp; }
  inline typename T::VertexType * cV( const int ) const {	assert(0);		static typename T::VertexType *vp=0; return vp;	}

	inline typename T::VertexType *       & FVp( const int i )				{	return this->V(i); }
	inline typename T::VertexType * const & FVp( const int i) const		{	return this->V(i); }
	inline typename T::VertexType * cFVp( const int i ) const					{	return this->cV(i); }

  inline typename T::CoordType & P( const int ) 	    {	assert(0);		static typename T::CoordType coord(0, 0, 0); return coord;	}
  inline const typename T::CoordType & P( const int ) const {	assert(0);		static typename T::CoordType coord(0, 0, 0); return coord;	}
  inline const typename T::CoordType &cP( const int ) const	{	assert(0);		static typename T::CoordType coord(0, 0, 0); return coord;	}
	template <class RightF>
	void ImportData(const RightF & rightF) {T::ImportData(rightF);}
	inline void Alloc(const int & ns){T::Alloc(ns);}
	inline void Dealloc(){T::Dealloc();}
  static bool HasVertexRef()   { return false; }
	static bool HasFVAdjacency()   { return false; }
	static void Name(std::vector<std::string> & name){T::Name(name);}

};
template <class T> class VertexRef: public T {
public:
	VertexRef(){
		v[0]=0;
		v[1]=0;
		v[2]=0;
	}

  typedef typename T::VertexType::CoordType CoordType;
  typedef typename T::VertexType::ScalarType ScalarType;

  inline typename T::VertexType *       & V( const int j ) 	     { assert(j>=0 && j<3); return v[j]; }
  inline typename T::VertexType * const & V( const int j ) const { assert(j>=0 && j<3); return v[j]; }
	inline typename T::VertexType *  cV( const int j ) const { assert(j>=0 && j<3);	return v[j]; }

	// Shortcut per accedere ai punti delle facce
  inline       CoordType & P( const int j ) 	    {	assert(j>=0 && j<3);		return v[j]->P();	}
  inline const CoordType & P( const int j ) const	{	assert(j>=0 && j<3);		return v[j]->cP(); }
  inline const CoordType &cP( const int j ) const	{	assert(j>=0 && j<3);		return v[j]->cP(); }

	/** Return the pointer to the ((j+1)%3)-th vertex of the face.
		@param j Index of the face vertex.
	 */
  inline typename T::VertexType *      &  V0( const int j )       { return V(j);}
  inline typename T::VertexType *      &  V1( const int j )       { return V((j+1)%3);}
  inline typename T::VertexType *      &  V2( const int j )       { return V((j+2)%3);}
  inline typename T::VertexType * const   V0( const int j ) const { return V(j);}
  inline typename T::VertexType * const   V1( const int j ) const { return V((j+1)%3);}
  inline typename T::VertexType * const   V2( const int j ) const { return V((j+2)%3);}
  inline typename T::VertexType *  cV0( const int j ) const { return cV(j);}
  inline typename T::VertexType *  cV1( const int j ) const { return cV((j+1)%3);}
  inline typename T::VertexType *  cV2( const int j ) const { return cV((j+2)%3);}

	/// Shortcut per accedere ai punti delle facce
  inline       CoordType &  P0( const int j )       { return V(j)->P();}
  inline       CoordType &  P1( const int j )       { return V((j+1)%3)->P();}
  inline       CoordType &  P2( const int j )       { return V((j+2)%3)->P();}
  inline const CoordType &  P0( const int j ) const { return V(j)->P();}
  inline const CoordType &  P1( const int j ) const { return V((j+1)%3)->P();}
  inline const CoordType &  P2( const int j ) const { return V((j+2)%3)->P();}
  inline const CoordType & cP0( const int j ) const { return cV(j)->P();}
  inline const CoordType & cP1( const int j ) const { return cV((j+1)%3)->P();}
  inline const CoordType & cP2( const int j ) const { return cV((j+2)%3)->P();}

	inline       typename T::VertexType *       & UberV( const int j )	      { assert(j>=0 && j<3); return v[j]; }
	inline const typename T::VertexType * const & UberV( const int j ) const	{ assert(j>=0 && j<3);	return v[j];	}

    // Small comment about the fact that the pointers are zero filled.
    // The importLocal is meant for copyng stuff between very different meshes, so copying the pointers would be meaningless.
		// if you are using ImportData for copying internally simplex you have to set up all the pointers by hand.
	template <class RightF>
	void ImportData(const RightF & rightF){  T::ImportData(rightF);}
	inline void Alloc(const int & ns){T::Alloc(ns);}
	inline void Dealloc(){T::Dealloc();}

  static bool HasVertexRef()   { return true; }
	static bool HasFVAdjacency()   { return true; }

	static void Name(std::vector<std::string> & name){name.push_back(std::string("VertexRef"));T::Name(name);}


  private:
  typename T::VertexType *v[3];
};



/*-------------------------- NORMAL ----------------------------------------*/ 

template <class T> class EmptyNormal: public T {
public:
  //typedef vcg::Point3s NormalType;
  typedef typename T::VertexType::NormalType NormalType;
  NormalType &N() { static NormalType dummy_normal(0, 0, 0);  assert(0); return dummy_normal; }
  const NormalType &cN() const { static NormalType dummy_normal(0, 0, 0); return dummy_normal; }
  NormalType &WN(int) { static NormalType dummy_normal(0, 0, 0);  assert(0); return dummy_normal; }
  const NormalType cWN(int) const { static NormalType dummy_normal(0, 0, 0); return dummy_normal; }

	template <class RightF>
	void ImportData(const RightF & rightF){ T::ImportData(rightF);}
  static bool HasWedgeNormal()   { return false; }
  static bool HasFaceNormal()   { return false; }
  static bool HasWedgeNormalOcc()   { return false; }
  static bool HasFaceNormalOcc()   { return false; }
  static bool HasWedgeNormalOcf()   { return false; }
  static bool HasFaceNormalOcf()   { return false; }
//  void ComputeNormal() {assert(0);}
//  void ComputeNormalizedNormal() {assert(0);}
	static void Name(std::vector<std::string> & name){ T::Name(name);}

};
template <class T> class NormalFromVert: public T {
public:
  typedef typename T::VertexType::NormalType NormalType;
  NormalType &N() { return _norm; }
  NormalType &cN() const { return _norm; }
	template <class RightF>
	void ImportData(const RightF & rightF){ N() = rightF.cN(); T::ImportData(rightF);}
 	inline void Alloc(const int & ns){T::Alloc(ns);}
	inline void Dealloc(){T::Dealloc();}
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
	template <class RightF>
	void ImportData(const RightF & rightF)
	{		
		N().Import(rightF.cN());
		T::ImportData( rightF);
	}

 	inline void Alloc(const int & ns){T::Alloc(ns);}
	inline void Dealloc(){T::Dealloc();}
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
	template <class RightF>
	void ImportData(const RightF & rightF){ for (int i=0; i<3; ++i) { WN(i) = rightF.cWN(i); } T::ImportData(rightF);}
 	inline void Alloc(const int & ns){T::Alloc(ns);}
	inline void Dealloc(){T::Dealloc();}
	static bool HasWedgeNormal()   { return true; }
  static void Name(std::vector<std::string> & name){name.push_back(std::string("WedgeNormal"));T::Name(name);}

private:
  NormalType _wnorm[3];    
};

template <class A, class T> class WedgeRealNormal: public T {
public:
  typedef A NormalType;
  NormalType &WN(const int i) { return _wn[i]; }
  NormalType const &cWN(const int i) const { return _wn[i]; }
	template <class RightF>
	void ImportData(const RightF & rightF){ for (int i=0; i<3; ++i) { WN(i) = rightF.cWN(i); } T::ImportData(rightF);}
	inline void Alloc(const int & ns){T::Alloc(ns);}
	inline void Dealloc(){T::Dealloc();}
  static bool HasWedgeNormal()   { return true; }
  static void Name(std::vector<std::string> & name){name.push_back(std::string("WedgeRealNormal"));T::Name(name);}

private:
  NormalType _wn[3];
};

template <class TT> class WedgeRealNormal3s: public WedgeRealNormal<vcg::Point3s, TT> {
public:  static void Name(std::vector<std::string> & name){name.push_back(std::string("WedgeRealNormal2s"));TT::Name(name);}};

template <class TT> class WedgeRealNormal3f: public WedgeRealNormal<vcg::Point3f, TT> {
public:  static void Name(std::vector<std::string> & name){name.push_back(std::string("WedgeRealNormal2f"));TT::Name(name);}};

template <class TT> class WedgeRealNormal3d: public WedgeRealNormal<vcg::Point3d, TT> {
public:  static void Name(std::vector<std::string> & name){name.push_back(std::string("WedgeRealNormal2d"));TT::Name(name);}};


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

template <class T> class EmptyWedgeTexCoord: public T {
public:
	typedef int WedgeTexCoordType;
  typedef vcg::TexCoord2<float,1> TexCoordType;
  TexCoordType &WT(const int) { static TexCoordType dummy_texture;  assert(0); return dummy_texture;}
  TexCoordType const &cWT(const int) const { static TexCoordType dummy_texture; return dummy_texture;}
	template <class RightF>
	void ImportData(const RightF & rightF){ T::ImportData(rightF);}
	inline void Alloc(const int & ns){T::Alloc(ns);}
	inline void Dealloc(){T::Dealloc();}
  static bool HasWedgeTexCoord()   { return false; }
  static bool HasWedgeTexCoordOcc()   { return false; }
  static void Name(std::vector<std::string> & name){T::Name(name);}

};
template <class A, class T> class WedgeTexCoord: public T {
public:
	typedef int WedgeTexCoordType;
  typedef A TexCoordType;
  TexCoordType &WT(const int i) { return _wt[i]; }
  TexCoordType const &cWT(const int i) const { return _wt[i]; }
	template <class RightF>
	void ImportData(const RightF & rightF){ for (int i=0; i<3; ++i) { WT(i) = rightF.cWT(i); } T::ImportData(rightF);}
	inline void Alloc(const int & ns){T::Alloc(ns);}
	inline void Dealloc(){T::Dealloc();}
  static bool HasWedgeTexCoord()   { return true; }
  static void Name(std::vector<std::string> & name){name.push_back(std::string("WedgeTexCoord"));T::Name(name);}

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
  int &Flags() { static int dummyflags(0);  assert(0); return dummyflags; }
  int Flags() const { return 0; }
	template <class RightF>
	void ImportData(const RightF & rightF){ T::ImportData(rightF);}
	inline void Alloc(const int & ns){T::Alloc(ns);}
	inline void Dealloc(){T::Dealloc();}
  static bool HasFlags()   { return false; }
  static bool HasFlagsOcc()   { return false; }
	static void Name(std::vector<std::string> & name){T::Name(name);}

};

template <class T> class BitFlags:  public T {
public:
  BitFlags(){_flags=0;}
   int &Flags() {return _flags; }
  int Flags() const {return _flags; }
  const int & cFlags() const {return _flags; }
	template <class RightF>
	void ImportData(const RightF & rightF){ Flags() = rightF.cFlags();T::ImportData(rightF);}
	inline void Alloc(const int & ns){T::Alloc(ns);}
	inline void Dealloc(){T::Dealloc();}
  static bool HasFlags()   { return true; }
  static void Name(std::vector<std::string> & name){name.push_back(std::string("BitFlags"));T::Name(name);}


private:
  int  _flags;    
};

/*-------------------------- COLOR ----------------------------------*/ 

template <class T> class EmptyColorMarkQuality: public T {
public:
	typedef int MarkType;
	inline void InitIMark()    {  }
	inline int & IMark()       { assert(0); static int tmp=-1; return tmp;}
	inline int IMark() const {return 0;}

	typedef float QualityType;
	typedef vcg::Color4b ColorType;
	ColorType &C() { static ColorType dumcolor(vcg::Color4b::White);  assert(0); return dumcolor; }
	const ColorType &cC() const { static ColorType dumcolor(vcg::Color4b::White);  assert(0); return dumcolor; }
	ColorType &WC(const int) { static ColorType dumcolor(vcg::Color4b::White);  assert(0); return dumcolor; }
	const ColorType &cWC(const int) const { static ColorType dumcolor(vcg::Color4b::White);  assert(0); return dumcolor; }
	QualityType &Q() { static QualityType dummyQuality(0);  assert(0); return dummyQuality; }
	const QualityType &cQ() const { static QualityType dummyQuality(0);  assert(0); return dummyQuality; }
  
	static bool HasFaceColor()   { return false; }
	static bool HasFaceColorOcc() { return false;}
	static bool HasFaceColorOcf() { return false;}

	static bool HasWedgeColor()   { return false; }
	static bool HasWedgeColorOcc()   { return false; }
	static bool HasWedgeColorOcf()   { return false; }

	static bool HasFaceQuality()   { return false; }
	static bool HasFaceQualityOcc()   { return false; }
	static bool HasFaceQualityOcf() { return false;}

	static bool HasMark()   { return false; }
	static bool HasMarkOcc()   { return false; }
	static bool HasMarkOcfb()   { return false; }
  
	static void Name(std::vector<std::string> & name){T::Name(name);}
	template <class RightF>
	void ImportData(const RightF & rightF){ T::ImportData(rightF);}
	inline void Alloc(const int & ns){T::Alloc(ns);}
	inline void Dealloc(){T::Dealloc();}

};
template <class A, class T> class Color: public T {
public:
  typedef A ColorType;
  Color():_color(vcg::Color4b::White) {}
  ColorType &C() { return _color; }
  const ColorType &cC() const { return _color; }
	template <class RightF>
	void ImportData(const RightF & rightF){ C() = rightF.cC();T::ImportData(rightF);}
	inline void Alloc(const int & ns){T::Alloc(ns);}
	inline void Dealloc(){T::Dealloc();}
  static bool HasFaceColor()   { return true; }
  static void Name(std::vector<std::string> & name){name.push_back(std::string("Color"));T::Name(name);}

private:
  ColorType _color;    
};

template <class A, class T> class WedgeColor: public T {
public:
  typedef A ColorType;
  ColorType &WC(const int i) { return _color[i]; }
  const ColorType &WC(const int i) const { return _color[i]; }
  const ColorType &cWC(const int i) const { return _color[i]; }

	template <class RightF>
	void ImportData(const RightF & rightF){ if (RightF::HasWedgeColor()) { for (int i=0; i<3; ++i) { WC(i) = rightF.cWC(i); } T::ImportData(rightF); } }
  static bool HasWedgeColor()   { return true; }
  static void Name(std::vector<std::string> & name){name.push_back(std::string("WedgeColor"));T::Name(name);}

private:
  ColorType _color[3];    
};

template <class T> class WedgeColor4b: public WedgeColor<vcg::Color4b, T> {
public: static void Name(std::vector<std::string> & name){name.push_back(std::string("WedgeColor4b"));T::Name(name);}
};

template <class T> class WedgeColor4f: public WedgeColor<vcg::Color4f, T> {
public: static void Name(std::vector<std::string> & name){name.push_back(std::string("WedgeColor4f"));T::Name(name);}
};

template <class T> class Color4b: public Color<vcg::Color4b, T> { public:
public: static void Name(std::vector<std::string> & name){name.push_back(std::string("Color4b"));T::Name(name);}
};

/*-------------------------- Quality  ----------------------------------*/ 

template <class A, class T> class Quality: public T {
public:
  typedef A QualityType;
  QualityType &Q() { return _quality; }
  const QualityType &cQ() const { return _quality; }
	template <class RightF>
	void ImportData(const RightF & rightF){ if(RightF::HasFaceQuality()) Q() = rightF.cQ();T::ImportData(rightF);}
	inline void Alloc(const int & ns){T::Alloc(ns);}
	inline void Dealloc(){T::Dealloc();}
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

template <class T> class Mark: public T {
public:
  static bool HasMark()      { return true; }
  static bool HasMarkOcc()   { return true; }
  inline void InitIMark()    { _imark = 0; }
  inline int & IMark()       { return _imark;}
  inline const int & IMark() const {return _imark;}
	template <class RightF>
	void ImportData(const RightF & rightF){ IMark() = rightF.IMark();T::ImportData(rightF);}
  static void Name(std::vector<std::string> & name){name.push_back(std::string("Mark"));T::Name(name);}
    
 private:
	int _imark;
};


/*----------------------------- VFADJ ------------------------------*/ 


template <class T> class EmptyAdj: public T {
public:
	typedef int VFAdjType;
  typename T::FacePointer       &VFp(const int)       { static typename T::FacePointer fp=0;  assert(0); return fp; }
  typename T::FacePointer const cVFp(const int) const { static typename T::FacePointer const fp=0; return fp; }
  typename T::FacePointer       &FFp(const int)       { static typename T::FacePointer fp=0;  assert(0); return fp; }
  typename T::FacePointer const cFFp(const int) const { static typename T::FacePointer const fp=0; return fp; }
  typename T::EdgePointer       &FEp(const int)       { static typename T::EdgePointer fp=0;  assert(0); return fp; }
  typename T::EdgePointer const cFEp(const int) const { static typename T::EdgePointer const fp=0; return fp; }
	typename T::HEdgePointer       &FHp()       { static typename T::HEdgePointer fp=0;  assert(0); return fp; }
	typename T::HEdgePointer const cFHp() const { static typename T::HEdgePointer const fp=0; assert(0);return fp; }
	char &VFi(const int j){(void)j; static char z=0;  assert(0); return z;};
  char &FFi(const int j){(void)j; static char z=0;  assert(0); return z;};
  const char &cVFi(const int j){(void)j; static char z=0; return z;};
  const char &cFFi(const int j) const {(void)j; static char z=0; return z;};
	template <class RightF>
	void ImportData(const RightF & rightF){ T::ImportData(rightF);}
	inline void Alloc(const int & ns){T::Alloc(ns);}
	inline void Dealloc(){T::Dealloc();}
  static bool HasVFAdjacency()   {   return false; }
  static bool HasFFAdjacency()   {   return false; }
  static bool HasFEAdjacency()   {   return false; }
	static bool HasFHAdjacency()   {   return false; }

  static bool HasFFAdjacencyOcc()   {   return false; }
  static bool HasVFAdjacencyOcc()   {   return false; }
  static bool HasFEAdjacencyOcc()   {   return false; }
	static bool HasFHAdjacencyOcc()   {   return false; }

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
	template <class RightF>
	void ImportData(const RightF & rightF){T::ImportData(rightF);}
	inline void Alloc(const int & ns){T::Alloc(ns);}
	inline void Dealloc(){T::Dealloc();}
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

	template <class RightF>
	void ImportData(const RightF & rightF){T::ImportData(rightF);}
	inline void Alloc(const int & ns){T::Alloc(ns);}
	inline void Dealloc(){T::Dealloc();}
  static bool HasFFAdjacency()      {   return true; }
  static bool HasFFAdjacencyOcc()   {   return false; }
  static void Name(std::vector<std::string> & name){name.push_back(std::string("FFAdj"));T::Name(name);}

private:
  typename T::FacePointer _ffp[3] ;    
  char _ffi[3] ;    
};


/*----------------------------- FEADJ ------------------------------*/ 

template <class T> class FEAdj: public T {
public:
	FEAdj(){
		_fep[0]=0;
		_fep[1]=0;
		_fep[2]=0;
	}
  typename T::FacePointer       &FEp(const int j)        { assert(j>=0 && j<3);  return _fep[j]; }
  typename T::FacePointer const  FEp(const int j) const  { assert(j>=0 && j<3);  return _fep[j]; }
  typename T::FacePointer const cFEp(const int j) const  { assert(j>=0 && j<3);  return _fep[j]; }
  char        &FEi(const int j)       { return _fei[j]; }
  const char &cFEi(const int j) const { return _fei[j]; }

  typename T::FacePointer        &FEp1( const int j )       { return FEp((j+1)%3);}
	typename T::FacePointer        &FEp2( const int j )       { return FEp((j+2)%3);}
	typename T::FacePointer  const  FEp1( const int j ) const { return FEp((j+1)%3);}
	typename T::FacePointer  const  FEp2( const int j ) const { return FEp((j+2)%3);}

	template <class RightF>
	void ImportData(const RightF & rightF){T::ImportData(rightF);}
	inline void Alloc(const int & ns){T::Alloc(ns);}
	inline void Dealloc(){T::Dealloc();}
  static bool HasFEAdjacency()      {   return true; }
  static void Name(std::vector<std::string> & name){name.push_back(std::string("FEAdj"));T::Name(name);}

private:
  typename T::FacePointer _fep[3] ;    
  char _fei[3] ;    
};


/*----------------------------- FHADJ ------------------------------*/

template <class T> class FHAdj: public T {
public:
	FHAdj(){_fh=0;}
	typename T::HEdgePointer       &FHp( )        { return _fh; }
	typename T::HEdgePointer const cFHp( ) const  { return _fh; }

	template <class RightF>
	void ImportData(const RightF & rightF){T::ImportData(rightF);}
	inline void Alloc(const int & ns){T::Alloc(ns);}
	inline void Dealloc(){T::Dealloc();}
	static bool HasFHAdjacency()      {   return true; }
	static void Name(std::vector<std::string> & name){name.push_back(std::string("FHAdj"));T::Name(name);}

private:
        typename T::HEdgePointer _fh ;
};
  } // end namespace face
}// end namespace vcg
#endif
