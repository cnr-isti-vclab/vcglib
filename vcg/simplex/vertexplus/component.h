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
Revision 1.11  2005/11/01 18:17:52  cignoni
Added an assert(0) in all the accesses to empty components

Revision 1.10  2005/10/15 16:24:10  ganovelli
Working release (compilata solo su MSVC), component_occ è migrato da component_opt

Revision 1.9  2005/10/14 13:30:07  cignoni
Added constant access functions and reflective functions (HasSomething stuff)
to all the components This is the first really working version...

Revision 1.8  2005/10/07 15:19:54  cignoni
minor updates to keep it in line with the rest of the library

Revision 1.7  2004/05/10 13:50:32  cignoni
Updated names of adj functions to the new standards

Revision 1.6  2004/04/05 11:53:06  cignoni
addend constant access funcs

Revision 1.5  2004/04/03 13:35:51  cignoni
minor changes

Revision 1.4  2004/03/31 13:15:28  cignoni
Added optional cpmponent

Revision 1.3  2004/03/31 12:28:37  ganovelli
*** empty log message ***

Revision 1.2  2004/03/29 14:26:38  cignoni
Error in color

Revision 1.1  2004/03/29 08:36:26  cignoni
First working version!


****************************************************************************/
#ifndef __VCG_VERTEX_PLUS_COMPONENT
#define __VCG_VERTEX_PLUS_COMPONENT
#include <vcg/space/point3.h>
#include <vcg/space/tcoord2.h>
#include <vcg/space/color4.h>

namespace vcg {
  namespace vert {
/*
Some naming Rules
All the Components that can be added to a vertex should be defined in the namespace vert:

*/

/*------------------------- COORD -----------------------------------------*/ 
template <class T> class EmptyCoord: public T {
public:
  typedef vcg::Point3f CoordType;
  typedef CoordType::ScalarType      ScalarType;

  CoordType &P() { static CoordType coord(0, 0, 0); return coord; }
  const CoordType &P() const { static CoordType coord(0, 0, 0);  assert(0); return coord; }
  const CoordType &cP() const { static CoordType coord(0, 0, 0);  assert(0); return coord; }
  CoordType &UberP() { static CoordType coord(0, 0, 0); return coord; }
  static bool HasCoord()   { return false; }

};

template <class A, class T> class Coord: public T {
public:
  typedef A CoordType;
  typedef typename CoordType::ScalarType      ScalarType;
  CoordType &P() { return _coord; }
  const CoordType &P() const { return _coord; }
  const CoordType &cP() const { return _coord; }
  CoordType &UberP() { return _coord; }
  
  static bool HasCoord()   { return true; }
private:
  CoordType _coord;    
};
template <class T> class Coord3f: public Coord<vcg::Point3f, T> {};
template <class T> class Coord3d: public Coord<vcg::Point3d, T> {};

/*-------------------------- NORMAL ----------------------------------------*/ 

template <class T> class EmptyNormal: public T {
public:
  typedef vcg::Point3s NormalType;
  NormalType &N() { static NormalType dummy_normal(0, 0, 0);  assert(0); return dummy_normal; }
  const NormalType cN()const { static NormalType dummy_normal(0, 0, 0);  assert(0); return dummy_normal; }
  static bool HasNormal()   { return false; }
  static bool HasNormalOpt()   { return false; }
};
template <class A, class T> class Normal: public T {
public:
  typedef A NormalType;
  NormalType &N() { return _norm; }
  const NormalType cN() const { return _norm; }
  static bool HasNormal()   { return true; }
private:
  NormalType _norm;    
};

template <class T> class Normal3s: public Normal<vcg::Point3s, T> {};
template <class T> class Normal3f: public Normal<vcg::Point3f, T> {};
template <class T> class Normal3d: public Normal<vcg::Point3d, T> {};

/*-------------------------- NORMAL ----------------------------------------*/ 

template <class TT> class EmptyTexture: public TT {
public:
  typedef vcg::TCoord2<float,1> TextureType;
  TextureType &T() { static TextureType dummy_texture;  assert(0); return dummy_texture; }
  static bool HasTexture()   { return false; }
  static bool HasTextureOpt()   { return false; }

};
template <class A, class TT> class Texture: public TT {
public:
  typedef A TextureType;
  TextureType &T() { return _t; }
  static bool HasTexture()   { return true; }

private:
  TextureType _t;    
};

template <class TT> class Texture2s: public Texture<TCoord2<short,1>, TT> {};
template <class TT> class Texture2f: public Texture<TCoord2<float,1>, TT> {};
template <class TT> class Texture2d: public Texture<TCoord2<double,1>, TT> {};

/*------------------------- FLAGS -----------------------------------------*/ 
template <class T> class EmptyFlag: public T {
public:
	typedef int FlagType;
  /// Return the vector of Flags(), senza effettuare controlli sui bit
  int &Flags() { static int dummyflags(0);  assert(0); return dummyflags; }
  const int Flags() const { return 0; }
  static bool HasFlags()   { return false; }

};

template <class T> class Flag:  public T {
public:
	typedef int FlagType;
  int &Flags() {return _flags; }
  const int Flags() const {return _flags; }
  static bool HasFlag()   { return true; }

private:
  int  _flags;    
};

/*-------------------------- COLOR ----------------------------------*/ 

template <class T> class EmptyColor: public T {
public:
  typedef vcg::Color4b ColorType;
  ColorType &C() { static ColorType dumcolor(vcg::Color4b::White); assert(0); return dumcolor; }
  static bool HasColor()   { return false; }
};
template <class A, class T> class Color: public T {
public:
  typedef A ColorType;
  ColorType &C() { return _color; }
  static bool HasColor()   { return true; }
private:
  ColorType _color;    
};

template <class T> class Color4b: public Color<vcg::Color4b, T> {};

/*-------------------------- Quality  ----------------------------------*/ 

template <class T> class EmptyQuality: public T {
public:
  typedef float QualityType;
  QualityType &Q() { static QualityType dummyQuality(0);  assert(0); return dummyQuality; }
  static bool HasQuality()   { return false; }
};
template <class A, class T> class Quality: public T {
public:
  typedef A QualityType;
  QualityType &Q() { return _quality; }
  static bool HasQuality()   { return true; }

private:
  QualityType _quality;    
};

template <class T> class Qualitys: public Quality<short, T> {};
template <class T> class Qualityf: public Quality<float, T> {};
template <class T> class Qualityd: public Quality<double, T> {};

/*----------------------------- VFADJ ------------------------------*/ 


template <class T> class EmptyVFAdj: public T {
public:
  typename T::FacePointer &VFp() { static typename T::FacePointer fp=0;  assert(0); return fp; }
  typename T::FacePointer cVFp() { static typename T::FacePointer fp=0;  assert(0); return fp; }
  int &VFi(){static int z=0; return z;};
  static bool HasVFAdjacency()   {   return false; }
  static bool HasVFAdjacencyOpt()   {   return false; }
};

template <class T> class VFAdj: public T {
public:
  typename T::FacePointer &VFp() {return _fp; }
  typename T::FacePointer cVFp() {return _fp; }
  int &VFi() {return _zp; }
  static bool HasVFAdjacency()   {   return true; }
  static bool HasVFAdjacencyOpt()   {   return false; }
private:
  typename T::FacePointer _fp ;    
  int _zp ;    
};

  } // end namespace vert
}// end namespace vcg
#endif
