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
Revision 1.2  2004/03/29 14:26:38  cignoni
Error in color

Revision 1.1  2004/03/29 08:36:26  cignoni
First working version!


****************************************************************************/
#ifndef __VCG_VERTEX_PLUS_COMPONENT
#define __VCG_VERTEX_PLUS_COMPONENT

#include <vcg/traced_vector.h>


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
  CoordType &UberP() { static CoordType coord(0, 0, 0); return coord; }
};

template <class A, class T> class Coord: public T {
public:
  typedef A CoordType;
  typedef typename CoordType::ScalarType      ScalarType;
  CoordType &P() { return _coord; }
  CoordType &UberP() { return _coord; }
private:
  CoordType _coord;    
};
template <class T> class Coord3f: public Coord<vcg::Point3f, T> {};
template <class T> class Coord3d: public Coord<vcg::Point3d, T> {};

template <class A, class T> class CoordOpt: public T {
public:
  typedef A CoordType;
  typedef typename CoordType::ScalarType      ScalarType;
	CoordType &P() { return CAT< TVector<VertType>,CoordType>::Get((VertType*)this); }
  CoordType &UberP() { return CAT< TVector<VertType>,CoordType>::Get((VertType*)this); }
};
template <class T> class Coord3fOpt: public CoordOpt<vcg::Point3f, T> {};
template <class T> class Coord3dOpt: public CoordOpt<vcg::Point3d, T> {};


/*-------------------------- NORMAL ----------------------------------------*/ 

template <class T> class EmptyNormal: public T {
public:
  typedef vcg::Point3s NormalType;
  NormalType &N() { static NormalType dummy_normal(0, 0, 0); return dummy_normal; }
};
template <class A, class T> class Normal: public T {
public:
  typedef A NormalType;
  NormalType &N() { return _norm; }
private:
  NormalType _norm;    
};

template <class T> class Normal3s: public Normal<vcg::Point3s, T> {};
template <class T> class Normal3f: public Normal<vcg::Point3f, T> {};
template <class T> class Normal3d: public Normal<vcg::Point3d, T> {};

/*-------------------------- NORMAL ----------------------------------------*/ 

template <class T> class EmptyTexture: public T {
public:
  typedef vcg::TCoord2<float,1> TextureType;
  TextureType &T() { static TextureType dummy_texture; return dummy_texture; }
  static bool HasTexture()   { return false; }
};
template <class A, class T> class Texture: public T {
public:
  typedef A TextureType;
  TextureType &T() { return _t; }
  static bool HasTexture()   { return true; }

private:
  TextureType _t;    
};

template <class T> class Texture2s: public Texture<TCoord2<short,1>, T> {};
template <class T> class Texture2f: public Texture<TCoord2<float,1>, T> {};
template <class T> class Texture2d: public Texture<TCoord2<double,1>, T> {};

/*------------------------- FLAGS -----------------------------------------*/ 
template <class T> class EmptyFlag: public T {
public:
	/// Return the vector of Flags(), senza effettuare controlli sui bit
  int &Flags() { static int dummyflags(0); return dummyflags; }
  const int Flags() const { return 0; }
};

template <class T> class Flag:  public T {
public:
   int &Flags() {return _flags; }
   const int Flags() const {return _flags; }

private:
  int  _flags;    
};

/*-------------------------- COLOR ----------------------------------*/ 

template <class T> class EmptyColor: public T {
public:
  typedef vcg::Point3f ColorType;
  ColorType &C() { static ColorType color(0, 0, 0); return color; }
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
  QualityType &Q() { static QualityType dummyQuality(0); return dummyQuality; }
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
  typename T::FacePointer &Fp() { static typename T::FacePointer fp=0; return fp; }
  int &Zp(){static int z=0; return z;};
  static bool HasVFAdjacency()   {   return false; }
};

template <class T> class VFAdj: public T {
public:
  //typedef A ColorType;
  typename T::FacePointer &Fp() {return fp; }
  typename T::FacePointer &Zp() {return Zp; }
  static bool HasVFAdjacency()   {   return true; }
private:
  typename T::FacePointer fp ;    
  int zp ;    
};

  } // end namespace vert
}// end namespace vcg
#endif