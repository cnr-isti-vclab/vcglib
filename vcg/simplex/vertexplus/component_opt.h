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
Revision 1.3  2004/03/31 14:16:40  ganovelli
Data structure to handle temporary attributes. First version

Revision 1.2  2004/03/31 13:15:28  cignoni
Added optional cpmponent

Revision 1.1  2004/03/31 12:46:53  cignoni
First working version!


****************************************************************************/
#ifndef __VCG_VERTEX_PLUS_COMPONENT_OPT
#define __VCG_VERTEX_PLUS_COMPONENT_OPT

#include <vcg/simplex/vertexplus/component.h>
#include <vcg/container/traced_vector.h>


namespace vcg {
  namespace vert {
/*
Some naming Rules
All the Components that can be added to a vertex should be defined in the namespace vert:

*/

/*------------------------- COORD -----------------------------------------*/ 

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

template <class A, class T> class NormalOpt: public T {
public:
  typedef A NormalType;
  NormalType &N() {return CAT< TVector<VertType>,NormalType>::Get((VertType*)this); }
private:
  NormalType _norm;    
};

template <class T> class Normal3sOpt: public NormalOpt<vcg::Point3s, T> {};
template <class T> class Normal3fOpt: public NormalOpt<vcg::Point3f, T> {};
template <class T> class Normal3dOpt: public NormalOpt<vcg::Point3d, T> {};

/*-------------------------- TEXTURE ----------------------------------------*/ 

template <class A, class T> class TextureOpt: public T {
public:
  typedef A TextureType;
  TextureType &T() {return CAT< TVector<VertType>,TextureType>::Get((VertType*)this); }
  static bool HasTexture()   { return true; }

private:
  TextureType _t;    
};

template <class T> class Texture2sOpt: public TextureOpt<TCoord2<short,1>, T> {};
template <class T> class Texture2fOpt: public TextureOpt<TCoord2<float,1>, T> {};
template <class T> class Texture2dOpt: public TextureOpt<TCoord2<double,1>, T> {};

///*------------------------- FLAGS -----------------------------------------*/ 

template <class T> class FlagOpt:  public T {
public:
   int &Flags() {return CAT< TVector<VertType>,int>::Get((VertType*)this); }
   const int Flags() const {return _flags; }

private:
  int  _flags;    
};

///*-------------------------- COLOR ----------------------------------*/ 

template <class A, class T> class ColorOpt: public T {
public:
  typedef A ColorType;
  ColorType &C() { return CAT< TVector<VertType>,ColorType>::Get((VertType*)this); }
  static bool HasColor()   { return true; }
private:
  ColorType _color;    
};

template <class T> class Color4bOpt: public ColorOpt<vcg::Color4b, T> {};

///*-------------------------- Quality  ----------------------------------*/ 

template <class A, class T> class QualityOpt: public T {
public:
  typedef A QualityType;
  QualityType &Q() { return CAT< TVector<VertType>,QualityType>::Get((VertType*)this);}
  static bool HasQuality()   { return true; }

private:
  QualityType _quality;    
};

template <class T> class QualitysOpt: public QualityOpt<short, T> {};
template <class T> class QualityfOpt: public QualityOpt<float, T> {};
template <class T> class QualitydOpt: public QualityOpt<double, T> {};
//
///*----------------------------- VFADJ ------------------------------*/ 

template <class T> class VFAdjOpt: public T {
public:
  typename T::FacePointer &Fp() {return CAT< TVector<VertType>,T::FacePointer>::Get((VertType*)this); }
  int &Zp() {return _zp; }
  static bool HasVFAdjacency()   {   return true; }
private:
  typename T::FacePointer _fp ;    
  int _zp ;    
};

  } // end namespace vert
}// end namespace vcg
#endif