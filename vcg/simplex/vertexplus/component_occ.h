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
Revision 1.1  2005/10/15 16:24:10  ganovelli
Working release (compilata solo su MSVC), component_occ è migrato da component_opt



****************************************************************************/
#ifndef __VCG_VERTEX_PLUS_COMPONENT_OCC
#define __VCG_VERTEX_PLUS_COMPONENT_OCC

#include <vcg/simplex/vertexplus/component.h>
#include <vcg/container/vector_occ.h>


namespace vcg {
  namespace vert {
/*
Some naming Rules
All the Components that can be added to a vertex should be defined in the namespace vert:

*/

/*------------------------- COORD -----------------------------------------*/ 

template <class A, class T> class CoordOcc: public T {
public:
  typedef A CoordType;
  typedef typename CoordType::ScalarType      ScalarType;
	CoordType &P() { return CAT< vector_occ<VertType>,CoordType>::Instance()->Get((VertType*)this); }
  CoordType &UberP() { return CAT< vector_occ<VertType>,CoordType>::Instance()->Get((VertType*)this); }
};
template <class T> class Coord3fOcc: public CoordOcc<vcg::Point3f, T> {};
template <class T> class Coord3dOcc: public CoordOcc<vcg::Point3d, T> {};


/*-------------------------- NORMAL ----------------------------------------*/ 

template <class A, class T> class NormalOcc: public T {
public:
  typedef A NormalType;
  NormalType &N() {return CAT< vector_occ<VertType>,NormalType>::Instance()->Get((VertType*)this); }
/*private:
  NormalType _norm;   */ 
};

template <class T> class Normal3sOcc: public NormalOcc<vcg::Point3s, T> {};
template <class T> class Normal3fOcc: public NormalOcc<vcg::Point3f, T> {};
template <class T> class Normal3dOcc: public NormalOcc<vcg::Point3d, T> {};

/*-------------------------- TEXCOORD ----------------------------------------*/ 

template <class A, class T> class TexCoordOcc: public T {
public:
  typedef A TexCoordType;
  TexCoordType &T() {return CAT< vector_occ<VertType>,TexCoordType>::Instance()->Get((VertType*)this); }
  static bool HasTexCoord()   { return true; }

/* private:
  TexCoordType _t;   */ 
};

template <class T> class TexCoord2sOcc: public TexCoordOcc<TexCoord2<short,1>, T> {};
template <class T> class TexCoord2fOcc: public TexCoordOcc<TexCoord2<float,1>, T> {};
template <class T> class TexCoord2dOcc: public TexCoordOcc<TexCoord2<double,1>, T> {};

///*------------------------- FLAGS -----------------------------------------*/ 

template <class T> class FlagOcc:  public T {
public:
   int &Flags() {return CAT< vector_occ<VertType>,int>::Instance()->Get((VertType*)this); }
   const int Flags() const {return _flags; }
/*
private:
  int  _flags;  */ 
};

///*-------------------------- COLOR ----------------------------------*/ 

template <class A, class T> class ColorOcc: public T {
public:
  typedef A ColorType;
  ColorType &C() { return CAT< vector_occ<VertType>,ColorType>::Instance()->Get((VertType*)this); }
  static bool HasColor()   { return true; }
/*private:
  ColorType _color;   */ 
};

template <class T> class Color4bOcc: public ColorOcc<vcg::Color4b, T> {};

///*-------------------------- Quality  ----------------------------------*/ 

template <class A, class T> class QualityOcc: public T {
public:
  typedef A QualityType;
  QualityType &Q() { return CAT< vector_occ<VertType>,QualityType>::Instance()->Get((VertType*)this);}
  static bool HasQuality()   { return true; }

/*private:
  QualityType _quality;  */  
};

template <class T> class QualitysOcc: public QualityOcc<short, T> {};
template <class T> class QualityfOcc: public QualityOcc<float, T> {};
template <class T> class QualitydOcc: public QualityOcc<double, T> {};
//
///*----------------------------- VFADJ ------------------------------*/ 

template <class T> class VFAdjOcc: public T {
public:
  typename T::FacePointer &Fp() {return CAT< vector_occ<VertType>,T::FacePointer>::Instance()->Get((VertType*)this); }
  int &Zp() {return _zp; }
  static bool HasVFAdjacency()   {   return true; }
private:
  typename T::FacePointer _fp ;    
  int _zp ;    
};

  } // end namespace vert
}// end namespace vcg
#endif
