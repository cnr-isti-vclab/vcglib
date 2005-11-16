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
Revision 1.2  2005/10/22 13:16:46  cignoni
Added a missing ';' in  FFAdjOcf (thanks to Mario Latronico).

Revision 1.1  2005/10/14 15:07:58  cignoni
First Really Working version


****************************************************************************/

/* 
Note
OCF = Optional Component Fast (hopefully)
compare with OCC(Optional Component Compact)

Mainly the trick here is to store a base pointer in each simplex...

****************************************************************************/
#ifndef __VCG_FACE_PLUS_COMPONENT_OCF
#define __VCG_FACE_PLUS_COMPONENT_OCF

#include <vcg/simplex/faceplus/component.h>
#include <vector>


namespace vcg {
  namespace face {
/*
All the Components that can be added to a faceex should be defined in the namespace face:

*/
template <class VALUE_TYPE>
class vector_ocf: public std::vector<VALUE_TYPE> {
  typedef std::vector<VALUE_TYPE> BaseType; 
	typedef typename vector_ocf<VALUE_TYPE>::iterator ThisTypeIterator;
  
public:
	vector_ocf():std::vector<VALUE_TYPE>(){
  ColorEnabled=false;
  NormalEnabled=false;
  WedgeTexEnabled=false;
  VFAdjacencyEnabled=false;
  FFAdjacencyEnabled=false;
  }
  
  // override di tutte le funzioni che possono spostare 
	// l'allocazione in memoria del container
	void push_back(const VALUE_TYPE & v)
  {
    ThisTypeIterator oldbegin=begin();
    ThisTypeIterator oldend=end();
    BaseType::push_back(v);
    if(oldbegin!=begin()) _update(begin(),end());
    else _update(oldend,end());
  }
	void pop_back();
  void resize(const unsigned int & _size) 
  {
    ThisTypeIterator oldbegin=begin();
    ThisTypeIterator oldend=end();
    BaseType::resize(_size);
    if(oldbegin!=begin()) _update(begin(),end());
    else _update(oldend,end());
    if(ColorEnabled)       CV.resize(_size);
    if(NormalEnabled)      NV.resize(_size);
    if(VFAdjacencyEnabled) AV.resize(_size);
    if(FFAdjacencyEnabled) AF.resize(_size);
    if (WedgeTexEnabled) WTV.resize(_size);
    
   }
  void reserve(const unsigned int & _size)
  {
    ThisTypeIterator oldbegin=begin();
    BaseType::reserve(_size);
    if (ColorEnabled)       CV.reserve(_size);
    if (NormalEnabled)      NV.reserve(_size);
    if (VFAdjacencyEnabled) AV.reserve(_size);
    if (FFAdjacencyEnabled) AF.reserve(_size);
    if (WedgeTexEnabled) WTV.reserve(_size);
    if(oldbegin!=begin()) _update(begin(),end());
  }

 void _update(ThisTypeIterator lbegin, ThisTypeIterator lend)
{
    ThisTypeIterator vi;
    //for(vi=lbegin;vi!=lend;++vi)
    for(vi=begin();vi!=end();++vi)
        (*vi).EV=this;
 }
////////////////////////////////////////
// Enabling Eunctions

void EnableColor() {
  assert(VALUE_TYPE::HasColorOcf());
  ColorEnabled=true;
  CV.resize(size());
}

void DisableColor() {
  assert(VALUE_TYPE::HasColorOcf());
  ColorEnabled=false;
  CV.clear();
}

void EnableNormal() {
  assert(VALUE_TYPE::HasNormalOcf());
  NormalEnabled=true;
  NV.resize(size());
}

void DisableNormal() {
  assert(VALUE_TYPE::HasNormalOcf());
  NormalEnabled=false;
  NV.clear();
}
  
void EnableVFAdjacency() {
  assert(VALUE_TYPE::HasVFAdjacencyOcf());
  VFAdjacencyEnabled=true;
  AV.resize(size());
}

void DisableVFAdjacency() {
  assert(VALUE_TYPE::HasVFAdjacencyOcf());
  VFAdjacencyEnabled=false;
  AV.clear();
}


void EnableFFAdjacency() {
  assert(VALUE_TYPE::HasFFAdjacencyOcf());
  FFAdjacencyEnabled=true;
  AF.resize(size());
}

void DisableFFAdjacency() {
  assert(VALUE_TYPE::HasFFAdjacencyOcf());
  FFAdjacencyEnabled=false;
  AF.clear();
}


void EnableWedgeTex() {
  assert(VALUE_TYPE::HasWedgeTexture());
  WedgeTexEnabled=true;
  WTV.resize(size());
}

void DisableWedgeTex() {
  assert(VALUE_TYPE::HasWedgeTexture());
  WedgeTexEnabled=false;
  WTV.clear();
}

struct AdjType {
  typename VALUE_TYPE::FacePointer _fp[3] ;    
  char _zp[3] ;    
  };
  
public:
  std::vector<typename VALUE_TYPE::ColorType> CV;
  std::vector<typename VALUE_TYPE::NormalType> NV;
  std::vector<struct AdjType> AV;
  std::vector<struct AdjType> AF;
  std::vector<typename VALUE_TYPE::TexCoordType> WTV;

  bool ColorEnabled;
  bool NormalEnabled;
  bool WedgeTexEnabled;
  bool VFAdjacencyEnabled;
  bool FFAdjacencyEnabled;
};


//template<>	void EnableAttribute<typename VALUE_TYPE::NormalType>(){	NormalEnabled=true;}

/*------------------------- COORD -----------------------------------------*/ 
/*----------------------------- VFADJ ------------------------------*/ 


template <class T> class VFAdjOcf: public T {
public:
  typename T::FacePointer &VFp(const int j) {
    assert(Base().VFAdjacencyEnabled); 
    return Base().AV[Index()]._fp[j]; 
  }

  typename T::FacePointer cVFp(const int j) const {
    if(! Base().VFAdjacencyEnabled ) return 0; 
    else return Base().AV[Index()]._fp[j]; 
  }

  char &VFi(const int j) {
    assert(Base().VFAdjacencyEnabled); 
    return Base().AV[Index()]._zp[j]; 
  }
  static bool HasVFAdjacency()   {   return true; }
  static bool HasVFAdjacencyOcf()   { return true; }

private:
};
/*----------------------------- FFADJ ------------------------------*/ 


template <class T> class FFAdjOcf: public T {
public:
  typename T::FacePointer &FFp(const int j) {
    assert(Base().FFAdjacencyEnabled); 
    return Base().AF[Index()]._fp[j]; 
  }

  typename T::FacePointer const  FFp(const int j) const { return cFFp(j);}
  typename T::FacePointer const cFFp(const int j) const {
    if(! Base().FFAdjacencyEnabled ) return 0; 
    else return Base().AF[Index()]._fp[j]; 
  }

  char &FFi(const int j) {
    assert(Base().FFAdjacencyEnabled); 
    return Base().AF[Index()]._zp[j]; 
  }
  static bool HasFFAdjacency()   {   return true; }
  static bool HasFFAdjacencyOcf()   { return true; }

private:
};

/*------------------------- Normal -----------------------------------------*/ 

template <class A, class T> class NormalOcf: public T {
public:
  typedef A NormalType;
  static bool HasNormal()   { return true; }
  static bool HasNormalOcf()   { return true; }

  NormalType &N() { 
    // you cannot use Normals before enabling them with: yourmesh.face.EnableNormal()
    assert(Base().NormalEnabled); 
    return Base().NV[Index()];  }
};

template <class T> class Normal3sOcf: public NormalOcf<vcg::Point3s, T> {};
template <class T> class Normal3fOcf: public NormalOcf<vcg::Point3f, T> {};
template <class T> class Normal3dOcf: public NormalOcf<vcg::Point3d, T> {};

///*-------------------------- COLOR ----------------------------------*/ 

template <class A, class T> class ColorOcf: public T {
public:
  typedef A ColorType;
  ColorType &C() { assert(Base().NormalEnabled); return Base().CV[Index()]; }
  static bool HasColor()   { return true; }
  static bool HasColorOcf()   { return true; }
};

template <class T> class Color4bOcf: public ColorOcf<vcg::Color4b, T> {};


///*-------------------------- WEDGE TEXCOORD  ----------------------------------*/ 

template <class A, class TT> class WedgeTextureOcf: public TT {
public:
  typedef A TexCoordType;
  TexCoordType &WT(const int i)              { assert(Base().WedgeTexEnabled); return Base().WTV[Index()]; }
  TexCoordType const &cWT(const int i) const { assert(Base().WedgeTexEnabled); return Base().WTV[Index()]; }
  static bool HasWedgeTexture()   { return true; }
};

template <class T> class WedgeTexturefOcf: public WedgeTextureOcf<TCoord2<float,1>, T> {};

///*-------------------------- InfoOpt  ----------------------------------*/ 

template < class T> class InfoOcf: public T {
public:
  vector_ocf<typename T::FaceType> &Base() const { return *EV;}

  inline int Index() const {
    typename T::FaceType const *tp=static_cast<typename T::FaceType const *>(this); 
    int tt2=tp- &*(EV->begin());
    return tt2;
  } 
public:
  vector_ocf<typename T::FaceType> *EV;
};


  } // end namespace face
}// end namespace vcg
#endif
