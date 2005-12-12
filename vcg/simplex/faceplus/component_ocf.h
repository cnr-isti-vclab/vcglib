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
Revision 1.5  2005/11/26 00:16:44  cignoni
Corrected a lot of bugs about the use of enabled entities

Revision 1.4  2005/11/21 21:46:20  cignoni
Changed HasColor -> HasFaceColor and HasNormal ->HasFaceNormal

Revision 1.3  2005/11/16 22:43:36  cignoni
Added WedgeTexture component

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
  
// Auxiliary types to build internal vectors
struct AdjTypePack {
  typename VALUE_TYPE::FacePointer _fp[3] ;    
  char _zp[3] ;    
  };
  
//template <class TexCoordType>
class WedgeTexTypePack {
public:
  WedgeTexTypePack() { 
    wt[0].u()=.5;wt[0].v()=.5;
    wt[1].u()=.5;wt[1].v()=.5;
    wt[2].u()=.5;wt[2].v()=.5;
    wt[0].n()=0; 
    wt[1].n()=0; 
    wt[2].n()=0; 
  }

  typename VALUE_TYPE::TexCoordType wt[3];
};




  // override di tutte le funzioni che possono spostare 
	// l'allocazione in memoria del container
	void push_back(const VALUE_TYPE & v)
  {
    ThisTypeIterator oldbegin=begin();
    ThisTypeIterator oldend=end();
    BaseType::push_back(v);
    if(oldbegin!=begin()) _updateOVP(begin(),end());
                     else _updateOVP(oldend, end());
  }
	void pop_back();
  void resize(const unsigned int & _size) 
  {
    ThisTypeIterator oldbegin=begin();
    ThisTypeIterator oldend=end();
    BaseType::resize(_size);
    if(oldbegin!=begin()) _updateOVP(begin(),end());
                     else _updateOVP(oldend, end());
    if(ColorEnabled)       CV.resize(_size);
    if(NormalEnabled)      NV.resize(_size);
    if(VFAdjacencyEnabled) AV.resize(_size);
    if(FFAdjacencyEnabled) AF.resize(_size);
    if (WedgeTexEnabled) WTV.resize(_size,WedgeTexTypePack());
    
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
    if(oldbegin!=begin()) _updateOVP(begin(),end());
  }

 void _updateOVP(ThisTypeIterator lbegin, ThisTypeIterator lend)
{
    ThisTypeIterator fi;
    //for(fi=begin();vi!=end();++vi)
    for(fi=lbegin;fi!=lend;++fi)
        (*fi)._ovp=this;
 }
////////////////////////////////////////
// Enabling Eunctions

bool IsColorEnabled() const {return ColorEnabled;}
void EnableColor() {
  assert(VALUE_TYPE::HasFaceColorOcf());
  ColorEnabled=true;
  CV.resize(size());
}

void DisableColor() {
  assert(VALUE_TYPE::HasFaceColorOcf());
  ColorEnabled=false;
  CV.clear();
}

bool IsNormalEnabled() const {return NormalEnabled;}
void EnableNormal() {
  assert(VALUE_TYPE::HasFaceNormalOcf());
  NormalEnabled=true;
  NV.resize(size());
}

void DisableNormal() {
  assert(VALUE_TYPE::HasFaceNormalOcf());
  NormalEnabled=false;
  NV.clear();
}
  
bool IsVFAdjacencyEnabled() const {return VFAdjacencyEnabled;}
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


bool IsFFAdjacencyEnabled() const {return FFAdjacencyEnabled;}
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

bool IsWedgeTexEnabled() const {return WedgeTexEnabled;}
void EnableWedgeTex() {
  assert(VALUE_TYPE::HasWedgeTextureOcf());
  WedgeTexEnabled=true;
  WTV.resize(size(),WedgeTexTypePack());
}

void DisableWedgeTex() {
  assert(VALUE_TYPE::HasWedgeTextureOcf());
  WedgeTexEnabled=false;
  WTV.clear();
}

public:
  std::vector<typename VALUE_TYPE::ColorType> CV;
  std::vector<typename VALUE_TYPE::NormalType> NV;
  std::vector<struct AdjTypePack> AV;
  std::vector<struct AdjTypePack> AF;
  std::vector<class WedgeTexTypePack> WTV;

  bool ColorEnabled;
  bool NormalEnabled;
  bool WedgeTexEnabled;
  bool VFAdjacencyEnabled;
  bool FFAdjacencyEnabled;
}; // end class vector_ocf


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
  static bool HasFaceNormal()   { return true; }
  static bool HasFaceNormalOcf()   { return true; }

  NormalType &N() { 
    // you cannot use Normals before enabling them with: yourmesh.face.EnableNormal()
    assert(Base().NormalEnabled); 
    return Base().NV[Index()];  }
  const NormalType &cN() const { 
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
  static bool HasFaceColor()   { return true; }
  static bool HasFaceColorOcf()   { return true; }
};

template <class T> class Color4bOcf: public ColorOcf<vcg::Color4b, T> {};


///*-------------------------- WEDGE TEXCOORD  ----------------------------------*/ 

template <class A, class TT> class WedgeTextureOcf: public TT {
public:
  WedgeTextureOcf(){ }
  typedef A TexCoordType;
  TexCoordType &WT(const int i)              { assert(Base().WedgeTexEnabled); return Base().WTV[Index()].wt[i]; }
  TexCoordType const &cWT(const int i) const { assert(Base().WedgeTexEnabled); return Base().WTV[Index()].wt[i]; }
  static bool HasWedgeTexture()   { return true; }
  static bool HasWedgeTextureOcf()   { return true; }
};

template <class T> class WedgeTexturefOcf: public WedgeTextureOcf<TCoord2<float,1>, T> {};

///*-------------------------- InfoOpt  ----------------------------------*/ 

template < class T> class InfoOcf: public T {
public:
  vector_ocf<typename T::FaceType> &Base() const { return *_ovp;}

  inline int Index() const {
    typename T::FaceType const *tp=static_cast<typename T::FaceType const *>(this); 
    int tt2=tp- &*(_ovp->begin());
    return tt2;
  } 
public:
  // ovp Optional Vector Pointer
  // Pointer to the base vector where each face element is stored. 
  // used to access to the vectors of the other optional members.
  vector_ocf<typename T::FaceType> *_ovp;
};

  } // end namespace face
  
  template < class, class > class TriMesh;

  namespace tri
  {
    template < class VertContainerType, class FaceType >
      bool HasPerWedgeTexture (const TriMesh < VertContainerType , face::vector_ocf< FaceType > > & m) 
    {
      if(FaceContainerType::HasWedgeTextureOcf()) return m.face.WedgeTexEnabled();
      else return FaceContainerType::HasWedgeTexture();
    }

    template < class VertContainerType, class FaceType >
      bool HasPerFaceColor (const TriMesh < VertContainerType , face::vector_ocf< FaceType > > & m) 
    {
      if(FaceContainerType::HasPerFaceColorOcf()) return m.face.FaceColor();
      else return FaceContainerType::HasFaceColor();
    }
  }
}// end namespace vcg
#endif
