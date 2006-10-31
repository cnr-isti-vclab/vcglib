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
Revision 1.16  2006/10/27 14:15:10  ganovelli
added overrides to HasFFAddAdjacency and HasVFAddAdjacency

Revision 1.15  2006/10/16 08:49:29  cignoni
Better managment of resize overloading when reducing the size of a vector

Revision 1.14  2006/10/09 20:20:55  cignoni
Added some missing Add***Ocf() for the default case.

Revision 1.13  2006/05/25 09:39:09  cignoni
missing std and other gcc detected syntax errors

Revision 1.12  2006/05/03 21:37:02  cignoni
Added Optional Mark

Revision 1.11  2006/02/28 11:59:39  ponchio
g++ compliance:

begin() -> (*this).begin() and for end(), size(), Base(), Index()

Revision 1.10  2006/01/30 08:47:40  cignoni
Corrected HasPerWedgeTexture

Revision 1.9  2006/01/05 15:46:06  cignoni
Removed a syntax error (double >) in HasPerWedgeTexture/HasPerFaceColor

Revision 1.8  2006/01/04 18:46:25  cignoni
Corrected push_back (did not worked at all!)
added missing cFFi

Revision 1.7  2006/01/03 10:54:21  cignoni
Corrected HasPerFaceColor and HasPerWedgeTexture to comply gcc

Revision 1.6  2005/12/12 11:17:32  cignoni
Corrected update function, now only the needed simplexes should be updated.

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
  MarkEnabled=false;
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
    BaseType::push_back(v);
	BaseType::back()._ovp = this;
    if (ColorEnabled)       CV.push_back(vcg::Color4b(vcg::Color4b::White));
    if (MarkEnabled)        MV.push_back(0);
    if (NormalEnabled)      NV.push_back(typename VALUE_TYPE::NormalType());
    if (VFAdjacencyEnabled) AV.push_back(AdjTypePack());
    if (FFAdjacencyEnabled) AF.push_back(AdjTypePack());
    if (WedgeTexEnabled)   WTV.push_back(WedgeTexTypePack());
  }
	void pop_back();
  void resize(const unsigned int & _size) 
  {
	int oldsize = BaseType::size();
    BaseType::resize(_size);
	if(oldsize<_size){
		ThisTypeIterator firstnew = BaseType::begin();
		advance(firstnew,oldsize);
		_updateOVP(firstnew,(*this).end());
	}  
    if (ColorEnabled)       CV.resize(_size);
    if (MarkEnabled)        MV.resize(_size);
    if (NormalEnabled)      NV.resize(_size);
    if (VFAdjacencyEnabled) AV.resize(_size);
    if (FFAdjacencyEnabled) AF.resize(_size);
    if (WedgeTexEnabled) WTV.resize(_size,WedgeTexTypePack());
    
   }
  void reserve(const unsigned int & _size)
  {
    ThisTypeIterator oldbegin=(*this).begin();
    BaseType::reserve(_size);

    if (ColorEnabled)       CV.reserve(_size);
    if (MarkEnabled)        MV.reserve(_size);
    if (NormalEnabled)      NV.reserve(_size);
    if (VFAdjacencyEnabled) AV.reserve(_size);
    if (FFAdjacencyEnabled) AF.reserve(_size);
    if (WedgeTexEnabled)   WTV.reserve(_size);

    if(oldbegin!=(*this).begin()) _updateOVP((*this).begin(),(*this).end());
  }

 void _updateOVP(ThisTypeIterator lbegin, ThisTypeIterator lend)
{
    ThisTypeIterator fi;
    //for(fi=(*this).begin();vi!=(*this).end();++vi)
    for(fi=lbegin;fi!=lend;++fi)
        (*fi)._ovp=this;
 }
////////////////////////////////////////
// Enabling Functions

bool IsColorEnabled() const {return ColorEnabled;}
void EnableColor() {
  assert(VALUE_TYPE::HasFaceColorOcf());
  ColorEnabled=true;
  CV.resize((*this).size());
}

void DisableColor() {
  assert(VALUE_TYPE::HasFaceColorOcf());
  ColorEnabled=false;
  CV.clear();
}

bool IsMarkEnabled() const {return MarkEnabled;}
void EnableMark() {
  assert(VALUE_TYPE::HasFaceMarkOcf());
  MarkEnabled=true;
  MV.resize((*this).size());
}

void DisableMark() {
  assert(VALUE_TYPE::HasFaceMarkOcf());
  MarkEnabled=false;
  MV.clear();
}

bool IsNormalEnabled() const {return NormalEnabled;}
void EnableNormal() {
  assert(VALUE_TYPE::HasFaceNormalOcf());
  NormalEnabled=true;
  NV.resize((*this).size());
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
  AV.resize((*this).size());
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
  AF.resize((*this).size());
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
  WTV.resize((*this).size(),WedgeTexTypePack());
}

void DisableWedgeTex() {
  assert(VALUE_TYPE::HasWedgeTextureOcf());
  WedgeTexEnabled=false;
  WTV.clear();
}

public:
  std::vector<typename VALUE_TYPE::ColorType> CV;
  std::vector<int> MV;
  std::vector<typename VALUE_TYPE::NormalType> NV;
  std::vector<struct AdjTypePack> AV;
  std::vector<struct AdjTypePack> AF;
  std::vector<class WedgeTexTypePack> WTV;

  bool ColorEnabled;
  bool MarkEnabled;
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
    assert((*this).Base().VFAdjacencyEnabled); 
    return (*this).Base().AV[(*this).Index()]._fp[j]; 
  }

  typename T::FacePointer cVFp(const int j) const {
    if(! (*this).Base().VFAdjacencyEnabled ) return 0; 
    else return (*this).Base().AV[(*this).Index()]._fp[j]; 
  }

  char &VFi(const int j) {
    assert((*this).Base().VFAdjacencyEnabled); 
    return (*this).Base().AV[(*this).Index()]._zp[j]; 
  }
  static bool HasVFAdjacency()   {   return true; }
  static bool HasVFAdjacencyOcf()   { return true; }

private:
};
/*----------------------------- FFADJ ------------------------------*/ 


template <class T> class FFAdjOcf: public T {
public:
  typename T::FacePointer &FFp(const int j) {
    assert((*this).Base().FFAdjacencyEnabled); 
    return (*this).Base().AF[(*this).Index()]._fp[j]; 
  }

  typename T::FacePointer const  FFp(const int j) const { return cFFp(j);}
  typename T::FacePointer const cFFp(const int j) const {
    if(! (*this).Base().FFAdjacencyEnabled ) return 0; 
    else return (*this).Base().AF[(*this).Index()]._fp[j]; 
  }

  char &FFi(const int j) {
    assert((*this).Base().FFAdjacencyEnabled); 
    return (*this).Base().AF[(*this).Index()]._zp[j]; 
  }
  const char cFFi(const int j) const {
    assert((*this).Base().FFAdjacencyEnabled); 
    return (*this).Base().AF[(*this).Index()]._zp[j]; 
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
    assert((*this).Base().NormalEnabled); 
    return (*this).Base().NV[(*this).Index()];  }
  const NormalType &cN() const { 
    // you cannot use Normals before enabling them with: yourmesh.face.EnableNormal()
    assert((*this).Base().NormalEnabled); 
    return (*this).Base().NV[(*this).Index()];  }


};

template <class T> class Normal3sOcf: public NormalOcf<vcg::Point3s, T> {};
template <class T> class Normal3fOcf: public NormalOcf<vcg::Point3f, T> {};
template <class T> class Normal3dOcf: public NormalOcf<vcg::Point3d, T> {};

///*-------------------------- COLOR ----------------------------------*/ 

template <class A, class T> class ColorOcf: public T {
public:
  typedef A ColorType;
  ColorType &C() { 
    assert((*this).Base().ColorEnabled); 
    return (*this).Base().CV[(*this).Index()]; 
  }
  static bool HasFaceColor()   { return true; }
  static bool HasFaceColorOcf()   { return true; }
};

template <class T> class Color4bOcf: public ColorOcf<vcg::Color4b, T> {};

///*-------------------------- MARK  ----------------------------------*/ 

template <class T> class MarkOcf: public T {
public:
  inline int & IMark()       { 
    assert((*this).Base().MarkEnabled); 
    return (*this).Base().MV[(*this).Index()]; 
  }
 
  inline int IMark() const   { 
    assert((*this).Base().MarkEnabled); 
    return (*this).Base().MV[(*this).Index()]; 
  } ;
 
  static bool HasFaceMark()   { return true; }
  static bool HasFaceMarkOcf()   { return true; }
  inline void InitIMark()    { IMark() = 0; }
};

///*-------------------------- WEDGE TEXCOORD  ----------------------------------*/ 

template <class A, class TT> class WedgeTextureOcf: public TT {
public:
  WedgeTextureOcf(){ }
  typedef A TexCoordType;
  TexCoordType &WT(const int i)              { assert((*this).Base().WedgeTexEnabled); return (*this).Base().WTV[(*this).Index()].wt[i]; }
  TexCoordType const &cWT(const int i) const { assert((*this).Base().WedgeTexEnabled); return (*this).Base().WTV[(*this).Index()].wt[i]; }
  static bool HasWedgeTexture()   { return true; }
  static bool HasWedgeTextureOcf()   { return true; }
};

template <class T> class WedgeTexturefOcf: public WedgeTextureOcf<TCoord2<float,1>, T> {};

///*-------------------------- InfoOpt  ----------------------------------*/ 

template < class T> class InfoOcf: public T {
public:
  vector_ocf<typename T::FaceType> &Base() const { return *_ovp;}

  static bool HasFaceColorOcf()   { return false; }
  static bool HasFaceNormalOcf()   { return false; }
  static bool HasFaceMarkOcf()   { return false; }
  static bool HasWedgeTextureOcf()   { return false; }

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
	  bool HasVFAdjacency (const TriMesh < VertContainerType , face::vector_ocf< FaceType > > & m) 
	{
	  if(FaceType::HasVFAdjacencyOcf()) return m.face.IsVFAdjacencyEnabled();
	  else return FaceType::FaceType::HasVFAdjacency();
	}

	template < class VertContainerType, class FaceType >
	  bool HasFFAdjacency (const TriMesh < VertContainerType , face::vector_ocf< FaceType > > & m) 
	{
	  if(FaceType::HasFFAdjacencyOcf()) return m.face.IsFFAdjacencyEnabled();
	  else return FaceType::FaceType::HasFFAdjacency();
	}

	template < class VertContainerType, class FaceType >
      bool HasPerWedgeTexture (const TriMesh < VertContainerType , face::vector_ocf< FaceType > > & m) 
    {
      if(FaceType::HasWedgeTextureOcf()) return m.face.IsWedgeTexEnabled();
      else return FaceType::HasWedgeTexture();
    }

    template < class VertContainerType, class FaceType >
      bool HasPerFaceColor (const TriMesh < VertContainerType , face::vector_ocf< FaceType > > & m) 
    {
      if(FaceType::HasFaceColorOcf()) return m.face.IsColorEnabled();
      else return FaceType::HasFaceColor();
    }

    template < class VertContainerType, class FaceType >
      bool HasPerFaceMark (const TriMesh < VertContainerType , face::vector_ocf< FaceType > > & m) 
    {
      if(FaceType::HasFaceMarkOcf()) return m.face.IsMarkEnabled();
      else return FaceType::HasFaceMark();
    }
  }
}// end namespace vcg
#endif
