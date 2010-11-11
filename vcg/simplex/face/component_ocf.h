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

/* 
Note
OCF = Optional Component Fast (hopefully)
compare with OCC(Optional Component Compact)

Mainly the trick here is to store a base pointer in each simplex...

****************************************************************************/
#ifndef __VCG_FACE_PLUS_COMPONENT_OCF
#define __VCG_FACE_PLUS_COMPONENT_OCF

#include <vcg/simplex/face/component.h>
#include <vector>
#include <limits>


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
	QualityEnabled=false;
	MarkEnabled=false;
  NormalEnabled=false;
  WedgeTexEnabled=false;
  VFAdjacencyEnabled=false;
  FFAdjacencyEnabled=false;
  WedgeColorEnabled=false;
  WedgeNormalEnabled=false;
  }
  
// Auxiliary types to build internal vectors
struct AdjTypePack {
  typename VALUE_TYPE::FacePointer _fp[3] ;    
  char _zp[3] ;  

  // Default constructor. 
  // Needed because we need to know if adjacency is initialized or not
  // when resizing vectors and during an allocate face.
  AdjTypePack() { 
 		_fp[0]=0;
		_fp[1]=0;
		_fp[2]=0;
  }
  };
  
//template <class TexCoordType>
class WedgeTexTypePack {
public:
  WedgeTexTypePack() { 
    wt[0].U()=.5;wt[0].V()=.5;
    wt[1].U()=.5;wt[1].V()=.5;
    wt[2].U()=.5;wt[2].V()=.5;
    wt[0].N()=-1; 
    wt[1].N()=-1; 
    wt[2].N()=-1; 
  }

  typename VALUE_TYPE::TexCoordType wt[3];
};

class WedgeColorTypePack {
public:
  WedgeColorTypePack() {
  typedef typename VALUE_TYPE::ColorType::ScalarType WedgeColorScalarType;
	for (int i=0; i<3; ++i)
	{
    wc[i][0] = WedgeColorScalarType(255);
    wc[i][1] = WedgeColorScalarType(255);
    wc[i][2] = WedgeColorScalarType(255);
    wc[i][3] = WedgeColorScalarType(255);
	}
  }

  typename VALUE_TYPE::ColorType wc[3];
};

class WedgeNormalTypePack {
public:
  WedgeNormalTypePack() {
  typedef typename VALUE_TYPE::NormalType::ScalarType WedgeNormalScalarType;
	for (int i=0; i<3; ++i)
	{
    wn[i][0] = WedgeNormalScalarType(0);
    wn[i][1] = WedgeNormalScalarType(0);
    wn[i][2] = WedgeNormalScalarType(1);
	}
  }

  typename VALUE_TYPE::NormalType wn[3];
};


  // override di tutte le funzioni che possono spostare 
	// l'allocazione in memoria del container
	void push_back(const VALUE_TYPE & v)
  {
    BaseType::push_back(v);
	BaseType::back()._ovp = this;
    if (QualityEnabled)     QV.push_back(0);
	if (ColorEnabled)       CV.push_back(vcg::Color4b(vcg::Color4b::White));
	if (MarkEnabled)        MV.push_back(0);
    if (NormalEnabled)      NV.push_back(typename VALUE_TYPE::NormalType());
	if (VFAdjacencyEnabled) AV.push_back(AdjTypePack());
    if (FFAdjacencyEnabled) AF.push_back(AdjTypePack());
	if (WedgeTexEnabled)    WTV.push_back(WedgeTexTypePack());
	if (WedgeColorEnabled)  WCV.push_back(WedgeColorTypePack());
	if (WedgeNormalEnabled) WNV.push_back(WedgeNormalTypePack());
  }
	void pop_back();
  void resize(const unsigned int & _size) 
  {
	  unsigned int oldsize = BaseType::size();
    BaseType::resize(_size);
	  if(oldsize<_size){
		  ThisTypeIterator firstnew = BaseType::begin();
		  advance(firstnew,oldsize);
		  _updateOVP(firstnew,(*this).end());
	  }  
    if (QualityEnabled)     QV.resize(_size);
	if (ColorEnabled)       CV.resize(_size);
	if (MarkEnabled)        MV.resize(_size);
	if (NormalEnabled)      NV.resize(_size);
	if (VFAdjacencyEnabled) AV.resize(_size);
    if (FFAdjacencyEnabled) AF.resize(_size);
    if (WedgeTexEnabled)    WTV.resize(_size,WedgeTexTypePack());
	if (WedgeColorEnabled)  WCV.resize(_size);
	if (WedgeNormalEnabled) WNV.resize(_size);
   }
  void reserve(const unsigned int & _size)
  {
    BaseType::reserve(_size);

    if (QualityEnabled)     QV.reserve(_size);
	if (ColorEnabled)       CV.reserve(_size);
	if (MarkEnabled)        MV.reserve(_size);
	if (NormalEnabled)      NV.reserve(_size);
	if (VFAdjacencyEnabled) AV.reserve(_size);
    if (FFAdjacencyEnabled) AF.reserve(_size);
	if (WedgeTexEnabled)    WTV.reserve(_size);
	if (WedgeColorEnabled)  WCV.reserve(_size);
	if (WedgeNormalEnabled) WNV.reserve(_size);

	if( BaseType::empty()) return ;

    ThisTypeIterator oldbegin=(*this).begin();
    if(oldbegin!=(*this).begin()) _updateOVP((*this).begin(),(*this).end());
  }

 void _updateOVP(ThisTypeIterator lbegin, ThisTypeIterator lend)
{
    ThisTypeIterator fi;
    //for(fi=(*this).begin();vi!=(*this).end();++vi)
    for(fi=lbegin;fi!=lend;++fi)
        (*fi)._ovp=this;
 }
 
 
  
// this function is called by the specialized Reorder function, that is called whenever someone call the allocator::CompactVertVector
void ReorderFace(std::vector<size_t> &newFaceIndex )
{
	size_t i=0;
	if (QualityEnabled)     assert( QV.size() == newFaceIndex.size() );
	if (ColorEnabled)       assert( CV.size() == newFaceIndex.size() );
	if (MarkEnabled)        assert( MV.size() == newFaceIndex.size() );
	if (NormalEnabled)      assert( NV.size() == newFaceIndex.size() );
	if (VFAdjacencyEnabled) assert( AV.size() == newFaceIndex.size() );
	if (FFAdjacencyEnabled) assert( AF.size() == newFaceIndex.size() );
	if (WedgeTexEnabled)    assert(WTV.size() == newFaceIndex.size() );
	if (WedgeColorEnabled)  assert(WCV.size() == newFaceIndex.size() );
	if (WedgeNormalEnabled) assert(WNV.size() == newFaceIndex.size() );

	for(i=0;i<newFaceIndex.size();++i)
		{
			if(newFaceIndex[i] != std::numeric_limits<size_t>::max() )
				{
					assert(newFaceIndex[i] <= i);
					if (QualityEnabled)      QV[newFaceIndex[i]] =  QV[i]; 
					if (ColorEnabled)        CV[newFaceIndex[i]] =  CV[i]; 
					if (MarkEnabled)         MV[newFaceIndex[i]] =  MV[i];
					if (NormalEnabled)       NV[newFaceIndex[i]] =  NV[i];
					if (VFAdjacencyEnabled)  AV[newFaceIndex[i]] =  AV[i];
					if (FFAdjacencyEnabled)  AF[newFaceIndex[i]] =  AF[i];
					if (WedgeTexEnabled)    WTV[newFaceIndex[i]] = WTV[i];
					if (WedgeColorEnabled)  WCV[newFaceIndex[i]] = WCV[i];
					if (WedgeNormalEnabled) WNV[newFaceIndex[i]] = WNV[i];
				}
		}
	
	if (QualityEnabled)      QV.resize(BaseType::size());
	if (ColorEnabled)        CV.resize(BaseType::size());
	if (MarkEnabled)         MV.resize(BaseType::size());
	if (NormalEnabled)       NV.resize(BaseType::size());
	if (VFAdjacencyEnabled)  AV.resize(BaseType::size());
	if (FFAdjacencyEnabled)  AF.resize(BaseType::size());
	if (WedgeTexEnabled)    WTV.resize(BaseType::size());
	if (WedgeColorEnabled)  WCV.resize(BaseType::size());
	if (WedgeNormalEnabled) WNV.resize(BaseType::size());
}

////////////////////////////////////////
// Enabling Functions
	
bool IsQualityEnabled() const {return QualityEnabled;}
void EnableQuality() {
	assert(VALUE_TYPE::HasFaceQualityOcf());
	QualityEnabled=true;
	QV.resize((*this).size());
}

void DisableQuality() {
	assert(VALUE_TYPE::HasFaceQualityOcf());
	QualityEnabled=false;
	QV.clear();
}
	
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
  assert(VALUE_TYPE::HasWedgeTexCoordOcf());
  WedgeTexEnabled=true;
  WTV.resize((*this).size(),WedgeTexTypePack());
}

void DisableWedgeTex() {
  assert(VALUE_TYPE::HasWedgeTexCoordOcf());
  WedgeTexEnabled=false;
  WTV.clear();
}

bool IsWedgeColorEnabled() const {return WedgeColorEnabled;}
void EnableWedgeColor() {
  assert(VALUE_TYPE::HasWedgeColorOcf());
  WedgeColorEnabled=true;
  WCV.resize((*this).size(),WedgeColorTypePack());
}

void DisableWedgeColor() {
  assert(VALUE_TYPE::HasWedgeColorOcf());
  WedgeColorEnabled=false;
  WCV.clear();
}

bool IsWedgeNormalEnabled() const {return WedgeNormalEnabled;}
void EnableWedgeNormal() {
  assert(VALUE_TYPE::HasWedgeNormalOcf());
  WedgeNormalEnabled=true;
  WNV.resize((*this).size(),WedgeNormalTypePack());
}

void DisableWedgeNormal() {
  assert(VALUE_TYPE::HasWedgeNormalOcf());
  WedgeNormalEnabled=false;
  WNV.clear();
}

public:
  std::vector<float> QV;
  std::vector<typename VALUE_TYPE::ColorType> CV;
  std::vector<int> MV;
  std::vector<typename VALUE_TYPE::NormalType> NV;
  std::vector<struct AdjTypePack> AV;
  std::vector<struct AdjTypePack> AF;
  std::vector<class WedgeTexTypePack> WTV;
  std::vector<class WedgeColorTypePack> WCV;
  std::vector<class WedgeNormalTypePack> WNV;

  bool QualityEnabled;
  bool ColorEnabled;
  bool MarkEnabled;
  bool NormalEnabled;
  bool WedgeTexEnabled;
  bool VFAdjacencyEnabled;
  bool FFAdjacencyEnabled;
  bool WedgeColorEnabled;
  bool WedgeNormalEnabled;
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

	template <class LeftF>
	void ImportData(const LeftF & leftF){
		T::ImportData(leftF);
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
  char cFFi(const int j) const {
    assert((*this).Base().FFAdjacencyEnabled); 
    return (*this).Base().AF[(*this).Index()]._zp[j]; 
  }
	
	typename T::FacePointer        &FFp1( const int j )       { return FFp((j+1)%3);}
	typename T::FacePointer        &FFp2( const int j )       { return FFp((j+2)%3);}
	typename T::FacePointer  const  FFp1( const int j ) const { return FFp((j+1)%3);}
	typename T::FacePointer  const  FFp2( const int j ) const { return FFp((j+2)%3);}

	template <class LeftF>
	void ImportData(const LeftF & leftF){
		T::ImportData(leftF);
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

	template <class LeftF>
	void ImportData(const LeftF & leftF){
		if((*this).Base().NormalEnabled && leftF.Base().NormalEnabled)
			N() = leftF.cN(); 
		T::ImportData(leftF);
	}

};

template <class T> class Normal3sOcf: public NormalOcf<vcg::Point3s, T> {};
template <class T> class Normal3fOcf: public NormalOcf<vcg::Point3f, T> {};
template <class T> class Normal3dOcf: public NormalOcf<vcg::Point3d, T> {};
		
///*-------------------------- QUALITY ----------------------------------*/ 

template <class A, class T> class QualityOcf: public T {
public:
  typedef A QualityType;
  QualityType &Q() { 
    assert((*this).Base().QualityEnabled); 
    return (*this).Base().QV[(*this).Index()]; 
  }
  const QualityType  Q() const  { 
    assert((*this).Base().QualityEnabled); 
    return (*this).Base().QV[(*this).Index()]; 
  }
  const QualityType  cQ() const  { 
    assert((*this).Base().QualityEnabled); 
    return (*this).Base().QV[(*this).Index()]; 
  }

	template <class LeftF>
	void ImportData(const LeftF & leftF){
		//if((*this).Base().QualityEnabled && leftF.Base().QualityEnabled)// WRONG I do not know anything about leftV!
		if((*this).Base().QualityEnabled)
				Q() = leftF.cQ(); 
		T::ImportData(leftF);
	}
  static bool HasFaceQuality()   { return true; }
  static bool HasFaceQualityOcf()   { return true; }
};

template <class T> class QualityfOcf: public QualityOcf<float, T> {};

///*-------------------------- COLOR ----------------------------------*/ 

template <class A, class T> class ColorOcf: public T {
public:
  typedef A ColorType;
  ColorType &C() { 
    assert((*this).Base().ColorEnabled); 
    return (*this).Base().CV[(*this).Index()]; 
  }
	const ColorType  C() const  { 
    assert((*this).Base().ColorEnabled); 
    return (*this).Base().CV[(*this).Index()]; 
  }
  const ColorType  cC() const  { 
    assert((*this).Base().ColorEnabled); 
    return (*this).Base().CV[(*this).Index()]; 
  }

	template <class LeftF>
	void ImportData(const LeftF & leftF){
		//if((*this).Base().ColorEnabled && leftF.Base().ColorEnabled)// WRONG I do not know anything about leftV!
		if((*this).Base().ColorEnabled )
				C() = leftF.cC(); 
		T::ImportData(leftF);
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

	template <class LeftF>
	void ImportData(const LeftF & leftF){
		//if((*this).Base().MarkEnabled && leftF.Base().MarkEnabled)// WRONG I do not know anything about leftV!
		if((*this).Base().MarkEnabled)
			IMark() = leftF.IMark(); 
		T::ImportData(leftF);
	}
  static bool HasFaceMark()   { return true; }
  static bool HasFaceMarkOcf()   { return true; }
  inline void InitIMark()    { IMark() = 0; }
};

///*-------------------------- WEDGE TEXCOORD  ----------------------------------*/ 

template <class A, class TT> class WedgeTexCoordOcf: public TT {
public:
  WedgeTexCoordOcf(){ }
  typedef A TexCoordType;
  TexCoordType &WT(const int i)              { assert((*this).Base().WedgeTexEnabled); return (*this).Base().WTV[(*this).Index()].wt[i]; }
  TexCoordType const &cWT(const int i) const { assert((*this).Base().WedgeTexEnabled); return (*this).Base().WTV[(*this).Index()].wt[i]; }
	template <class LeftF>
	void ImportData(const LeftF & leftF){
		//if(this->Base().WedgeTexEnabled && leftF.Base().WedgeTexEnabled)  // WRONG I do not know anything about leftV!
		if(this->Base().WedgeTexEnabled)
		{ WT(0) = leftF.cWT(0); WT(1) = leftF.cWT(1); WT(2) = leftF.cWT(2); }
		TT::ImportData(leftF);
	}
  static bool HasWedgeTexCoord()   { return true; }
  static bool HasWedgeTexCoordOcf()   { return true; }
};

template <class T> class WedgeTexCoordfOcf: public WedgeTexCoordOcf<TexCoord2<float,1>, T> {};

///*-------------------------- WEDGE COLOR  ----------------------------------*/

template <class A, class TT> class WedgeColorOcf: public TT {
public:
  WedgeColorOcf(){ }
  typedef A ColorType;
  ColorType &WC(const int i)              { assert((*this).Base().WedgeColorEnabled); return (*this).Base().WCV[(*this).Index()].wc[i]; }
  const ColorType cWC(const int i) const { assert((*this).Base().WedgeColorEnabled); return (*this).Base().WCV[(*this).Index()].wc[i]; }
	template <class LeftF>
	void ImportData(const LeftF & leftF){
		//if(this->Base().WedgeColorEnabled && leftF.Base().WedgeColorEnabled)  // WRONG I do not know anything about leftV!
		if(this->Base().WedgeColorEnabled)
		{ WC(0) = leftF.cWC(0); WC(1) = leftF.cWC(1); WC(2) = leftF.cWC(2); }
		TT::ImportData(leftF);
	}
  static bool HasWedgeColor()   { return true; }
  static bool HasWedgeColorOcf()   { return true; }
};

template <class T> class WedgeColor4bOcf: public WedgeColorOcf<vcg::Color4b, T> {};

///*-------------------------- WEDGE NORMAL ----------------------------------*/

template <class A, class TT> class WedgeNormalOcf: public TT {
public:
  WedgeNormalOcf(){ }
  typedef A NormalType;
  NormalType &WN(const int i)              { assert((*this).Base().WedgeNormalEnabled); return (*this).Base().WNV[(*this).Index()].wn[i]; }
  NormalType const &cWN(const int i) const { assert((*this).Base().WedgeNormalEnabled); return (*this).Base().WNV[(*this).Index()].wn[i]; }
	template <class LeftF>
	void ImportData(const LeftF & leftF){
		//if(this->Base().WedgeNormalEnabled && leftF.Base().WedgeNormalEnabled)  // WRONG I do not know anything about leftV!
		if(this->Base().WedgeNormalEnabled)
		{ WN(0) = leftF.cWN(0); WN(1) = leftF.cWN(1); WN(2) = leftF.cWN(2); }
		TT::ImportData(leftF);
	}
  static bool HasWedgeNormal()   { return true; }
  static bool HasWedgeNormalOcf()   { return true; }
};

template <class T> class WedgeNormal3sOcf: public WedgeNormalOcf<vcg::Point3s, T> {};
template <class T> class WedgeNormal3fOcf: public WedgeNormalOcf<vcg::Point3f, T> {};
template <class T> class WedgeNormal3dOcf: public WedgeNormalOcf<vcg::Point3d, T> {};

///*-------------------------- InfoOpt  ----------------------------------*/

template < class T> class InfoOcf: public T {
public:
    // You should never ever try to copy a vertex that has OCF stuff.
		// use ImportData function.
    inline InfoOcf &operator=(const InfoOcf & /*other*/) {
        assert(0); return *this;
    }


  vector_ocf<typename T::FaceType> &Base() const { return *_ovp;}

	template <class LeftF>
	void ImportData(const LeftF & leftF){T::ImportData(leftF);}

  static bool HasFaceColorOcf()   { return false; }
  static bool HasFaceNormalOcf()   { return false; }
  static bool HasFaceMarkOcf()   { return false; }
  static bool HasWedgeTexCoordOcf()   { return false; }
  static bool HasFFAdjacencyOcf()   { return false; }
  static bool HasVFAdjacencyOcf()   { return false; }
  //static bool HasFaceQualityOcf()   { return false; }

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
  
  template < class, class,class,class > class TriMesh;

  namespace tri
  {


	template < class VertContainerType, class FaceType, class Container1, class Container2  >
		bool HasPerFaceVFAdjacency (const TriMesh < VertContainerType , face::vector_ocf< FaceType >, Container1, Container2 > & m)
	{
	  if(FaceType::HasVFAdjacencyOcf()) return m.face.IsVFAdjacencyEnabled();
	  else return FaceType::FaceType::HasVFAdjacency();
	}

	template < class VertContainerType, class FaceType, class Container1, class Container2 >
		bool HasFFAdjacency (const TriMesh < VertContainerType , face::vector_ocf< FaceType >, Container1, Container2  > & m)
	{
	  if(FaceType::HasFFAdjacencyOcf()) return m.face.IsFFAdjacencyEnabled();
	  else return FaceType::FaceType::HasFFAdjacency();
	}

	template < class VertContainerType, class FaceType, class Container1, class Container2  >
			bool HasPerWedgeTexCoord (const TriMesh < VertContainerType , face::vector_ocf< FaceType >, Container1, Container2  > & m)
    {
      if(FaceType::HasWedgeTexCoordOcf()) return m.face.IsWedgeTexEnabled();
      else return FaceType::HasWedgeTexCoord();
    }

		template < class VertContainerType, class FaceType, class Container1, class Container2  >
			bool HasPerFaceColor (const TriMesh < VertContainerType , face::vector_ocf< FaceType >, Container1, Container2  > & m)
    {
      if(FaceType::HasFaceColorOcf()) return m.face.IsColorEnabled();
      else return FaceType::HasFaceColor();
    }
		
		template < class VertContainerType, class FaceType, class Container1, class Container2  >
		bool HasPerFaceQuality (const TriMesh < VertContainerType , face::vector_ocf< FaceType >, Container1, Container2  > & m)
    {
      if(FaceType::HasFaceQualityOcf()) return m.face.IsQualityEnabled();
      else return FaceType::HasFaceQuality();
    }
		
		template < class VertContainerType, class FaceType, class Container1, class Container2  >
			bool HasPerFaceMark (const TriMesh < VertContainerType , face::vector_ocf< FaceType >, Container1, Container2  > & m)
    {
      if(FaceType::HasFaceMarkOcf()) return m.face.IsMarkEnabled();
      else return FaceType::HasFaceMark();
    }

		template < class FaceType >
			void ReorderFace( std::vector<size_t>  &newFaceIndex, face::vector_ocf< FaceType > &faceVec)
		{
			faceVec.ReorderFace(newFaceIndex);
		}

  }
}// end namespace vcg
#endif
