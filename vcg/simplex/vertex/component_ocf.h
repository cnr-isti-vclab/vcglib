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
Revision 1.16  2008/04/03 23:15:40  cignoni
added optional mark and cleaned up some nasty type bugs.

Revision 1.15  2008/03/17 11:39:15  ganovelli
added curvature and curvatruredir (compiled .net 2005 and gcc)

Revision 1.14  2008/03/11 09:22:07  cignoni
Completed the garbage collecting functions CompactVertexVector and CompactFaceVector.

Revision 1.13  2008/02/05 20:42:43  cignoni
Other small typos

Revision 1.12  2008/02/04 21:26:49  ganovelli
added ImportLocal which imports all local attributes into vertexplus and faceplus.
A local attribute is everything (N(), C(), Q()....) except pointers to other simplices
(i.e. FFAdj, VFAdj, VertexRef) which are set to NULL.
Added some function for const attributes

Revision 1.11  2007/12/11 18:25:31  cignoni
added missing include limits

Revision 1.10  2007/12/11 11:36:03  cignoni
Added the CompactVertexVector garbage collecting function.

Revision 1.9  2006/12/11 23:42:00  ganovelli
bug Index()() instead of Index()

Revision 1.8  2006/12/04 11:17:42  ganovelli
added forward declaration of TriMesh

Revision 1.7  2006/11/07 17:22:52  cignoni
many gcc compiling issues

Revision 1.6  2006/11/07 15:13:57  zifnab1974
Necessary changes for compilation with gcc 3.4.6. Especially the hash function is a problem

Revision 1.5  2006/11/07 11:29:24  cignoni
Corrected some errors in the reflections Has*** functions

Revision 1.4  2006/10/31 16:02:59  ganovelli
vesione 2005 compliant

Revision 1.3  2006/02/28 11:59:55  ponchio
g++ compliance:

begin() -> (*this).begin() and for end(), size(), Base(), Index()

Revision 1.2  2005/12/12 11:17:32  cignoni
Corrected update function, now only the needed simplexes should be updated.

Revision 1.1  2005/10/14 15:07:59  cignoni
First Really Working version


****************************************************************************/

/*
Note
OCF = Optional Component Fast (hopefully)
compare with OCC(Optional Component Compact)

Mainly the trick here is to store a base pointer in each simplex...

****************************************************************************/
#ifndef __VCG_VERTEX_PLUS_COMPONENT_OCF
#define __VCG_VERTEX_PLUS_COMPONENT_OCF

#include <vcg/simplex/vertex/component.h>

#include <vector>
#include <limits>

namespace vcg {
	namespace vertex {
/*
All the Components that can be added to a vertex should be defined in the namespace vert:

*/
template <class VALUE_TYPE>
class vector_ocf: public std::vector<VALUE_TYPE> {
	typedef std::vector<VALUE_TYPE> BaseType;
	typedef typename vector_ocf<VALUE_TYPE>::iterator ThisTypeIterator;

public:
	vector_ocf():std::vector<VALUE_TYPE>(){
		ColorEnabled = false;
    CurvatureEnabled = false;
    CurvatureDirEnabled = false;
    MarkEnabled = false;
		NormalEnabled = false;
    QualityEnabled = false;
    RadiusEnabled = false;
    TexCoordEnabled = false;
		VFAdjacencyEnabled = false;
	}

	// override di tutte le funzioni che possono spostare
	// l'allocazione in memoria del container
	void push_back(const VALUE_TYPE & v)
	{
		BaseType::push_back(v);
		BaseType::back()._ovp = this;
		if (ColorEnabled)         CV.push_back(vcg::Color4b(vcg::Color4b::White));
		if (MarkEnabled)          MV.push_back(0);
		if (NormalEnabled)        NV.push_back(typename VALUE_TYPE::NormalType());
		if (TexCoordEnabled)      TV.push_back(typename VALUE_TYPE::TexCoordType());
		if (VFAdjacencyEnabled)   AV.push_back(VFAdjType());
		if (CurvatureEnabled)     CuV.push_back(typename VALUE_TYPE::CurvatureType());
		if (CurvatureDirEnabled)  CuDV.push_back(typename VALUE_TYPE::CurvatureDirType());
		if (RadiusEnabled)        RadiusV.push_back(typename VALUE_TYPE::RadiusType());
	}
	void pop_back();
	void resize(const unsigned int & _size)
	{
		const unsigned int oldsize = BaseType::size();
			BaseType::resize(_size);
		if(oldsize<_size){
			ThisTypeIterator firstnew = BaseType::begin();
			advance(firstnew,oldsize);
			_updateOVP(firstnew,(*this).end());
		}
		if (ColorEnabled)         CV.resize(_size);
		if (MarkEnabled)          MV.resize(_size);
		if (NormalEnabled)        NV.resize(_size);
		if (TexCoordEnabled)      TV.resize(_size);
		if (VFAdjacencyEnabled)   AV.resize(_size);
		if (CurvatureEnabled)     CuV.resize(_size);
		if (CurvatureDirEnabled)  CuDV.resize(_size);
		if (RadiusEnabled)        RadiusV.resize(_size);
	}

	void reserve(const unsigned int & _size)
	{
		BaseType::reserve(_size);
		if (ColorEnabled)        CV.reserve(_size);
		if (MarkEnabled)         MV.reserve(_size);
		if (NormalEnabled)       NV.reserve(_size);
		if (TexCoordEnabled)     TV.reserve(_size);
		if (VFAdjacencyEnabled)  AV.reserve(_size);
		if (CurvatureEnabled)    CuV.reserve(_size);
		if (CurvatureDirEnabled) CuDV.reserve(_size);
		if (RadiusEnabled)       RadiusV.reserve(_size);
	}

	void _updateOVP(ThisTypeIterator lbegin, ThisTypeIterator lend)
	{
		ThisTypeIterator vi;
		for(vi=lbegin;vi!=lend;++vi)
				(*vi)._ovp=this;
	}

// this function is called by the specialized Reorder function, that is called whenever someone call the allocator::CompactVertVector
void ReorderVert(std::vector<size_t> &newVertIndex )
{
	size_t i=0;
	assert( (!ColorEnabled)         || ( CV.size() == newVertIndex.size() ) );
	assert( (!MarkEnabled)          || ( MV.size() == newVertIndex.size() ) );
	assert( (!NormalEnabled)        || ( NV.size() == newVertIndex.size() ) );
	assert( (!VFAdjacencyEnabled)   || ( AV.size() == newVertIndex.size() ) );
	assert( (!CurvatureEnabled)     || ( CuV.size() == newVertIndex.size() ) );
	assert( (!CurvatureDirEnabled)  || ( CuDV.size() == newVertIndex.size() ) );
	assert( (!RadiusEnabled)        || ( RadiusV.size() == newVertIndex.size() ) );
	assert( (!TexCoordEnabled)       || ( TV.size() == newVertIndex.size() ) );

	for(i=0;i<newVertIndex.size();++i)
		{
			if(newVertIndex[i] != std::numeric_limits<size_t>::max() )
				{
					assert(newVertIndex[i] <= i);
					if (ColorEnabled)			    CV[newVertIndex[i]] =  CV[i];
					if (MarkEnabled)				  MV[newVertIndex[i]] =  MV[i];
					if (NormalEnabled)        NV[newVertIndex[i]] =  NV[i];
					if (VFAdjacencyEnabled)   AV[newVertIndex[i]] =  AV[i];
					if (CurvatureEnabled)     CuV[newVertIndex[i]] = CuV[i];
					if (CurvatureDirEnabled)  CuDV[newVertIndex[i]] =CuDV[i];
					if (RadiusEnabled)        RadiusV[newVertIndex[i]] = RadiusV[i];
					if (TexCoordEnabled)        TV[newVertIndex[i]] = TV[i];
				}
		}

	if (ColorEnabled)         CV.resize(BaseType::size());
	if (MarkEnabled)          MV.resize(BaseType::size());
	if (NormalEnabled)        NV.resize(BaseType::size());
	if (VFAdjacencyEnabled)   AV.resize(BaseType::size());
	if (CurvatureEnabled)     CuV.resize(BaseType::size());
	if (CurvatureDirEnabled)  CuDV.resize(BaseType::size());
	if (RadiusEnabled)  RadiusV.resize(BaseType::size());
	if (TexCoordEnabled)  TV.resize(BaseType::size());
}



////////////////////////////////////////
// Enabling Eunctions

bool IsQualityEnabled() const {return QualityEnabled;}
void EnableQuality() {
	assert(VALUE_TYPE::HasQualityOcf());
	QualityEnabled=true;
	QV.resize((*this).size());
}

void DisableQuality() {
	assert(VALUE_TYPE::HasQualityOcf());
	QualityEnabled=false;
	QV.clear();
}

bool IsColorEnabled() const {return ColorEnabled;}
void EnableColor() {
	assert(VALUE_TYPE::HasColorOcf());
	ColorEnabled=true;
	CV.resize((*this).size());
}

void DisableColor() {
	assert(VALUE_TYPE::HasColorOcf());
	ColorEnabled=false;
	CV.clear();
}

bool IsMarkEnabled() const {return MarkEnabled;}
void EnableMark() {
	assert(VALUE_TYPE::HasMarkOcf());
	MarkEnabled=true;
	MV.resize((*this).size());
}

void DisableMark() {
	assert(VALUE_TYPE::HasMarkOcf());
	MarkEnabled=false;
	MV.clear();
}

bool IsNormalEnabled() const {return NormalEnabled;}
void EnableNormal() {
	assert(VALUE_TYPE::HasNormalOcf());
	NormalEnabled=true;
	NV.resize((*this).size());
}

void DisableNormal() {
	assert(VALUE_TYPE::HasNormalOcf());
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

bool IsCurvatureEnabled() const {return CurvatureEnabled;}
void EnableCurvature() {
	assert(VALUE_TYPE::HasCurvatureOcf());
	CurvatureEnabled=true;
	CuV.resize((*this).size());
}

void DisableCurvature() {
	assert(VALUE_TYPE::HasCurvatureOcf());
	CurvatureEnabled=false;
	CuV.clear();
}

bool IsCurvatureDirEnabled() const {return CurvatureDirEnabled;}
void EnableCurvatureDir() {
	assert(VALUE_TYPE::HasCurvatureDirOcf());
	CurvatureDirEnabled=true;
	CuDV.resize((*this).size());
}

void DisableCurvatureDir() {
	assert(VALUE_TYPE::HasCurvatureDirOcf());
	CurvatureDirEnabled=false;
	CuDV.clear();
}

bool IsRadiusEnabled() const {return RadiusEnabled;}
void EnableRadius() {
	assert(VALUE_TYPE::HasRadiusOcf());
	RadiusEnabled=true;
	RadiusV.resize((*this).size());
}

void DisableRadius() {
	assert(VALUE_TYPE::HasRadiusOcf());
	RadiusEnabled=false;
	RadiusV.clear();
}


bool IsTexCoordEnabled() const {return TexCoordEnabled;}
void EnableTexCoord() {
	assert(VALUE_TYPE::HasTexCoordOcf());
	TexCoordEnabled=true;
	TV.resize((*this).size());
}

void DisableTexCoord() {
	assert(VALUE_TYPE::HasTexCoordOcf());
	TexCoordEnabled=false;
	TV.clear();
}
struct VFAdjType {
	typename VALUE_TYPE::FacePointer _fp ;
	int _zp ;
	};

public:
  std::vector<typename VALUE_TYPE::ColorType> CV;
  std::vector<typename VALUE_TYPE::CurvatureType> CuV;
	std::vector<typename VALUE_TYPE::CurvatureDirType> CuDV;
  std::vector<int> MV;
  std::vector<typename VALUE_TYPE::NormalType> NV;
  std::vector<typename VALUE_TYPE::QualityType> QV;
	std::vector<typename VALUE_TYPE::RadiusType> RadiusV;
	std::vector<typename VALUE_TYPE::TexCoordType> TV;
	std::vector<struct VFAdjType> AV;

  bool ColorEnabled;
	bool CurvatureEnabled;
	bool CurvatureDirEnabled;
  bool MarkEnabled;
  bool NormalEnabled;
  bool QualityEnabled;
	bool RadiusEnabled;
  bool TexCoordEnabled;
  bool VFAdjacencyEnabled;
};


//template<>	void EnableAttribute<typename VALUE_TYPE::NormalType>(){	NormalEnabled=true;}

/*------------------------- COORD -----------------------------------------*/
/*----------------------------- VFADJ ------------------------------*/


template <class T> class VFAdjOcf: public T {
public:
	typename T::FacePointer &VFp() {
		assert((*this).Base().VFAdjacencyEnabled);
		return (*this).Base().AV[(*this).Index()]._fp;
	}

	typename T::FacePointer cVFp() const {
		if(! (*this).Base().VFAdjacencyEnabled ) return 0;
		else return (*this).Base().AV[(*this).Index()]._fp;
	}

    int &VFi() {
        assert((*this).Base().VFAdjacencyEnabled);
        return (*this).Base().AV[(*this).Index()]._zp;
    }
    int cVFi() const {
        if(! (*this).Base().VFAdjacencyEnabled ) return -1;
        return (*this).Base().AV[(*this).Index()]._zp;
    }
    template <class LeftV>
	void ImportLocal(const LeftV & leftV)
	{
		if((*this).Base().VFAdjacencyEnabled) // init the data only if they are enabled!
			{
				VFp() = NULL;
				VFi() = -1;
			}
		T::ImportLocal(leftV);
	}

  static bool HasVFAdjacency()   {   return true; }
  static bool HasVFAdjacencyOcf()   {   return true; }
	static bool IsVFAdjacencyEnabled(const typename T::VertexType *vp)   {return vp->Base().VFAdjacencyEnabled;}

private:
};

/*------------------------- Normal -----------------------------------------*/

template <class A, class T> class NormalOcf: public T {
public:
	typedef A NormalType;
	static bool HasNormal()   { return true; }
	static bool HasNormalOcf()   { return true; }

	NormalType &N() {
		// you cannot use Normals before enabling them with: yourmesh.vert.EnableNormal()
		assert((*this).Base().NormalEnabled);
		return (*this).Base().NV[(*this).Index()];  }
	const NormalType &N() const {
		// you cannot use Normals before enabling them with: yourmesh.vert.EnableNormal()
		assert((*this).Base().NormalEnabled);
		return (*this).Base().NV[(*this).Index()];  }

	const NormalType &cN() const {
		// you cannot use Normals before enabling them with: yourmesh.vert.EnableNormal()
		assert((*this).Base().NormalEnabled);
		return (*this).Base().NV[(*this).Index()];  }

	template <class LeftV>
	void ImportLocal(const LeftV & leftV){
			if((*this).Base().NormalEnabled && leftV.Base().NormalEnabled ) // copy the data only if they are enabled in both vertices
				N().Import(leftV.cN());
			T::ImportLocal(leftV);}
};

template <class T> class Normal3sOcf: public NormalOcf<vcg::Point3s, T> {};
template <class T> class Normal3fOcf: public NormalOcf<vcg::Point3f, T> {};
template <class T> class Normal3dOcf: public NormalOcf<vcg::Point3d, T> {};

///*-------------------------- COLOR ----------------------------------*/

template <class A, class T> class ColorOcf: public T {
public:
	typedef A ColorType;
	ColorType &C() { assert((*this).Base().ColorEnabled); return (*this).Base().CV[(*this).Index()]; }
	const ColorType &cC() const { assert((*this).Base().ColorEnabled); return (*this).Base().CV[(*this).Index()]; }
	template <class LeftV>
	void ImportLocal(const LeftV & leftV)
		{
			if((*this).Base().ColorEnabled && leftV.Base().ColorEnabled ) // copy the data only if they are enabled in both vertices
					C() = leftV.cC();
			T::ImportLocal(leftV);
		}

	static bool HasColor()   { return true; }
	static bool HasColorOcf()   { assert(!T::HasColorOcf()); return true; }
};

template <class T> class Color4bOcf: public ColorOcf<vcg::Color4b, T> {};

///*-------------------------- QUALITY  ----------------------------------*/

template <class A, class T> class QualityOcf: public T {
public:
	typedef A QualityType;
	QualityType &Q() { assert((*this).Base().QualityEnabled); return (*this).Base().QV[(*this).Index()]; }
	template <class LeftV>
	void ImportLocal(const LeftV & leftV)
		{
//			if((*this).Base().QualityEnabled && leftV.Base().QualityEnabled ) // copy the data only if they are enabled in both vertices
      if((*this).Base().QualityEnabled && leftV.HasQuality() ) // copy the data only if they are enabled in both vertices
						Q() = leftV.cQ();
			T::ImportLocal(leftV);
		}
	static bool HasQuality()   { return true; }
	static bool HasQualityOcf()   { assert(!T::HasQualityOcf()); return true; }
};

template <class T> class QualityfOcf: public QualityOcf<float, T> {};


///*-------------------------- TEXTURE  ----------------------------------*/

template <class A, class TT> class TexCoordOcf: public TT {
public:
  typedef A TexCoordType;
  TexCoordType &T() {  assert((*this).Base().TexCoordEnabled); return (*this).Base().TV[(*this).Index()]; }
  const TexCoordType &cT() const {  assert((*this).Base().TexCoordEnabled); return (*this).Base().TV[(*this).Index()]; }
	template < class LeftV>
	void ImportLocal(const LeftV & leftV)
		{
			//if((*this).Base().TexCoordEnabled && leftV.Base().TexCoordEnabled ) // WRONG I do not know anything about leftV!
		if((*this).Base().TexCoordEnabled) // copy the data only if they are enabled in both vertices
						T() = leftV.cT();
			TT::ImportLocal(leftV);
		}
	static bool HasTexCoord()   { return true; }
	static bool HasTexCoordOcf()   { assert(!TT::HasTexCoordOcf()); return true; }
};

template <class T> class TexCoordfOcf: public TexCoordOcf<TexCoord2<float,1>, T> {};

///*-------------------------- MARK  ----------------------------------*/

template <class T> class MarkOcf: public T {
public:
    typedef int MarkType;
	inline int & IMark()       {
		assert((*this).Base().MarkEnabled);
		return (*this).Base().MV[(*this).Index()];
	}

	inline int IMark() const   {
		assert((*this).Base().MarkEnabled);
		return (*this).Base().MV[(*this).Index()];
	}

	template <class LeftV>
	void ImportLocal(const LeftV & leftV)
	{
		//if((*this).Base().MarkEnabled && leftV.Base().MarkEnabled ) // WRONG I do not know anything about leftV!
		if((*this).Base().MarkEnabled) // copy the data only if they are enabled in both vertices
				IMark() = leftV.IMark();
		T::ImportLocal(leftV);
	}
	static bool HasMark()   { return true; }
	static bool HasMarkOcf()   { return true; }
	inline void InitIMark()    { IMark() = 0; }
};


///*-------------------------- CURVATURE  ----------------------------------*/

template <class A, class TT> class CurvatureOcf: public TT {
public:
	typedef Point2<A> CurvatureType;
	typedef typename CurvatureType::ScalarType ScalarType;

	ScalarType  &Kh(){  assert((*this).Base().CurvatureEnabled); return (*this).Base().CuV[(*this).Index()][0];}
	ScalarType  &Kg(){  assert((*this).Base().CurvatureEnabled); return (*this).Base().CuV[(*this).Index()][1];}
	const ScalarType &cKh() const { assert((*this).Base().CurvatureEnabled); return (*this).Base().CuV[(*this).Index()][0];}
	const ScalarType &cKg() const { assert((*this).Base().CurvatureEnabled); return (*this).Base().CuV[(*this).Index()][1];}

	template <class LeftV>
	void ImportLocal(const LeftV & leftV){
//		if((*this).Base().CurvatureEnabled && leftV.Base().CurvatureEnabled ) // WRONG I do not know anything about leftV!
      if((*this).Base().CurvatureEnabled && LeftV::IsCurvatureEnabled(&leftV))
			{
				(*this).Base().CuV[(*this).Index()][0] = leftV.cKh();
				(*this).Base().CuV[(*this).Index()][1] = leftV.cKg();
			}
		TT::ImportLocal(leftV);
	}

  static bool HasCurvature() { return true; }
	static bool IsCurvatureEnabled(const typename TT::VertexType *v)   { return v->Base().CurvatureEnabled; }

	static bool HasCurvatureOcf()   { return true; }
	static void Name(std::vector<std::string> & name){name.push_back(std::string("CurvatureOcf"));TT::Name(name);}

private:
};

template <class T> class CurvaturefOcf: public CurvatureOcf<float, T> {};
template <class T> class CurvaturedOcf: public CurvatureOcf<double, T> {};


///*-------------------------- CURVATURE DIR ----------------------------------*/

template <class S>
struct CurvatureDirTypeOcf{
	typedef Point3<S> VecType;
	typedef  S   ScalarType;
	CurvatureDirTypeOcf () {}
	Point3<S>max_dir,min_dir;
	S k1,k2;
};


template <class A, class TT> class CurvatureDirOcf: public TT {
public:
	typedef A CurvatureDirType;
	typedef typename CurvatureDirType::VecType VecType;
	typedef typename CurvatureDirType::ScalarType ScalarType;

	VecType &PD1(){ assert((*this).Base().CurvatureDirEnabled); return (*this).Base().CuDV[(*this).Index()].max_dir;}
	VecType &PD2(){ assert((*this).Base().CurvatureDirEnabled); return (*this).Base().CuDV[(*this).Index()].min_dir;}
  const VecType &cPD1() const {assert((*this).Base().CurvatureDirEnabled); return (*this).Base().CuDV[(*this).Index()].max_dir;}
  const VecType &cPD2() const {assert((*this).Base().CurvatureDirEnabled); return (*this).Base().CuDV[(*this).Index()].min_dir;}

	ScalarType &K1(){ assert((*this).Base().CurvatureDirEnabled); return (*this).Base().CuDV[(*this).Index()].k1;}
	ScalarType &K2(){ assert((*this).Base().CurvatureDirEnabled); return (*this).Base().CuDV[(*this).Index()].k2;}
	const ScalarType &cK1() const {assert((*this).Base().CurvatureDirEnabled); return (*this).Base().CuDV[(*this).Index()].k1;}
	const ScalarType &cK2()const  {assert((*this).Base().CurvatureDirEnabled); return (*this).Base().CuDV[(*this).Index()].k2;}

  template <class LeftV>
  void ImportLocal(const LeftV & leftV){
//		if((*this).Base().CurvatureEnabled && leftV.Base().CurvatureEnabled ) // WRONG I do not know anything about leftV!
    if((*this).Base().CurvatureDirEnabled && LeftV::IsCurvatureDirEnabled(&leftV))
      {
        (*this).PD1() = leftV.cPD1();
        (*this).PD2() = leftV.cPD2();
        (*this).K1() = leftV.cK1();
        (*this).K2() = leftV.cK2();
      }
    TT::ImportLocal(leftV);
  }
  static bool HasCurvatureDir()  { return true; }
  static bool HasCurvatureDirOcf()  { return true; }
  static void Name(std::vector<std::string> & name){name.push_back(std::string("CurvatureDirOcf"));TT::Name(name);}

private:
};


template <class T> class CurvatureDirfOcf: public CurvatureDirOcf<CurvatureDirTypeOcf<float>, T> {
public:	static void Name(std::vector<std::string> & name){name.push_back(std::string("CurvatureDirfOcf"));T::Name(name);}
};
template <class T> class CurvatureDirdOcf: public CurvatureDirOcf<CurvatureDirTypeOcf<double>, T> {
public:	static void Name(std::vector<std::string> & name){name.push_back(std::string("CurvatureDirdOcf"));T::Name(name);}
};


///*-------------------------- RADIUS  ----------------------------------*/

template <class A, class TT> class RadiusOcf: public TT {
public:
	typedef A RadiusType;
	typedef RadiusType ScalarType;

	RadiusType &R(){  assert((*this).Base().RadiusEnabled); return (*this).Base().RadiusV[(*this).Index()];}
	const RadiusType &cR() const { assert((*this).Base().RadiusEnabled); return (*this).Base().RadiusV[(*this).Index()];}

	template <class LeftV>
	void ImportLocal(const LeftV & leftV)
	{
		//if ((*this).Base().RadiusEnabled && leftV.Base().RadiusEnabled ) // WRONG I do not know anything about leftV!
		if ((*this).Base().RadiusEnabled) 
			(*this).Base().RadiusV[(*this).Index()] = leftV.cR();
		TT::ImportLocal(leftV);
	}

	static bool HasRadius()     { return true; }
	static bool HasRadiusOcf()  { return true; }
	static void Name(std::vector<std::string> & name){name.push_back(std::string("RadiusOcf")); TT::Name(name);}

private:
};

template <class T> class RadiusfOcf: public RadiusOcf<float, T> {};
template <class T> class RadiusdOcf: public RadiusOcf<double, T> {};


///*-------------------------- InfoOpt  ----------------------------------*/

template < class T> class InfoOcf: public T {
public:
    // You should never ever try to copy a vertex that has OCF stuff.
    // use ImportLocal function.
    inline InfoOcf &operator=(const InfoOcf & /*other*/) {
        assert(0); return *this;
    }

		vector_ocf<typename T::VertexType> &Base() const { return *_ovp;}

	inline int Index() const {
		typename  T::VertexType const *tp=static_cast<typename T::VertexType const*>(this);
		int tt2=tp- &*(_ovp->begin());
		return tt2;
	}
public:
	vector_ocf<typename T::VertexType> *_ovp;

	static bool HasQualityOcf()   { return false; }
	static bool HasTexCoordOcf()   { return false; }
	static bool HasVFAdjacencyOcf()   { return false; }
};


} // end namespace vert


namespace tri
{

	template < class, class,class, class> class TriMesh;

	template < class VertexType, class ContainerType0, class ContainerType1 ,class ContainerType2 >
	bool HasPerVertexRadius (const TriMesh < vertex::vector_ocf< VertexType > , ContainerType0,   ContainerType1,  ContainerType2 > & m)
	{
		if(VertexType::HasRadiusOcf()) return m.vert.IsRadiusEnabled();
		else return VertexType::HasRadius();
	}
	
	template < class VertexType, class ContainerType0, class ContainerType1 ,class ContainerType2>
		bool HasPerVertexQuality (const TriMesh < vertex::vector_ocf< VertexType > , ContainerType0,   ContainerType1,  ContainerType2> & m)
	{
		if(VertexType::HasQualityOcf()) return m.vert.IsQualityEnabled();
		else return VertexType::HasQuality();
	}
	
	template < class VertexType, class ContainerType0, class ContainerType1 ,class ContainerType2 >
		bool HasPerVertexTexCoord (const TriMesh < vertex::vector_ocf< VertexType > , ContainerType0,   ContainerType1,  ContainerType2> & m)
	{
		if(VertexType::HasTexCoordOcf()) return m.vert.IsTexCoordEnabled();
		else return VertexType::HasTexCoord();
	}

	template < class VertexType, class ContainerType0, class ContainerType1 ,class ContainerType2 >
		bool HasPerVertexNormal (const TriMesh < vertex::vector_ocf< VertexType > , ContainerType0,   ContainerType1,  ContainerType2> & m)
	{
		if(VertexType::HasNormalOcf()) return m.vert.IsNormalEnabled();
		else return VertexType::HasNormal();
	}

	template < class VertexType, class ContainerType0, class ContainerType1 ,class ContainerType2>
		bool HasPerVertexColor (const TriMesh < vertex::vector_ocf< VertexType > , ContainerType0,   ContainerType1,  ContainerType2> & m)
	{
		if(VertexType::HasColorOcf()) return m.vert.IsColorEnabled();
		else return VertexType::HasColor();
	}

	template < class VertexType, class ContainerType0, class ContainerType1 ,class ContainerType2 >
		bool HasPerVertexCurvature (const TriMesh < vertex::vector_ocf< VertexType > , ContainerType0,   ContainerType1,  ContainerType2 > & m)
	{
		if(VertexType::HasCurvatureOcf()) return m.vert.IsCurvatureEnabled();
		else return VertexType::HasCurvature();
	}

	template < class VertexType, class ContainerType0, class ContainerType1 ,class ContainerType2>
		bool HasPerVertexCurvatureDir (const TriMesh < vertex::vector_ocf< VertexType > , ContainerType0,   ContainerType1,  ContainerType2 > & m)
	{
		if(VertexType::HasCurvatureDirOcf()) return m.vert.IsCurvatureDirEnabled();
		else return VertexType::HasCurvatureDir();
	}

	template < class VertexType >
	void ReorderVert( std::vector<size_t>  &newVertIndex, vertex::vector_ocf< VertexType > &vertVec)
		{
			vertVec.ReorderVert(newVertIndex);
		}
}
}// end namespace vcg
#endif
