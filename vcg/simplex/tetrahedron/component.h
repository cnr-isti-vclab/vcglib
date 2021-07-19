/****************************************************************************
* VCGLib                                                            o o     *
* Visual and Computer Graphics Library                            o     o   *
*                                                                _   O  _   *
* Copyright(C) 2004-2016                                           \/)\/    *
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
Revision 1.1  2007/05/09 10:31:53  ganovelli
added



****************************************************************************/
#ifndef __VCG_TETRAHEDRON_PLUS_COMPONENT
#define __VCG_TETRAHEDRON_PLUS_COMPONENT

#include <vector>

#include <vcg/space/color4.h>
#include <vcg/space/tetra3.h>

namespace vcg {
namespace tetrahedron {
/*
Some naming Rules
All the Components that can be added to a tetra should be defined in the namespace tetra:
*/
template <class T> class EmptyCore : public T {
public:
    //Empty vertexref
    inline typename T::VertexType *       & V( const int  )      {	assert(0);		static typename T::VertexType *vp=0; return vp; }
    inline typename T::VertexType * const & V( const int ) const {	assert(0);		static typename T::VertexType *vp=0; return vp; }
    inline const typename T::VertexType *  cV( const int ) const {	assert(0);		static typename T::VertexType *vp=0; return vp;	}
    inline       typename T::CoordType & P( const int )          {	assert(0);		static typename T::CoordType coord(0, 0, 0); return coord;	}
    inline const typename T::CoordType & P( const int )    const {	assert(0);		static typename T::CoordType coord(0, 0, 0); return coord;	}
    inline const typename T::CoordType &cP( const int )    const {	assert(0);		static typename T::CoordType coord(0, 0, 0); return coord;	}

    static bool HasVertexRef()   { return false; }
    static bool HasTVAdjacency() { return false; }

    // //Empty normals
    // typedef typename T::VertexType::NormalType NormalType;
    // NormalType &N(const int & ){ static NormalType dummynormal(0, 0, 0); assert(0); return dummynormal; }
    // const NormalType cN(const int & ) const { static NormalType dummynormal(0, 0, 0); assert(0); return dummynormal; }

    // static bool HasFaceNormal()      { return false; }
    // static bool HasFaceNormalOcc()   { return false; }

    //Empty color
    typedef vcg::Color4b ColorType;
    ColorType &C()       { static ColorType dummycolor(vcg::Color4b::White); assert(0); return dummycolor; }
    ColorType cC() const { static ColorType dummycolor(vcg::Color4b::White); assert(0); return dummycolor; }

    static bool HasColor() { return false; }
    static bool IsColorEnabled() { return T::TetraType::HasColor(); }

    //Empty Quality
    typedef float QualityType;
    typedef vcg::Point3f Quality3Type;
    QualityType &Q()       { static QualityType dummyquality(0); assert(0); return dummyquality; }
    QualityType cQ() const { static QualityType dummyquality(0); assert(0); return dummyquality; }
    Quality3Type &Q3()       { static Quality3Type dummyQuality3(0,0,0);  assert(0); return dummyQuality3; }
    Quality3Type cQ3() const { static Quality3Type dummyQuality3(0,0,0);  assert(0); return dummyQuality3; }

    static bool HasQuality() { return false; }
    static bool HasQuality3() { return false; }
    inline bool IsQualityEnabled() const { return T::TetraType::HasQuality(); }
    inline bool IsQuality3Enabled() const { return T::TetraType::HasQuality3(); }

    //Empty flags
    int &Flags() { static int dummyflags(0); assert(0); return dummyflags; }
    int cFlags() const { return 0; }

    static bool HasFlags()      { return false; }
    static bool HasFlagsOcc()   { return false; }

    //Empty IMark
    typedef int MarkType;

    inline void InitIMark()   {  }
    inline int & IMark()      { assert(0); static int tmp=-1; return tmp;}
    inline int cIMark() const {return 0;}

    inline bool IsMarkEnabled() const { return T::TetraType::HasMark(); }
    
    static bool HasMark()     { return false; }
    static bool HasMarkOcc()  { return false; }

    //Empty Adjacency
    typedef int VTAdjType;
    typename T::TetraPointer & VTp     ( const int )       { static typename T::TetraPointer tp=0; assert(0); return tp; }
    typename T::TetraPointer const cVTp( const int ) const { static typename T::TetraPointer const tp=0; assert(0); return tp; }

    typename T::TetraPointer & TTp     ( const int )       { static typename T::TetraPointer tp=0; assert(0); return tp; }
    typename T::TetraPointer const cTTp( const int ) const { static typename T::TetraPointer const tp=0; assert(0); return tp; }

    char & VTi( const int )       { static char z=0; assert(0); return z; }
    char   VTi( const int ) const { static char z=0; assert(0); return z; }
    char  cVTi( const int ) const { static char z=0; assert(0); return z; }
    char & TTi( const int )       { static char z=0; assert(0); return z; }
    char  cTTi( const int ) const { static char z=0; assert(0); return z; }

    bool IsVTInitialized(const int j) const {return  static_cast<const typename T::TetraType *>(this)->cVTi(j)!=-1;}
    void VTClear(int j) {
        if(IsVTInitialized(j)) {
            static_cast<typename T::TetraPointer>(this)->VTp(j)=0;
            static_cast<typename T::TetraPointer>(this)->VTi(j)=-1;
        }
    }

    static bool HasVTAdjacency()    { return false; }
    static bool HasTTAdjacency()    { return false; }
    static bool HasTTAdjacencyOcc() { return false; }
    static bool HasVTAdjacencyOcc() { return false; }

    template <class RightValuteType>
    void ImportData(const RightValuteType & ) {}

    static void Name(std::vector<std::string> & name) { T::Name(name); }
};
/*-------------------------- VERTEX ----------------------------------------*/
// template <class T> class EmptyVertexRef: public T {
// public:
//  // typedef typename T::VertexType VertexType;
//  // typedef typename T::CoordType CoordType;
//   inline typename T::VertexType *       & V( const int j ) 	    {	assert(0);		static typename T::VertexType *vp=0; return vp; }
//   inline typename T::VertexType * const & V( const int j ) const {	assert(0);		static typename T::VertexType *vp=0; return vp; }
// 	inline typename T::VertexType * const  cV( const int j ) const {	assert(0);		static typename T::VertexType *vp=0; return vp;	}
// 	inline       typename T::CoordType & P( const int j ) 	    {	assert(0);		static typename T::CoordType coord(0, 0, 0); return coord;	}
// 	inline const typename T::CoordType & P( const int j ) const {	assert(0);		static typename T::CoordType coord(0, 0, 0); return coord;	}
// 	inline const typename T::CoordType &cP( const int j ) const	{	assert(0);		static typename T::CoordType coord(0, 0, 0); return coord;	}

// 	static bool HasVertexRef()   { return false; }
// };

template <class T> class VertexRef: public T {
public:
    VertexRef(){
        v[0]=0;
        v[1]=0;
        v[2]=0;
        v[3]=0;

        /******* vertex and faces indices scheme*********
 *
 * 						/2\`
 * 					 /   \  `
 * 					/     \    `
 * 				 /       \  _ 3`
 *				/     _   \    '
 * 			 / _         \  '
 *      /0___________1\'
 *
 */
        findices[0][0] = 0; findices[0][1] = 1; findices[0][2] = 2;
        findices[1][0] = 0; findices[1][1] = 3; findices[1][2] = 1;
        findices[2][0] = 0; findices[2][1] = 2; findices[2][2] = 3;
        findices[3][0] = 1; findices[3][1] = 3; findices[3][2] = 2;
    }

    typedef typename T::VertexType::CoordType CoordType;
    typedef typename T::VertexType::ScalarType ScalarType;

    inline       typename T::VertexType *       & V( const int j )  { assert(j>=0 && j<4); return v[j]; }
    inline const typename T::VertexType *        cV( const int j )  { assert(j>=0 && j<4); return v[j]; }

    inline  size_t cFtoVi (const int f, const int j) const { assert(f >= 0 && f < 4); assert(j >= 0 && j < 3); return findices[f][j]; }

    // Shortcut for tetra points
    inline       CoordType & P( const int j )       { assert(j>=0 && j<4); return v[j]->P();	}
    inline const CoordType &cP( const int j ) const { assert(j>=0 && j<4); return v[j]->P(); }

    /** Return the pointer to the ((j+1)%4)-th vertex of the tetra.
                @param j Index of the face vertex.
         */
    inline       typename T::VertexType *       &  V0( const int j )       { return V(j);}
    inline       typename T::VertexType *       &  V1( const int j )       { return V((j+1)%4);}
    inline       typename T::VertexType *       &  V2( const int j )       { return V((j+2)%4);}
    inline       typename T::VertexType *       &  V3( const int j )       { return V((j+3)%4);}
    inline const typename T::VertexType * const &  V0( const int j ) const { return V(j);}
    inline const typename T::VertexType * const &  V1( const int j ) const { return V((j+1)%4);}
    inline const typename T::VertexType * const &  V2( const int j ) const { return V((j+2)%4);}
    inline const typename T::VertexType * const &  V3( const int j ) const { return V((j+3)%4);}
    inline const typename T::VertexType * const & cV0( const int j ) const { return cV(j);}
    inline const typename T::VertexType * const & cV1( const int j ) const { return cV((j+1)%4);}
    inline const typename T::VertexType * const & cV2( const int j ) const { return cV((j+2)%4);}
    inline const typename T::VertexType * const & cV3( const int j ) const { return cV((j+3)%4);}

    /// Shortcut to get vertex values
    inline       CoordType &P0 (const int j)       { return V(j)->P(); }
    inline       CoordType &P2 (const int j)       { return V((j + 2) % 4)->P(); }
    inline       CoordType &P3 (const int j)       { return V((j + 3) % 4)->P(); }
    inline       CoordType &P1 (const int j)       { return V((j + 1) % 4)->P(); }
    inline const CoordType &P0 (const int j) const { return V(j)->P(); }
    inline const CoordType &P1 (const int j) const { return V((j + 1) % 4)->P(); }
    inline const CoordType &P2 (const int j) const { return V((j + 2) % 4)->P(); }
    inline const CoordType &P3 (const int j) const { return V((j + 3) % 4)->P(); }
    inline const CoordType &cP0(const int j) const { return cV(j)->P(); }
    inline const CoordType &cP1(const int j) const { return cV((j + 1) % 4)->P(); }
    inline const CoordType &cP2(const int j) const { return cV((j + 2) % 4)->P(); }
    inline const CoordType &cP3(const int j) const { return cV((j + 3) % 4)->P(); }

    static bool HasVertexRef()   { return true; }
    static bool HasTVAdjacency() { return true; }

    static void Name(std::vector<std::string> & name){name.push_back(std::string("VertexRef"));T::Name(name);}

    template <class RightValueType>
    void ImportData(const RightValueType & rTetra) { T::ImportData(rTetra); }

private:
    typename T::VertexType *v[4];
    size_t findices[4][3];
};


/*------------------------- FACE NORMAL -----------------------------------------*/
// template <class A, class T> class EmptyFaceNormal: public T {
// public:
// 	typedef ::vcg::Point3<A> NormalType;
// 	/// Return the vector of Flags(), senza effettuare controlli sui bit
// 	NormalType N(const int & ){ static int dummynormal(0); return dummynormal; }
//   const NormalType cN(const int & ) const { return 0; }
//   static bool HasFaceNormal()   { return false; }
//   static bool HasFaceNormalOcc()   { return false; }
// 	static void Name(std::vector<std::string> & name){T::Name(name);}

// };

// template <class A, class T> class FaceNormal:  public T {
// public:
// 	typedef A NormalType;

// 	inline NormalType N(const int & i){  assert((i>=0)&&(i < 4)); return _facenormals[i]; }
//   inline NormalType cN(const int & i) const { assert((i>=0)&&(i < 4)); return _facenormals[i]; }
//   static bool HasFaceNormals()   { return true; }
//   static bool HasFaceNormalOcc()   { return false; }

//   template <class RightValueType>
//   void ImportData(const RightValueType & rightT)
//   {
//     if(rightT.IsNormalEnabled()) N().Import(rightT.cN());
//     T::ImportData(rightT);
//   }

//   static void Name(std::vector<std::string> & name){name.push_back(std::string("FaceNormal"));T::Name(name);}

// private:
//   NormalType  _facenormals[4];
// };

//template <class T> class FaceNormal3f: public FaceNormal<vcg::Point3f, T>{
//public:static void Name(std::vector<std::string> & name){name.push_back(std::string("FaceNormal3f"));T::Name(name);} };

//template <class T> class FaceNormal3d: public FaceNormal<vcg::Point3d, T>{
//public:static void Name(std::vector<std::string> & name){name.push_back(std::string("FaceNormal3d"));T::Name(name);} };

/*------------------------- FLAGS -----------------------------------------*/
// template <class T> class EmptyBitFlags: public T {
// public:
// 	/// Return the vector of Flags(), senza effettuare controlli sui bit
//   int &Flags() { static int dummyflags(0); return dummyflags; }
//   const int Flags() const { return 0; }
//   static bool HasFlags()   { return false; }
//   static bool HasFlagsOcc()   { return false; }
// 	static void Name(std::vector<std::string> & name){T::Name(name);}

// };

template <class T> class BitFlags:  public T {
public:
    typedef int FlagType;
    BitFlags(){_flags=0;}
    
		inline int &Flags() {return _flags; }
    inline int cFlags() const {return _flags; }

    template <class RightValueType>
    void ImportData(const RightValueType & rightT){
        if(RightValueType::HasFlags())
            Flags() = rightT.cFlags();
        T::ImportData(rightT);
    }

    static bool HasFlags()   { return true; }
    static void Name(std::vector<std::string> & name){name.push_back(std::string("BitFlags"));T::Name(name);}


private:
    int  _flags;
};

/*-------------------------- QUALITY  ----------------------------------------*/
template <class A, class T> class Quality: public T {
public:
    typedef A QualityType;
    Quality():_quality(0) {}
    QualityType &Q()       { return _quality; }
    QualityType cQ() const { return _quality; }

    template <class RightValueType>
    void ImportData(const RightValueType & rightT){
        if(rightT.IsQualityEnabled())
            Q() = rightT.cQ();
        T::ImportData(rightT);
    }

    static bool HasQuality()   { return true; }
    static void Name(std::vector<std::string> & name){name.push_back(std::string("Quality"));T::Name(name);}
private:
    QualityType _quality;
};

template <class T> class Qualityf: public Quality<float, T> {
public:  static void Name(std::vector<std::string> & name){name.push_back(std::string("Qualityf"));T::Name(name);}
};
template <class T> class Qualityd: public Quality<double, T> {
public:  static void Name(std::vector<std::string> & name){name.push_back(std::string("Qualityd"));T::Name(name);}
};

/*-------------------------- Quality3  ----------------------------------*/
template <class A, class T> class Quality3: public T {
public:
    typedef vcg::Point3<A> Quality3Type;
    Quality3Type &Q3()       { return _quality; }
    Quality3Type cQ3() const { return _quality; }
    template <class RightValueType>
    void ImportData(const RightValueType & rightT){
        if(rightT.IsQuality3Enabled()) Q3() = rightT.cQ3();
        T::ImportData(rightT);
    }

    static bool HasQuality3()   { return true; }
    static void Name(std::vector<std::string> & name){name.push_back(std::string("Quality3"));T::Name(name);}
private:
    Quality3Type _quality;
};

template <class T> class Quality3f: public Quality3<float, T> {
public:  static void Name(std::vector<std::string> & name){name.push_back(std::string("Quality3f"));T::Name(name);}
};
template <class T> class Quality3d: public Quality3<double, T> {
public:  static void Name(std::vector<std::string> & name){name.push_back(std::string("Quality3d"));T::Name(name);}
};

/*-------------------------- COLOR  ----------------------------------------*/
template <class A, class T>
class Color: public T {
public:
    typedef A ColorType;
    Color():_color(vcg::Color4b::White) {}
    ColorType &C()       { return _color; }
    ColorType cC() const { return _color; }
    template <class RightValueType>
    void ImportData(const RightValueType & rightT){
        if(rightT.IsColorEnabled()) C() = rightT.cC();
        T::ImportData(rightT);
    }

    static bool HasColor()   { return true; }
    static void Name(std::vector<std::string> & name){name.push_back(std::string("Color"));T::Name(name);}

private:
    ColorType _color;
};

template <class T> class Color4b : public Color<vcg::Color4b, T> {
public: static void Name(std::vector<std::string> & name){name.push_back(std::string("Color4b"));T::Name(name); }
};
/*-------------------------- INCREMENTAL MARK  ----------------------------------------*/

// template <class T> class EmptyMark: public T {
// public:
// 	typedef int MarkType;
//   static bool HasMark()   { return false; }
//   static bool HasMarkOcc()   { return false; }
//   inline void InitIMark()    {  }
//   inline int & IMark()       { assert(0); static int tmp=-1; return tmp;}
//   inline const int IMark() const {return 0;}
//   static void Name(std::vector<std::string> & name){T::Name(name);}

// };

template <class T> class Mark: public T {
public:
    static bool HasMark()      { return true; }
    static bool HasMarkOcc()   { return false; }
    inline void InitIMark()    { _imark = 0; }
    inline int & IMark()       { return _imark;}
    inline int  cIMark() const {return _imark;}

    template <class RightValueType>
    void ImportData(const RightValueType & rightT){
        if(rightT.IsMarkEnabled()) IMark() = rightT.cIMark();
        T::ImportData(rightT);
    }

    static void Name(std::vector<std::string> & name){name.push_back(std::string("Mark"));T::Name(name);}

private:
    int _imark;
};


/*----------------------------- VTADJ ------------------------------*/

// template <class T> class EmptyAdj: public T {
// public:
// 	typedef int VFAdjType;
// 	typename T::TetraPointer & VTp( const int ) { static typename T::TetraPointer tp=0; return tp; }
// 	typename T::TetraPointer const cVTp( const int ) const { static typename T::TetraPointer const tp=0; return tp; }
// 	typename T::TetraPointer & TTp( const int ) { static typename T::TetraPointer tp=0; return tp; }
// 	typename T::TetraPointer const cTTp( const int ) const { static typename T::TetraPointer const tp=0; return tp; }
// 	char & VTi( const int j ) { static char z=0; return z; }
// 	char & TTi( const int j ) { static char z=0; return z; }
// 	static bool HasVTAdjacency() { return false; }
// 	static bool HasTTAdjacency() { return false; }
// 	static bool HasTTAdjacencyOcc() { return false; }
// 	static bool HasVTAdjacencyOcc() { return false; }
// 	static void Name( std::vector< std::string > & name ){ T::Name(name); }
// };

template <class T> class VTAdj: public T {
public:
    VTAdj() {
        _vtp[0]=0;
        _vtp[1]=0;
        _vtp[2]=0;
        _vtp[3]=0;
        _vti[0]=-1;
        _vti[1]=-1;
        _vti[2]=-1;
        _vti[3]=-1;
    }

    typename T::TetraPointer & VTp( const int j ) { assert( j >= 0 && j < 4 ); return _vtp[j]; }
    typename T::TetraPointer const VTp( const int j ) const { assert( j >= 0 && j < 4 ); return _vtp[j]; }
    typename T::TetraPointer const cVTp( const int j ) const { assert( j >= 0 && j < 4 ); return _vtp[j]; }

    char & VTi( const int j ) { return _vti[j]; }
    const char & cVTi( const int j ) const { return _vti[j]; }

    static bool HasVTAdjacency() { return true; }
    static bool HasVTAdjacencyOcc() { return false; }

    static void Name( std::vector< std::string > & name ) { name.push_back( std::string("VTAdj") ); T::Name(name); }

    template <class RightValueType>
    void ImportData(const RightValueType & rightT){T::ImportData(rightT);}
private:
    typename T::TetraPointer _vtp[4];
    char _vti[4];
};

/*----------------------------- TTADJ ------------------------------*/

template <class T> class TTAdj: public T {
public:
    TTAdj(){
        _ttp[0]=0;
        _ttp[1]=0;
        _ttp[2]=0;
        _ttp[3]=0;
    }
    typename T::TetraPointer       &TTp(const int j)        { assert(j>=0 && j<4);  return _ttp[j]; }
    typename T::TetraPointer const  TTp(const int j) const  { assert(j>=0 && j<4);  return _ttp[j]; }
    typename T::TetraPointer const cTTp(const int j) const  { assert(j>=0 && j<4);  return _ttp[j]; }
    char        &TTi(const int j)       { return _tti[j]; }
    const char &cTTi(const int j) const { return _tti[j]; }

    typename T::TetraPointer        &TTp1( const int j )       { return TTp((j+1)%4);}
    typename T::TetraPointer        &TTp2( const int j )       { return TTp((j+2)%4);}
    typename T::TetraPointer        &TTp3( const int j )       { return TTp((j+3)%4);}
    typename T::TetraPointer  const  TTp1( const int j ) const { return TTp((j+1)%4);}
    typename T::TetraPointer  const  TTp2( const int j ) const { return TTp((j+2)%4);}
    typename T::TetraPointer  const  TTp3( const int j ) const { return TTp((j+3)%4);}

    bool IsBorderF(const int & i)  const { assert( (i>=0) && (i < 4)); { return TTp(i) == this;}}

    static bool HasTTAdjacency()      {   return true; }
    static bool HasTTAdjacencyOcc()   {   return false; }
    static void Name(std::vector<std::string> & name){name.push_back(std::string("TTAdj"));T::Name(name);}

    template <class RightValueType>
    void ImportData(const RightValueType & rightT){T::ImportData(rightT);}

private:
    typename T::TetraPointer _ttp[4] ;
    char _tti[4] ;
};

} // end namespace vert
}// end namespace vcg
#endif
