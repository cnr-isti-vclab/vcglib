#ifndef RECURSE_H
#define RECURSE_H

#include <vcg/space/point3.h>


class Empty{};

/*------------------------------------------------------------------*/ 
template <class T> class EmptyCoord: public T {
public:
  typedef vcg::Point3f CoordType;
  CoordType &P() { static CoordType coord(0, 0, 0); return coord; }
};
template <class A, class T> class Coord: public T {
public:
  typedef A CoordType;
  CoordType &P() { return coord; }
private:
  CoordType coord;    
};
template <class T> class Coord3f: public Coord<vcg::Point3f, T> {};
template <class T> class Coord3d: public Coord<vcg::Point3d, T> {};

/*------------------------------------------------------------------*/ 

template <class T> class EmptyNormal: public T {
public:
  typedef vcg::Point3s NormalType;
  NormalType &N() { static NormalType normal(0, 0, 0); return normal; }
};
template <class A, class T> class NormalX: public T {
public:
  typedef A NormalType;
  NormalType &N() { return norm; }
private:
  NormalType norm;    
};
template <class T> class Normal3s: public NormalX<vcg::Point3s, T> {};
template <class T> class Normal3f: public NormalX<vcg::Point3f, T> {};
template <class T> class Normal3d: public NormalX<vcg::Point3d, T> {};

/*------------------------------------------------------------------*/ 

template <class T> class EmptyColor: public T {
public:
  typedef vcg::Point3s ColorType;
  ColorType &C() { static ColorType color(0, 0, 0); return color; }
};
template <class A, class T> class Color: public T {
public:
  typedef A ColorType;
  ColorType &C() { return color; }
private:
  ColorType color;    
};
template <class T> class Color3c: public Color<vcg::Point3<char>, T> {};


/*------------------------------------------------------------------*/ 


template <class T> class EmptyVFAdj: public T {
public:
  //typedef vcg::Point3s ColorType;
  typename T::FacePointer &Fp() { static typename T::FacePointer fp=0; return fp; }
  int &Zp(){static int z=0; return z;};
};

template <class T> class VFAdj: public T {
public:
  //typedef A ColorType;
  typename T::FacePointer &Fp() {return fp; }
  typename T::FacePointer &Zp() {return Zp; }
private:
  typename T::FacePointer fp ;    
  int zp ;    
};


/*------------------------------------------------------------------*/ 
template <class BVT, typename BFT>
class FVEmpty{
public:
typedef BVT *VertPointer;
typedef BFT *FacePointer;
};

template <class BVT, typename BFT>
class Vertex: public EmptyVFAdj<EmptyColor<EmptyNormal<EmptyCoord<FVEmpty <BVT,BFT> > > > >{
public:
};




template <class BVT, typename BFT,
          template <typename> class A> class Vertex1: public A<Vertex<BVT,BFT> > {};

template <class BVT, typename BFT, 
          template <typename> class A, 
          template <typename> class B> class Vertex2: public B<Vertex1<BVT,BFT, A> > {};

template <class BVT, typename BFT,
          template <typename> class A, 
          template <typename> class B, 
          template <typename> class C > class Vertex3: public C<Vertex2<BVT,BFT, A, B> > {};

template <class BVT, typename BFT,
          template <typename> class A, 
          template <typename> class B, 
          template <typename> class C, 
          template <typename> class D > class Vertex4: public D<Vertex3<BVT,BFT, A, B, C> > {};


template < typename T>
class Default : public T {};

template <typename BVT, typename BFT,
          template <typename> class A, 
          template <typename> class B = Default, 
          template <typename> class C = Default,
          template <typename> class D = Default> class VertexK: public Default<Vertex4<BVT,BFT, A, B, C, D> > {};


template < typename T=int>
class DefaultQ : public T {};

template <typename BVT, typename BFT,
          template <typename> class A, 
          template <typename> class B = DefaultQ, 
          template <typename> class C = DefaultQ,
          template <typename> class D = DefaultQ > class VertexQ: public Vertex4<BVT,BFT, A, B, C, D>  {};

template <typename BVT, typename BFT,
          template <typename> class A
          > class VertexQ<BVT, BFT, A, DefaultQ, DefaultQ, DefaultQ>: public Vertex1<BVT,BFT, A>  {};

template <typename BVT, typename BFT,
          template <typename> class A, 
          template <typename> class B
          > class VertexQ<BVT, BFT, A, B, DefaultQ, DefaultQ>: public Vertex2<BVT,BFT, A,B>  {};

template <typename BVT, typename BFT,
          template <typename> class A, 
          template <typename> class B, 
          template <typename> class C
          > class VertexQ<BVT, BFT, A, B, C, DefaultQ>: public Vertex3<BVT,BFT, A,B,C>  {};

//template <typename BFT,
//          template <typename> class A, 
//          template <typename> class B = Default, 
//          template <typename> class C = Default,
//          template <typename> class D = Default > class VertexZ: public Default<Vertex3< VertexZ, BFT, A, B, C> > {};

/* un vertice si definisce cosi':

Vertex2<Coord3f, Normal3f> myvertex;

*/

#endif
