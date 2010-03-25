#ifndef VCG_USED_TYPES_H
#define VCG_USED_TYPES_H

#include <string>
#include <vcg/space/point3.h>
#include <vcg/container/derivation_chain.h>

namespace vcg{

struct DummyTypes{
		typedef char VertexType; 		// simplex types
		typedef char EdgeType;
		typedef char FaceType;
		typedef char TetraType;
		typedef char HEdgeType; 		// connector types

		typedef vcg::Point3<bool> CoordType; 		 
		typedef char ScalarType;						 

		typedef VertexType * VertexPointer;
		typedef EdgeType *	EdgePointer		;
		typedef FaceType * FacePointer		;
		typedef TetraType * TetraPointer	;
		typedef HEdgeType * HEdgePointer	;

    static void Name(std::vector<std::string> & /*name*/){}
		template < class LeftV>
    void ImportLocal(const LeftV  & /*left*/ ) {}
};

struct AllTypes{
		struct AVertexType {};
		struct AEdgeType {};
		struct AFaceType {};
		struct AHEdgeType {};
};

template <template <typename> class A = DefaultDeriver, template <typename> class B = DefaultDeriver,
					template <typename> class C = DefaultDeriver, template <typename> class D = DefaultDeriver,
					template <typename> class E = DefaultDeriver, template <typename> class F = DefaultDeriver,
					template <typename> class G = DefaultDeriver, template <typename> class H = DefaultDeriver,
					template <typename> class I = DefaultDeriver, template <typename> class J = DefaultDeriver,
					template <typename> class K = DefaultDeriver, template <typename> class L = DefaultDeriver>
					class UsedTypes: public Arity12<DummyTypes, A, B, C, D, E, F, G, H, I, J, K, L>  {
};


template <class A>
	struct Use{
		template <class T> struct AsVertexType: public T{typedef A VertexType;	typedef VertexType * VertexPointer	;};
		template <class T> struct AsEdgeType: public T{typedef A EdgeType;			typedef EdgeType *	EdgePointer			;};
		template <class T> struct AsFaceType: public T{typedef A FaceType;			typedef FaceType * FacePointer			;};
		template <class T> struct AsTetraType: public T{typedef A TetraType;		typedef TetraType * TetraPointer		;};
		template <class T> struct AsHEdgeType: public T{typedef A HEdgeType;		typedef HEdgeType * HEdgePointer		;};
};

};

#endif // USED_TYPES_H
