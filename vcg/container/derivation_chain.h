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

#ifndef __VCG_DERIVATION_CHAIN
#define __VCG_DERIVATION_CHAIN

namespace vcg{
/*------------------------------------------------------------------*/ 

// Metaprogramming Core
template <template <typename,typename,typename,typename>
					class Base,class BVT, class BET, class BFT,class BTT,
          template <typename> class A> 
          class Arity1: public A<Base<BVT,BET,BFT,BTT> > {
          };

template <template <typename,typename,typename,typename>
					class Base,class BVT, class BET, typename BFT, class BTT,
          template <typename> class A, template <typename> class B> 
          class Arity2: public B<Arity1<Base,BVT,BET,BFT,BTT, A> > {};

template <template <typename,typename,typename,typename>
					class Base,class BVT, class BET, typename BFT,class BTT,
          template <typename> class A, template <typename> class B, 
          template <typename> class C > 
          class Arity3: public C<Arity2<Base,BVT,BET,BFT,BTT, A, B> > {};

template <template <typename,typename,typename,typename>
					class Base,class BVT, class BET, typename BFT,class BTT,
          template <typename> class A, template <typename> class B, 
          template <typename> class C, template <typename> class D> 
          class Arity4: public D<Arity3<Base,BVT,BET,BFT,BTT, A, B, C> > {};

template <template <typename,typename,typename,typename>
					class Base,class BVT, class BET, typename BFT,class BTT,
          template <typename> class A, template <typename> class B, 
          template <typename> class C, template <typename> class D,
          template <typename> class E > 
          class Arity5: public E<Arity4<Base,BVT,BET,BFT,BTT, A, B, C, D> > {};

template <template <typename,typename,typename,typename>
					class Base,class BVT, class BET, typename BFT,class BTT,
          template <typename> class A, template <typename> class B, 
          template <typename> class C, template <typename> class D,
          template <typename> class E, template <typename> class F > 
          class Arity6: public F<Arity5<Base,BVT,BET,BFT,BTT, A, B, C, D, E> > {};

template <template <typename,typename,typename,typename>
					class Base,class BVT, class BET, typename BFT,class BTT,
          template <typename> class A, template <typename> class B, 
          template <typename> class C, template <typename> class D,
          template <typename> class E, template <typename> class F,
					template <typename> class G> 
          class Arity7: public G<Arity6<Base,BVT,BET,BFT,BTT, A, B, C, D, E, F> > {};

template <template <typename,typename,typename,typename>
					class Base,class BVT, class BET, typename BFT,class BTT,
          template <typename> class A, template <typename> class B, 
          template <typename> class C, template <typename> class D,
          template <typename> class E, template <typename> class F,
					template <typename> class G, template <typename> class H> 
          class Arity8: public H<Arity7<Base,BVT,BET,BFT,BTT, A, B, C, D, E, F, G > > {};

template <template <typename,typename,typename,typename>
					class Base,class BVT, class BET, typename BFT,class BTT,
          template <typename> class A, template <typename> class B, 
          template <typename> class C, template <typename> class D,
          template <typename> class E, template <typename> class F,
					template <typename> class G, template <typename> class H,
					template <typename> class I>
          class Arity9: public I<Arity8<Base,BVT,BET,BFT,BTT, A, B, C, D, E, F, G, H > > {};
					
template <template <typename,typename,typename,typename>
					class Base,class BVT, class BET, typename BFT,class BTT,
          template <typename> class A, template <typename> class B, 
          template <typename> class C, template <typename> class D,
          template <typename> class E, template <typename> class F,
					template <typename> class G, template <typename> class H,
					template <typename> class I, template <typename> class J>
          class Arity10: public J<Arity9<Base,BVT,BET,BFT,BTT, A, B, C, D, E, F, G, H, I > > {};

template <template <typename,typename,typename,typename>
          class Base,class BVT, class BET, typename BFT,class BTT,
          template <typename> class A, template <typename> class B, 
          template <typename> class C, template <typename> class D,
          template <typename> class E, template <typename> class F,
          template <typename> class G, template <typename> class H,
          template <typename> class I, template <typename> class J,
          template <typename> class K>
          class Arity11: public K<Arity10<Base,BVT,BET,BFT,BTT, A, B, C, D, E, F, G, H, I, J> > {};


template <template <typename,typename,typename,typename>
          class Base,class BVT, class BET, typename BFT,class BTT,
          template <typename> class A, template <typename> class B, 
          template <typename> class C, template <typename> class D,
          template <typename> class E, template <typename> class F,
          template <typename> class G, template <typename> class H,
          template <typename> class I, template <typename> class J,
          template <typename> class K, template <typename> class L>
          class Arity12: public L<Arity11<Base,BVT,BET,BFT,BTT, A, B, C, D, E, F, G, H, I, J, K> > {};

template < typename T=int>
class DefaultDeriver : public T {};
class DumClass {};

}// end namespace vcg
#endif
