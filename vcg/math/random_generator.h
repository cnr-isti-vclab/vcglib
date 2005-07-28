/*#***************************************************************************
 * RandomGenerator.h                                                o o      *
 *                                                                o     o    *
 * Visual Computing Group                                         _  O  _    *
 * IEI Institute, CNUCE Institute, CNR Pisa                        \/)\/     *
 *                                                                /\/|       *
 * Copyright(C) 1999 by Paolo Cignoni, Paolo Pingi, Claudio Rocchini |       *
 * All rights reserved.                                              \       *
 *                                                                           *
 * Permission  to use, copy, modify, distribute  and sell this  software and *
 * its documentation for any purpose is hereby granted without fee, provided *
 * that  the above copyright notice appear  in all copies and that both that *
 * copyright   notice  and  this  permission  notice  appear  in  supporting *
 * documentation. the author makes  no representations about the suitability *
 * of this software for any purpose. It is provided  "as is" without express *
 * or implied warranty.                                                      *
 *                                                                           *
 *****************************************************************************/
/****************************************************************************
  History
$Log: not supported by cvs2svn $         
 *****************************************************************************/

// RandomGenerator is derived from a STL extension of sgi: 
// it is based on the Subtractive Ring method. 
// See section 3.6 of Knuth for an implementation of the subtractive method in FORTRAN. 
// Section 3.2.2 of Knuth analyzes this class of algorithms. 
// (D. E. Knuth, The Art of Computer Programming. Volume 2: Seminumerical Algorithms, second edition. Addison-Wesley, 1981.) .
// Note: this code assumes that int is 32 bits.

#ifndef __VCG_RandomGenerator
#define __VCG_RandomGenerator
namespace vcg {
  namespace math {

    class RandomGenerator : public std::unary_function<unsigned int, unsigned int> {
private:
  unsigned int _M_table[55];
  size_t _M_index1;
  size_t _M_index2;
public:
  unsigned int operator()(unsigned int __limit) {
    _M_index1 = (_M_index1 + 1) % 55;
    _M_index2 = (_M_index2 + 1) % 55;
    _M_table[_M_index1] = _M_table[_M_index1] - _M_table[_M_index2];
    return _M_table[_M_index1] % __limit;
  }

  void _M_initialize(unsigned int __seed)
  {
    unsigned int __k = 1;
    _M_table[54] = __seed;
    size_t __i;
    for (__i = 0; __i < 54; __i++) {
        size_t __ii = (21 * (__i + 1) % 55) - 1;
        _M_table[__ii] = __k;
        __k = __seed - __k;
        __seed = _M_table[__ii];
    }
    for (int __loop = 0; __loop < 4; __loop++) {
        for (__i = 0; __i < 55; __i++)
            _M_table[__i] = _M_table[__i] - _M_table[(1 + __i + 30) % 55];
    }
    _M_index1 = 0;
    _M_index2 = 31;
  }
	
	RandomGenerator(unsigned int __seed) { _M_initialize(__seed); }
  RandomGenerator() { _M_initialize(161803398u); }

};
  } // end namespace math
} // end namespace vcg


#endif