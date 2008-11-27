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

#ifndef __VCG_RandomGenerator
#define __VCG_RandomGenerator

#include <functional>
namespace vcg {
  namespace math {

/**
 * RandomGenerator includes two Random Number Generation algortihms.
 *
 * The first one is derived from a STL extension of sgi: 
 * it is based on the Subtractive Ring method. 
 * Note: this code assumes that int is 32 bits.
 *
 * The second one is an improved Marsenne-Twister algorithm (MT19937)
 * Coded by Takuji Nishimura and Makoto Matsumoto (see copyright note below)
 * and successively modified to be a C++ class by Daniel Dunbar.
 *
 * References for Subtractive Ring:
 *
 * D. E. Knuth, The Art of Computer Programming. Volume 2: Seminumerical Algorithms, 2nd Edition. Addison-Wesley, 1981.
 *  (section 3.6 of Knuth for an implementation of the subtractive method in FORTRAN)
 *  (section 3.2.2 of Knuth analyzes this class of algorithms)
 *
 * References for improved Marsenne-Twister:
 *
 *   http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html
 *
 */
class RandomGenerator
{

// definitions (used by the improved Marsenne-Twister algorithm)
private:

	static const long N = 624;
	static const long M = 397;
	static const unsigned long MATRIX_A = 0x9908b0dfUL;   // constant vector a
	static const unsigned long UPPER_MASK = 0x80000000UL; // most significant w-r bits
	static const unsigned long LOWER_MASK = 0x7fffffffUL; // least significant r bits

// private data member
private:

	// Subtractive Ring RNG status variables
	unsigned int _M_table[55];
	size_t _M_index1;
	size_t _M_index2;

	// Improved Marsenne-Twister RNG status variables
	unsigned long mt[N]; // the array for the state vector
	int mti;

// public methods
public:

	/// Initialize Subtractive Ring RNG with a given seed.
	void initializeSubtractiveRing(unsigned int __seed = 161803398u)
	{
		unsigned int __k = 1;
		_M_table[54] = __seed;
		size_t __i;
		for (__i = 0; __i < 54; __i++) 
		{
			size_t __ii = (21 * (__i + 1) % 55) - 1;
			_M_table[__ii] = __k;
			__k = __seed - __k;
			__seed = _M_table[__ii];
		}
		for (int __loop = 0; __loop < 4; __loop++) 
		{
			for (__i = 0; __i < 55; __i++)
				_M_table[__i] = _M_table[__i] - _M_table[(1 + __i + 30) % 55];
		}
		_M_index1 = 0;
		_M_index2 = 31;
	}

	/// Return a random number in the given range (__limit) using the Subtractive Ring method.
	unsigned int generateWithSubtractiveRing(unsigned int __limit)
	{
		_M_index1 = (_M_index1 + 1) % 55;
		_M_index2 = (_M_index2 + 1) % 55;
		_M_table[_M_index1] = _M_table[_M_index1] - _M_table[_M_index2];
		return _M_table[_M_index1] % __limit;
	}

	/// Initialize Improved Marsenne Twister RNG with the given seed.
	void initializeImprovedMarsenneTwister(unsigned long seed = 5489UL)
	{
		mt[0]= seed & 0xffffffffUL;
		for (mti=1; mti<N; mti++) 
		{
			mt[mti] = (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti); 
			/* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
			/* In the previous versions, MSBs of the seed affect   */
			/* only MSBs of the array mt[].                        */
			/* 2002/01/09 modified by Makoto Matsumoto             */
			mt[mti] &= 0xffffffffUL;
			/* for >32 bit machines */
		}
	}

	/** 
	 * Initialize by an array with array-length.
	 *
	 * init_key is the array for initializing keys
	 * key_length is its length
	 */
	void initializeImprovedMarsenneTwister(unsigned long init_key[], int key_length)
	{
		int i, j, k;
		initializeImprovedMarsenneTwister(19650218UL);
		i=1; j=0;
		k = (N>key_length ? N : key_length);
		for (; k; k--) 
		{
			mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1664525UL)) + init_key[j] + j; /* non linear */
			mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
			i++; j++;
			
			if (i>=N) 
			{ 
				mt[0] = mt[N-1]; 
				i=1; 
			}

			if (j>=key_length) j=0;
		}

		for (k=N-1; k; k--) 
		{
			mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1566083941UL)) - i; /* non linear */
			mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
			i++;
			if (i>=N) 
			{ 
				mt[0] = mt[N-1]; 
				i=1; 
			}
		}

		mt[0] = 0x80000000UL; /* MSB is 1; assuring non-zero initial array */ 
	}

	/// Return a random number in the [0,0xffffffff] interval using the improved Marsenne Twister algorithm.
	unsigned long generateWithImprovedMarsenneTwister()
	{
		unsigned long y;
		static unsigned long mag01[2]={0x0UL, MATRIX_A};
		/* mag01[x] = x * MATRIX_A  for x=0,1 */

		if (mti >= N) // generate N words at one time
		{
			int kk;

			for (kk=0;kk<N-M;kk++) 
			{
				y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
				mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1UL];
			}

			for (;kk<N-1;kk++) 
			{
				y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
				mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
			}

			y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
			mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];

			mti = 0;
		}

		y = mt[mti++];

		/* Tempering */
		y ^= (y >> 11);
		y ^= (y << 7) & 0x9d2c5680UL;
		y ^= (y << 15) & 0xefc60000UL;
		y ^= (y >> 18);

		return y;
	}

	/// Generates a random number in the [0,1] real interval using the improved Marsenne-Twister.
	double generateDoubleLRwithImprovedMT()
	{
		return generateWithImprovedMarsenneTwister()*(1.0/4294967295.0); 
	}

	/// Generates a random number in the [0,1) real interval using the improved Marsenne-Twister.
	double generateDoubleLwithImprovedMT()
	{
		return generateWithImprovedMarsenneTwister()*(1.0/4294967296.0); 
	}

	/// Generates a random number in the (0,1) real interval using the improved Marsenne-Twister.
	double generateDoubleWithImprovedMT()
	{
		return (((double)generateWithImprovedMarsenneTwister()) + 0.5)*(1.0/4294967296.0);
	}


// construction
public:

	// ctor
	RandomGenerator(){}

};

} // end namespace math
} // end namespace vcg



/*
   Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
   All rights reserved.                          

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:

     1. Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.

     2. Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.

     3. The names of its contributors may not be used to endorse or promote 
        products derived from this software without specific prior written 
        permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/


#endif /* __VCG_RandomGenerator */
