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

#include <vcg/math/base.h>

namespace vcg {
namespace math {

/**
 * Common interface for random generation (with uniform distribution).
 *
 * Two RNGs are available: Subtractive Ring and an improved Marsenne-Twister.
 */
class RandomGenerator
{

// construction
public:

	RandomGenerator(){}

	virtual ~RandomGenerator()
	{}

// public methods
public:

	/// (Re-)initialize with a given seed.
	virtual void initialize(unsigned int seed)=0;

	/// Return a random number in the given range (note that not all the RNG can handle a given limit).
	virtual unsigned int generate(unsigned int limit)=0;

	/// Return a random number in the [0,1) real interval.
	virtual double generate01()=0;

	/// Returns a random number in the [0,1] real interval.
	virtual double generate01closed()=0;

	/// Generates a random number in the (0,1) real interval.
	virtual double generate01open()=0;
};

/**
 * Uniform RNG derived from a STL extension of sgi.
 *
 * It is based on the Subtractive Ring method. 
 * This implementation assumes that int is 32 bits.
 *
 * References
 *
 * D. E. Knuth, The Art of Computer Programming. Volume 2: Seminumerical Algorithms, 2nd Edition. Addison-Wesley, 1981.
 *  (section 3.6 of Knuth for an implementation of the subtractive method in FORTRAN)
 *  (section 3.2.2 of Knuth analyzes this class of algorithms)
 */
class SubtractiveRingRNG : public RandomGenerator
{

// private data member
private:

	// Subtractive Ring RNG status variables
	unsigned int _M_table[55];
	size_t _M_index1;
	size_t _M_index2;

// construction
public:

	// ctor
	SubtractiveRingRNG(int default_seed=161803398u)
	{
		initialize(default_seed);
	}

	virtual ~SubtractiveRingRNG()
	{}

// public methods
public:

	/// (Re-)initialize with a given seed.
	void initialize(unsigned int seed)
	{
		unsigned int __k = 1;
		_M_table[54] = seed;
		size_t __i;
		for (__i = 0; __i < 54; __i++) 
		{
			size_t __ii = (21 * (__i + 1) % 55) - 1;
			_M_table[__ii] = __k;
			__k = seed - __k;
			seed = _M_table[__ii];
		}
		for (int __loop = 0; __loop < 4; __loop++) 
		{
			for (__i = 0; __i < 55; __i++)
				_M_table[__i] = _M_table[__i] - _M_table[(1 + __i + 30) % 55];
		}
		_M_index1 = 0;
		_M_index2 = 31;
	}

	/// Return a random number in the given range (limit) using the Subtractive Ring method.
	unsigned int generate(unsigned int limit= 0xffffffffu)
	{
		_M_index1 = (_M_index1 + 1) % 55;
		_M_index2 = (_M_index2 + 1) % 55;
		_M_table[_M_index1] = _M_table[_M_index1] - _M_table[_M_index2];
		return _M_table[_M_index1] % limit;
	}

	/// Return a random number in the [0,1) real interval using the Subtractive Ring method.
	double generate01()
	{
		const unsigned int lmt = 0xffffffffu;
		unsigned int number = generate(lmt);
		return static_cast<double>(number) / static_cast<double>(lmt);
	}

	/// Returns a random number in the [0,1] real interval using the Subtractive Ring method.
	double generate01closed()
	{
		const unsigned int lmt = 0xffffffffu;
		unsigned int number = generate(lmt);
		return static_cast<double>(number) / static_cast<double>(0xfffffffEu);
	}

	/// Generates a random number in the (0,1) real interval using the Subtractive Ring method.
	double generate01open()
	{
		const unsigned int lmt = 0xffffffffu;
		unsigned int number = generate(lmt);
		return (static_cast<double>(number) + 0.5) * (1.0/static_cast<double>(lmt));
	}

};

/**
 * The second one is an improved Marsenne-Twister algorithm (MT19937)
 * Coded by Takuji Nishimura and Makoto Matsumoto (see copyright note below)
 * and successively modified to be a C++ class by Daniel Dunbar.
 *
 *
 * References for improved Marsenne-Twister:
 *
 *   http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html
 *
 */
class MarsenneTwisterRNG : public RandomGenerator
{

// definitions
private:

	static const int N = 624;
	static const int M = 397;
	static const unsigned int MATRIX_A = 0x9908b0dfu;   // constant vector a
	static const unsigned int UPPER_MASK = 0x80000000u; // most significant w-r bits
	static const unsigned int LOWER_MASK = 0x7fffffffu; // least significant r bits

// private data member
private:

	// Improved Marsenne-Twister RNG status variables
	unsigned int mt[N]; // the array for the state vector
	int mti;

// construction
public:

	// ctor
	MarsenneTwisterRNG()
	{
		initialize(5489u);
	}

	virtual ~MarsenneTwisterRNG()
	{}


// public methods
public:

	/// (Re-)initialize with the given seed.
	void initialize(unsigned int seed)
	{
		mt[0]= seed & 0xffffffffu;
		for (mti=1; mti<N; mti++) 
		{
			mt[mti] = (1812433253u * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti); 
			/* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
			/* In the previous versions, MSBs of the seed affect   */
			/* only MSBs of the array mt[].                        */
			/* 2002/01/09 modified by Makoto Matsumoto             */
			mt[mti] &= 0xffffffffu;
			/* for >32 bit machines */
		}
	}

	/** 
	 * Initialize by an array with array-length.
	 *
	 * init_key is the array for initializing keys
	 * key_length is its length
	 */
	void initializeByArray(unsigned int init_key[], int key_length)
	{
		int i, j, k;
		initialize(19650218u);
		i=1; j=0;
		k = (N>key_length ? N : key_length);
		for (; k; k--) 
		{
			mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1664525u)) + init_key[j] + j; /* non linear */
			mt[i] &= 0xffffffffu; /* for WORDSIZE > 32 machines */
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
			mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1566083941u)) - i; /* non linear */
			mt[i] &= 0xffffffffu; /* for WORDSIZE > 32 machines */
			i++;
			if (i>=N) 
			{ 
				mt[0] = mt[N-1]; 
				i=1; 
			}
		}

		mt[0] = 0x80000000u; /* MSB is 1; assuring non-zero initial array */ 
	}

	/** 
	 * Return a random number in the [0,0xffffffff] interval using the improved Marsenne Twister algorithm.
	 * 
	 * NOTE: Limit is not considered, the interval is fixed.
	 */
	unsigned int generate(unsigned int /*limit*/)
	{
		unsigned int y;
		static unsigned int mag01[2]={0x0u, MATRIX_A};
		/* mag01[x] = x * MATRIX_A  for x=0,1 */

		if (mti >= N) // generate N words at one time
		{
			int kk;

			for (kk=0;kk<N-M;kk++) 
			{
				y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
				mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1u];
			}

			for (;kk<N-1;kk++) 
			{
				y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
				mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1u];
			}

			y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
			mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1u];

			mti = 0;
		}

		y = mt[mti++];

		/* Tempering */
		y ^= (y >> 11);
		y ^= (y << 7) & 0x9d2c5680u;
		y ^= (y << 15) & 0xefc60000u;
		y ^= (y >> 18);

		return y;
	}

	/// Returns a random number in the [0,1] real interval using the improved Marsenne-Twister.
	double generate01closed()
	{
		return generate(0)*(1.0/4294967295.0); 
	}

	/// Returns a random number in the [0,1) real interval using the improved Marsenne-Twister.
	double generate01()
	{
		return generate(0)*(1.0/4294967296.0); 
	}

	/// Generates a random number in the (0,1) real interval using the improved Marsenne-Twister.
	double generate01open()
	{
		return (((double)generate(0)) + 0.5)*(1.0/4294967296.0);
	}

};


/* boxmuller
 * Implements the Polar form of the Box-Muller Transformation
 *  (c) Copyright 1994, Everett F. Carter Jr.
 *  Permission is granted by the author to use this software for any
 *  application provided this copyright notice is preserved.
 */
inline double box_muller(RandomGenerator &generator, double m, double s) /* normal random variate generator */
{ /* mean m, standard deviation s */
  double x1, x2, w, y1;
  static double y2;
  static int use_last = 0;
  static RandomGenerator *last_generator = 0;
  if(last_generator != &generator)
    use_last = 0;
  last_generator = &generator;
  if (use_last){ /* use value from previous call */
    y1 = y2;
    use_last = 0;
  } else {
    do {
	  x1 = 2.0 * generator.generate01closed() - 1.0;
	  x2 = 2.0 * generator.generate01closed() - 1.0;
	  w = x1 * x1 + x2 * x2;
	} while ( w >= 1.0 );
    w = sqrt( (-2.0 * log( w ) ) / w );
    y1 = x1 * w;
    y2 = x2 * w;
   use_last = 1;
  }
  return( m + y1 * s );
}
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
