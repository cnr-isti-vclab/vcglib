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

#ifndef __VCGLIB_CALLBACK
#define __VCGLIB_CALLBACK

#include <stdarg.h>  // For va_start, etc.
#include <memory>    // For std::unique_ptr
#include <iostream>
#include <string.h>

namespace vcg {
/*! \brief This function allow lenghty algorithms to report progress and status.
 *  
 * This function is usually called from inside a long processing algorithm 
 * with a reasonable frequency (e.g. not more that 10 times per sec) 
 * reporting the what the algorithm is doing and at what point of the processing (as a int percentage) 
 * we currently are. 
 * 
 * Users of the library usually rewrote this function to handle this reporting in the correct way 
 * (e.g. in a log window and in a progress bar).
 * 
 * The algorithms in the library should also check the return value 
 * to see if the user asked for an interruption, but they can ignore this request
 * (hint: most of algs in the lib ignore any request).
 */

typedef bool CallBackPos(const int pos, const char * str );

inline bool DummyCallBackPos(const int ,  const char * ) {return true;}
inline bool COutCallBackPos(const int ,  const char *str ) {
  std::cout << str; return true;
}
inline bool CErrCallBackPos(const int ,  const char *str ) {
  std::cerr << str; return true;
}

/*! \brief Helper function that formats a string a la printf and returns a simple char*
 * 
 * Quite useful for formatting the above callback... 
 */

inline char* StrFormat(const std::string fmt_str, ...) {
    int final_n, n = ((int)fmt_str.size()) * 2; /* Reserve two times as much as the length of the fmt_str */
    std::unique_ptr<char[]> formatted;
    va_list ap;
    while(1) {
        formatted.reset(new char[n]); /* Wrap the plain char array into the unique_ptr */
        strcpy(&formatted[0], fmt_str.c_str());
        va_start(ap, fmt_str);
        final_n = vsnprintf(&formatted[0], n, fmt_str.c_str(), ap);
        va_end(ap);
        if (final_n < 0 || final_n >= n)
            n += abs(final_n - n + 1);
        else
            break;
    }
    return formatted.get();
}

}// End namespace

#endif
