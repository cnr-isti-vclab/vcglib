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


#ifndef __VCGLIB_SIMPLE__
#define __VCGLIB_SIMPLE__

#include <vector>

namespace vcg {

template <class STL_CONT, class ATTR_TYPE>
class SimpleTempData{
public:
typedef typename ATTR_TYPE attr_type;

STL_CONT& c;
std::vector<ATTR_TYPE> data;

SimpleTempData(STL_CONT  &_c):c(_c){};

// access to data
ATTR_TYPE & operator[](const STL_CONT::value_type *v){return data[v-&*c.begin()];}
ATTR_TYPE & operator[](const int & i){return data[i];}

// start temporary attribute
void Start(){data.reserve(c.capacity());data.resize(c.size());}

// stop temporary attribute
void Stop(){data.clear();}

// update temproary data size 
bool UpdateSize(){
		if(data.size() != c.size())
			{
				data.resize(c.size());
				return false;
			}
		return true;
	}
};

}; // end namespace vcg

#endif