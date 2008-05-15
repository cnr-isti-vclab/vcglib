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
Revision 1.5  2007/02/02 00:01:54  tarini
overloaded operator "[]" (once more) to make it possible to index the temp. structure with an iterator

Revision 1.4  2005/07/11 13:12:34  cignoni
small gcc-related compiling issues (typenames,ending cr, initialization order)

Revision 1.3  2004/12/11 15:37:47  ganovelli
added one more  [], now it is polymorphic, added typenames

Revision 1.2  2004/03/31 22:36:44  ganovelli
First Working Release (with this comment)


****************************************************************************/

#ifndef __VCGLIB_SIMPLE__
#define __VCGLIB_SIMPLE__

#include <vector>
 
namespace vcg {

template <class STL_CONT>
class SimpleTempDataBase{
public:
	virtual void Enable()= 0;
	virtual void Disable()= 0;
	virtual void Resize(const int & sz) = 0;
	virtual void Reorder(std::vector<size_t> & newVertIndex)=0;
};

template <class STL_CONT, class ATTR_TYPE>
class SimpleTempData:public SimpleTempDataBase<STL_CONT>{

public:
typedef SimpleTempData<STL_CONT,ATTR_TYPE> SimpTempDataType;
typedef ATTR_TYPE AttrType;

STL_CONT& c;
std::vector<ATTR_TYPE> data;

SimpleTempData(STL_CONT  &_c):c(_c){data.reserve(c.capacity());data.resize(c.size());};
SimpleTempData(STL_CONT  &_c,ATTR_TYPE val):c(_c){
	data.reserve(c.capacity());data.resize(c.size());
	typename std::vector<ATTR_TYPE>::iterator i;
	for(i = data.begin(); i!= data.end(); ++i)
	*i = val;
};

~SimpleTempData(){data.clear();}

// access to data
ATTR_TYPE & operator[](const typename STL_CONT::value_type & v){return data[&v-&*c.begin()];}
ATTR_TYPE & operator[](const typename STL_CONT::value_type * v){return data[v-&*c.begin()];}
ATTR_TYPE & operator[](const typename STL_CONT::iterator & cont){return data[&(*cont)-&*c.begin()];}
ATTR_TYPE & operator[](const int & i){return data[i];}

// enable temporary attribute: REMOVED!! The memory is allocated by the contructor
void Enable(){assert(0);}

// disable and initialize temporary attribute: REMOVED!! The memory is allocated by the contructor
void Enable(ATTR_TYPE val){assert(0);}

// stop temporary attribute: REMOVED!! The memory is freed by the distructor
void Disable(){assert(0);}

// update temporary data size 
bool UpdateSize(){
		if(data.size() != c.size())
			{
				data.resize(c.size());
				return false;
			}
		return true;
	}

void Resize(const int & sz){
	data.resize(sz);
}

void Reorder(std::vector<size_t> & newVertIndex){
	for(unsigned int i = 0 ; i < data.size(); ++i){
		if( newVertIndex[i] != std::numeric_limits<size_t>::max())
			data[newVertIndex[i]] = data[i];
	}
}

};



}; // end namespace vcg

#endif
