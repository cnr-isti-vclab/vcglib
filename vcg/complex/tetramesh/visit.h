

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
Revision 1.1  2006/12/03 14:55:44  ganovelli
created


****************************************************************************/


#ifndef __VCGLIB_TETRASUBSET
#define __VCGLIB_TETRASUBSET

#include <vector>
namespace vcg {
namespace tetra {

/** \addtogroup tetramesh */
/*@{*/

	/// assumes TTTopology has been computed
	template <class TetraPtrContainer>
		void Component(TetraPtrContainer & src, TetraPtrContainer & conn_com){

			typename TetraPtrContainer::iterator ti;
			typedef typename TetraPtrContainer::value_type TetraPointer;

			for(ti = src.begin(); ti != src.end(); ++ti)
				(*ti)->SetS();

			while(!src.empty()){
				TetraPointer tp = src.back();
				src.pop_back();
				conn_com.push_back(tp);
				for(unsigned int i = 0; i < 4; ++i)
					if(!tp->TTp(i)->IsD())
						if(!tp->TTp(i)->IsS()){
							tp->TTp(i)->SetS();
							src.push_back(tp->TTp(i));
						}
			}
			for(ti = conn_com.begin(); ti != conn_com.end(); ++ti)
				(*ti)->ClearS();
		}

/*@}*/
}	// End namespace
}	// End namespace


#endif

