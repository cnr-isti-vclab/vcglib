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
Revision 1.1  2004/03/31 12:46:53  cignoni
First working version!


****************************************************************************/
#ifndef __VCG_VERTEX_PLUS_COMPONENT_OPT
#define __VCG_VERTEX_PLUS_COMPONENT_OPT

#include <vcg/component.h>
#include <vcg/traced_vector.h>


namespace vcg {
  namespace vert {
/*
Some naming Rules
All the Components that can be added to a vertex should be defined in the namespace vert:

*/

/*------------------------- COORD -----------------------------------------*/ 

template <class A, class T> class CoordOpt: public T {
public:
  typedef A CoordType;
  typedef typename CoordType::ScalarType      ScalarType;
	CoordType &P() { return CAT< TVector<VertType>,CoordType>::Get((VertType*)this); }
  CoordType &UberP() { return CAT< TVector<VertType>,CoordType>::Get((VertType*)this); }
};
template <class T> class Coord3fOpt: public CoordOpt<vcg::Point3f, T> {};
template <class T> class Coord3dOpt: public CoordOpt<vcg::Point3d, T> {};


/*-------------------------- NORMAL ----------------------------------------*/ 
/*-------------------------- TEXTURE ----------------------------------------*/ 

/*------------------------- FLAGS -----------------------------------------*/ 

/*-------------------------- COLOR ----------------------------------*/ 


/*-------------------------- Quality  ----------------------------------*/ 


/*----------------------------- VFADJ ------------------------------*/ 


  } // end namespace vert
}// end namespace vcg
#endif