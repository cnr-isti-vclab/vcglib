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


****************************************************************************/
#ifndef _CMESH_H
#define _CMESH_H

#pragma warning(disable:4786 4804 4666) 

#include <math.h>
#include <vcg/simplex/vertex/with/vcvq.h>
#include <vcg/simplex/face/with/rtfcfmfn.h>
#include <vcg/space/index/grid_static_ptr.h>
#include <vcg/complex/trimesh/base.h>



using namespace vcg;
using namespace std;

// Vertex, Face, Mesh and Grid definitions.
class MyEdge;
class CFace;
class CVertex   : public VertexVCVQ< double,MyEdge,CFace > {};
class CFace     : public FaceRTFCFMFN< CVertex,MyEdge,CFace > {};
class CMesh     : public tri::TriMesh< vector<CVertex>, vector<CFace> > {};

#endif
