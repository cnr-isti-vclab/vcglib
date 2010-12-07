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
Revision 1.6  2005/01/26 22:45:34  cignoni
Release 4.04
final updates for gcc compiling issues

Revision 1.5  2004/05/14 00:32:36  ganovelli
just color and quality on the vertex



****************************************************************************/
#ifndef _CMESH_H
#define _CMESH_H


#include <math.h>
#include <vcg/space/index/grid_static_ptr.h>
#include <vcg/complex/trimesh/base.h>

// Vertex, Face, Mesh and Grid definitions.
#include <vcg/simplex/vertex/base.h>
#include <vcg/simplex/vertex/component.h>
#include <vcg/simplex/face/base.h>
#include <vcg/simplex/face/component.h>
#include <vcg/simplex/face/component_rt.h>
class CFace;
class CVertex;
struct UsedTypes:public vcg::UsedTypes< vcg::Use<CFace>::AsFaceType, vcg::Use<CVertex>::AsVertexType>{};
class CVertex   : public vcg::Vertex<UsedTypes,vcg::vertex::Coord3d,vcg::vertex::Qualityf,vcg::vertex::Normal3d,vcg::vertex::Color4b,vcg::vertex::BitFlags> {};
class CFace     : public vcg::Face< UsedTypes,vcg::face::VertexRef, vcg::face::Normal3d, vcg::face::EdgePlane,vcg::face::Color4b,vcg::face::Mark,vcg::face::BitFlags> {};
class CMesh     : public vcg::tri::TriMesh< std::vector<CVertex>, std::vector<CFace> > {};

#endif
