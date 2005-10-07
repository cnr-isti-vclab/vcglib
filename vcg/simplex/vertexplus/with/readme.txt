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

This folders contains most common VERTEX configuration files.

The name of the file specifies the members that are added to the vertex
class. The name is a sequence of letter pairs, in strict alphabetical order. The
possible admitted letters pairs are:

Adjacency Info

AF - Vertex-Face adjacency
AE - Vertex-Edge adjacency
AT - Vertex-Tetra adjacency

Per-Vertex Data
VC - Color
VN - Normal
VM - Incremental Mark
VQ - Quality
VT - Texture Coords

E.g. 

#include<vcg/simplex/vertex/with/afvcvnvq.h> 

generate a type 

VertexAFVCVQ<VScalarType,FaceType> 

That can store V-F adjacency, color, normal and quality.
  

