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

This folders contains most common FACE configuration files.

The name of the file specifies the members that are added to the vertex
class. The name is a sequence of letter pairs, in strict alphabetical order. The
possible admitted letters pairs are:

AF - Face-Face adjacency
AS - Shared Vertex-Face and Face-Face Adjacency
AV - Vertex-face adjacency

FC - Per-Face Color
FM - Per-Face Incremental Mark
FN - Per-Face Normal
FQ - Per-Face Quality
RT - Data for Optimized Point-Face Distance and Ray-Tracing Stuff
WC - Per-Wedge Color
WN - Per-Wedge Normal
WQ - Per-Wedge Quality
WT - Per-Wedge Texture Coords

E.g.

#include<vcg/simplex/vertex/with/affnwt.h>

generate a type

VertexAFFNWT<VertexType>

that can store F-F adjacency, Per-face normal and color and per-wedge texture coords.
