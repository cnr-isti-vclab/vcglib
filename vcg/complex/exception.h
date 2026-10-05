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
#ifndef __VCG_EXCEPTION_H
#define __VCG_EXCEPTION_H

#include <stdexcept>

namespace vcg
{
// Every exception keeps its detail in what(), prefixed by its category ("Mesh does not
// satisfy precondition: There are faces with zero area"), so that a caller can show it.
// They used to print the detail to stdout and return only the fixed category from what().
class MissingComponentException : public std::runtime_error
{
public:
  MissingComponentException(const std::string &err):std::runtime_error("Missing component: " + err) {}
};

class MissingCompactnessException : public std::runtime_error
{
public:
  MissingCompactnessException(const std::string &err):std::runtime_error("Lack of compactness: " + err) {}
};

class MissingTriangularRequirementException : public std::runtime_error
{
public:
  MissingTriangularRequirementException(const std::string &err):std::runtime_error("Mesh has to be composed by triangle and not polygons: " + err) {}
};

class MissingPolygonalRequirementException : public std::runtime_error
{
public:
  MissingPolygonalRequirementException(const std::string &err):std::runtime_error("Mesh has to be composed by polygonal faces (not plain triangles): " + err) {}
};

class MissingTetrahedralRequirementException : public std::runtime_error
{
public:
  MissingTetrahedralRequirementException(const std::string &err):std::runtime_error("Mesh has to be composed by tetrahedras: " + err) {}
};

class MissingPreconditionException : public std::runtime_error
{
public:
  MissingPreconditionException(const std::string &err):std::runtime_error("Mesh does not satisfy precondition: " + err) {}
};

} // end namespace vcg
#endif // EXCEPTION_H
