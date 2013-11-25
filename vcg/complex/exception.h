/****************************************************************************
* VCGLib                                                            o o     *
* Visual and Computer Graphics Library                            o     o   *
*                                                                _   O  _   *
* Copyright(C) 2004-2012                                           \/)\/    *
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

namespace vcg
{
class MissingComponentException : public std::runtime_error
{
public:
  MissingComponentException(const std::string &err):std::runtime_error(err)
  {
    std::cout << "Missing Component Exception -" << err << "- \n";
  }
    virtual const char *what() const throw ()
    {
      static char buf[128]="Missing Component";
      return buf;
    }
};
}

#endif // EXCEPTION_H
