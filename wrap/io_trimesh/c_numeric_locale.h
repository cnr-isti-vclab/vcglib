/****************************************************************************
* VCGLib                                                            o o     *
* Visual and Computer Graphics Library                            o     o   *
*                                                                _   O  _   *
* Copyright(C) 2004-2026                                           \/)\/    *
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

#ifndef __VCGLIB_C_NUMERIC_LOCALE
#define __VCGLIB_C_NUMERIC_LOCALE

#include <clocale>
#if defined(_WIN32)
#include <string>
#else
#include <locale.h>
#if defined(__APPLE__)
#include <xlocale.h>
#endif
#endif

namespace vcg {
namespace tri {
namespace io {

/** Makes the calling thread use the "C" numeric locale for as long as the object lives.

The text importers and exporters read and write numbers with atof, fscanf and fprintf,
which follow the process's C locale. A host that sets it from the environment (Qt does,
with setlocale(LC_ALL, "") on Unix) makes them read and write a comma as the decimal
separator in an Italian or German session: "0.5" reads as 0 and files are written as
"0,5", which no other program reads. File formats always use a dot.

Declare one at the top of an importer or exporter. Only the calling thread is switched,
so other threads of the host keep their locale, and the previous state is restored on
exit, so the host never sees a change.
*/
class CNumericLocale
{
public:
#if defined(_WIN32)
  CNumericLocale() : prevMode(_configthreadlocale(_ENABLE_PER_THREAD_LOCALE))
  {
    const char *cur = std::setlocale(LC_NUMERIC, nullptr);
    prevLocale = cur ? cur : "C";
    std::setlocale(LC_NUMERIC, "C");
  }
  ~CNumericLocale()
  {
    std::setlocale(LC_NUMERIC, prevLocale.c_str());
    _configthreadlocale(prevMode);
  }
#else
  CNumericLocale() : prev(uselocale(cLocale())) {}
  ~CNumericLocale() { uselocale(prev); }
#endif
  CNumericLocale(const CNumericLocale &) = delete;
  CNumericLocale &operator=(const CNumericLocale &) = delete;

private:
#if defined(_WIN32)
  int prevMode;
  std::string prevLocale;
#else
  // Created once and never freed: it lives as long as the process. Its other categories
  // are "C" too, which nothing in the importers depends on.
  static locale_t cLocale()
  {
    static const locale_t loc = newlocale(LC_NUMERIC_MASK, "C", (locale_t)0);
    return loc;
  }
  locale_t prev;
#endif
};

} // end namespace io
} // end namespace tri
} // end namespace vcg

#endif
