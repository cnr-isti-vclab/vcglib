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
Revision 1.4  2004/11/28 04:10:59  ponchio
winsockapi include problem

Revision 1.3  2004/10/21 13:40:16  ponchio
Debugging.

Revision 1.2  2004/10/21 12:14:02  ponchio
Support for mfile (>4Gb)

Revision 1.1  2004/10/19 17:20:24  ponchio
renamed

Revision 1.4  2004/10/19 17:16:53  ponchio
Changed interface

Revision 1.3  2004/07/20 14:03:47  ponchio
Changed interface.

Revision 1.2  2004/07/05 15:49:39  ponchio
Windows (DevCpp, mingw) port.

Revision 1.1  2004/07/01 21:38:30  ponchio
First draft created.


****************************************************************************/

#ifdef WIN32
#ifndef _WINDOWS_
#define _WINSOCKAPI_
#include <windows.h>
#endif
#else
#include <sys/time.h>
#include <unistd.h>     
#endif

class Watch {
public:
  Watch();
  void Start();
  float Pause();
  void Continue();
  void Reset();
  float Time();
  int Usec();
private:
  double Diff();

#ifdef WIN32
  LARGE_INTEGER tstart, tend;
  LARGE_INTEGER freq;
#else  
  struct timeval tstart, tend;  
  struct timezone tz;   
#endif
  double elapsed;
};

class Report {
 public:
  Report(unsigned int tot = 1, float inter = 30.0f) { Init(tot, inter); }
  void Init(unsigned int tot, float inter = 30.0f);
  void Step(unsigned int count);
  void Finish();
 private:
  Watch watch;
  int tot;
  float last;
  float interval;
};
