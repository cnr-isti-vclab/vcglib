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
Revision 1.1  2004/07/01 21:38:30  ponchio
First draft created.


****************************************************************************/
#include "stopwatch.h"

StopWatch::StopWatch(): elapsed(0) {}

#ifdef WIN32
void StopWatch::Start(void) {
  static int first = 1;
  if(first) {
    QueryPerformanceFrequency(&freq);
    first = 0;
  }
  QueryPerformanceCounter(&tstart);
}

void StopWatch::Stop() {
  QueryPerformanceCounter(&tend);         
  elapsed += Diff();
}               

double StopWatch::Diff() {  
  return ((double)tend.QuadPart -
	  (double)tstart.QuadPart)/
    ((double)freq.QuadPart);
}         

#else

void StopWatch::Start() {
   gettimeofday(&tstart, &tz); 
}

void StopWatch::Stop() {
  gettimeofday(&tend, &tz); 
  elapsed += Diff();
}
 
double StopWatch::Diff() {
  double t1 =  (double)tstart.tv_sec + (double)tstart.tv_usec/(1000*1000);
  double t2 =  (double)tend.tv_sec + (double)tend.tv_usec/(1000*1000);
  return t2 - t1;   
}
#endif

void StopWatch::Reset() {
  elapsed = 0;
}

double StopWatch::Elapsed() {
  return elapsed;
}

int StopWatch::Usec() {
#ifdef WIN32
  return 0;
#else  
  struct timeval ttime;
  gettimeofday(&ttime, &tz);
  return ttime.tv_usec;
#endif
}
