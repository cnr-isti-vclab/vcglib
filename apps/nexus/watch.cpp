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
Revision 1.3  2004/10/21 13:40:16  ponchio
Debugging.

Revision 1.2  2004/10/21 12:14:02  ponchio
Support for mfile (>4Gb)

Revision 1.1  2004/10/19 17:20:24  ponchio
renamed

Revision 1.5  2004/10/19 17:16:52  ponchio
Changed interface

Revision 1.4  2004/07/30 12:44:14  ponchio
#ifdef corrected

Revision 1.3  2004/07/20 14:03:47  ponchio
Changed interface.

Revision 1.2  2004/07/05 15:49:39  ponchio
Windows (DevCpp, mingw) port.

Revision 1.1  2004/07/01 21:38:30  ponchio
First draft created.


****************************************************************************/
#include "watch.h"
#include <stdio.h>

#ifdef WIN32
Watch::Watch(): elapsed(0) {
    QueryPerformanceFrequency(&freq);
}


void Watch::Start(void) {
  QueryPerformanceCounter(&tstart);
  elapsed = 0;
}

double Watch::Pause() {
  QueryPerformanceCounter(&tend);         
  elapsed += Diff();
  return (double)elapsed;
}               

void Watch::Continue() {
  QueryPerformanceCounter(&tstart);
}

double Watch::Time() {
  QueryPerformanceCounter(&tend);         
  return (double)(elapsed + Diff());
}

double Watch::Diff() {  
  return ((double)tend.QuadPart -
	  (double)tstart.QuadPart)/
    ((double)freq.QuadPart);
}         

#else

Watch::Watch(): elapsed() {}

void Watch::Start() {
   gettimeofday(&tstart, &tz); 
   elapsed = 0;
}

double Watch::Pause() {
  gettimeofday(&tend, &tz); 
  elapsed += Diff();
  return (double)elapsed;
}

void Watch::Continue() {
  gettimeofday(&tstart, &tz); 
}

double Watch::Time() {
  gettimeofday(&tend, &tz); 
  return (double)(elapsed + Diff());
}
 
double Watch::Diff() {
  double t1 =  (double)tstart.tv_sec + (double)tstart.tv_usec/(1000*1000);
  double t2 =  (double)tend.tv_sec + (double)tend.tv_usec/(1000*1000);
  return t2 - t1;   
}
#endif

void Watch::Reset() {
  elapsed = 0;
}

int Watch::Usec() {
#ifdef WIN32
  return 0;
#else  
  struct timeval ttime;
  gettimeofday(&ttime, &tz);
  return ttime.tv_usec;
#endif
}

void Report::Init(unsigned int t, double inter) {
  watch.Start();
  tot = t;
  last = 0;
  interval = inter;
}

void Report::Step(unsigned int count) {
  if(count == 0) return;
  double now = watch.Time();
  if(now - last < interval) return;
  //estimate final time
  double tot_time = now * tot/(double)count;
  printf("%d/%d\telapsed: %.1f\tremaining: %.1f\ttotal: %.1f\n",
	 count, tot, now, tot_time - now, tot_time);
  last = now;
}

void Report::Finish() {
  double now = watch.Time();
  printf("Tot: %.1f\tN: %d\tN/sec: %.1f\n",
	 now, tot, tot/now);
}
