//*******************************************************************
// L3DChronometer.h -- Header for L3DChronometer.cpp
// Copyright (c) 2010 Jose Maria Noguera
// Aug 10, 2010
//
// Jose Maria Noguera Rozua http://wwwdi.ujaen.es/~jnoguera/
//
//*******************************************************************


#ifndef _L3D_CHRONOMETER_H_
#define _L3D_CHRONOMETER_H_

#include "../../L3DPlatform.h"

#ifdef DEVICE_X86_WIN
	#include <windows.h>
#else
	#include <sys/time.h>
#endif


#define TIMER_STOPPED 0
#define TIMER_RUNNING 1
#define TIMER_PAUSED 2

/** @brief Chronometer used to measure performance of algoritms */

class L3DChronometer{

public:
    /** Initializes the Chronometer */
    L3DChronometer();
    /** Starts to measure time */
    void start ();
    /** Pause without losing the measured time */
    void pause ();
    /** Resume a previously paused Chronometer */
    void resume ();
    /** Stops the Chronometer */
    void stop ();
    float elapsed() { stop(); float r = timeMSecs(); resume(); return r }

    /** Gets time in seconds. The Chronometer must be stop */
    float timeSecs ();
    /** Gets time in seconds. The Chronometer must be stop */
    float timeMSecs ();
    /** Gets time in seconds. The Chronometer must be stop*/
    float timeUSecs ();

private:

#ifdef DEVICE_X86_WIN
  /* Initial time */
  LARGE_INTEGER tini;
  /* Seconds accumulator */
  LARGE_INTEGER tacum;
  /* Frequency */
  LARGE_INTEGER freq;
#else
  /** Initial time */
  timeval tini;
  /* Seconds accumulator */
  unsigned sec;
  /* MicroSeconds accumulator */
  long usec;
#endif

  /* Chronometer state: stopped, running, paused */
  int state;
  /* Updates time accumulators */
  void updateTime ();

};

#endif
