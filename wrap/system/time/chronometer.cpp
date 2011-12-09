
#include "L3DChronometer.h"


/******************************************************************************

Implementacion de metodos multiplataforma

*******************************************************************************/

// -----------------------------------------------------------------------------

void L3DChronometer::stop ()
{
  if (state == TIMER_RUNNING || state == TIMER_PAUSED) {
    if (state == TIMER_RUNNING)
      updateTime ();
    state = TIMER_STOPPED;
  }
}

// -----------------------------------------------------------------------------

void L3DChronometer::pause ()
{
  if (state == TIMER_RUNNING) {
    updateTime ();
    state = TIMER_PAUSED;
  }
}

/******************************************************************************

Implementacion de metodos para Windows

*******************************************************************************/

#ifdef DEVICE_X86_WIN

// -----------------------------------------------------------------------------

L3DChronometer::L3DChronometer() { 
	 tacum.QuadPart = 0; 
	 QueryPerformanceFrequency(&freq); 
	 state = TIMER_STOPPED;
 }

// -----------------------------------------------------------------------------

float L3DChronometer::timeSecs () { 
	return (float) tacum.QuadPart / freq.QuadPart;
}

// -----------------------------------------------------------------------------

float L3DChronometer::timeMSecs () { 
	return (float) tacum.QuadPart / ((float)freq.QuadPart / 1000); 
}

// -----------------------------------------------------------------------------

float L3DChronometer::timeUSecs () { 
	return (float) tacum.QuadPart / ((float)freq.QuadPart / 1000000u); 
}

// -----------------------------------------------------------------------------

void L3DChronometer::updateTime ()
{
  LARGE_INTEGER tnow;

  QueryPerformanceCounter(&tnow);
  tacum.QuadPart += (tnow.QuadPart - tini.QuadPart);
  tini.QuadPart = tnow.QuadPart;
}

// -----------------------------------------------------------------------------

void L3DChronometer::start ()
{
  if (state == TIMER_STOPPED) {
    tacum.QuadPart = 0;
    QueryPerformanceCounter(&tini);
    state = TIMER_RUNNING;
  }
}

// -----------------------------------------------------------------------------

void L3DChronometer::resume ()
{
  if (state == TIMER_PAUSED) {
    QueryPerformanceCounter(&tini);
    state = TIMER_RUNNING;
  }
}

// -----------------------------------------------------------------------------



/******************************************************************************

Implementacion de metodos para sistemas operativos posix

*******************************************************************************/

#else

// -----------------------------------------------------------------------------

L3DChronometer::L3DChronometer() { 
	sec = 0;
	usec = 0;
	state = TIMER_STOPPED; 
}

// -----------------------------------------------------------------------------

float L3DChronometer::timeSecs () { 
	return (float) sec + (float) usec / 1000000.0f;
}

// -----------------------------------------------------------------------------

float L3DChronometer::timeMSecs () {
	return (float) sec * 1000.0f + (float) usec / 1000.0f;
}

// -----------------------------------------------------------------------------

float L3DChronometer::timeUSecs () {
	return (float) sec * 1000000.0f + (float) usec;
}

// -----------------------------------------------------------------------------

void L3DChronometer::updateTime ()
{
  timeval tnow;
  gettimeofday (&tnow, 0);
  
  sec += tnow.tv_sec - tini.tv_sec;
  usec += tnow.tv_usec - tini.tv_usec;
  if (usec < 0) {
    sec --;
    usec += 1000000u;
  }
}

// -----------------------------------------------------------------------------

void L3DChronometer::start () 
{ 
  if (state == TIMER_STOPPED) {
    sec = 0; usec = 0; 
    gettimeofday (&tini, 0); 
    state = TIMER_RUNNING;
  } 
}

// -----------------------------------------------------------------------------

void L3DChronometer::resume () 
{ 
  if (state == TIMER_PAUSED) { 
    gettimeofday (&tini, 0); 
    state = TIMER_RUNNING;
  }
}

// -----------------------------------------------------------------------------

#endif

