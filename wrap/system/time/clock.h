//*******************************************************************
// L3DTimer.h -- Header for L3DTimer.cpp
// Copyright (c) 2010 Jose Maria Noguera
// Aug 10, 2010
//
// Jose Maria Noguera Rozua http://wwwdi.ujaen.es/~jnoguera/
//
//*******************************************************************


#ifndef _L3D_TIMER_
#define _L3D_TIMER_

#include "../../L3DPlatform.h"

class L3DTimer
{

public:
	L3DTimer();
	~L3DTimer(){;}

	/**
	This function retrieves the system time, in milliseconds. The system time is the time elapsed since the first call
	to this function.
	*/
	static int GetCurrentSystemTime(void);

private:
/*
	static int base;
	static bool initialized;
	static int basetime;
*/
	#ifdef DEVICE_X86_WIN
		static int Win_GetCurrentSystemTime(void);
	#else
		static int Linux_GetCurrentSystemTime(void);
	#endif
};

#endif

