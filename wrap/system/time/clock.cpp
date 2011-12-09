#include "L3DTimer.h"


#ifdef DEVICE_X86_WIN
	#include <Windows.h>
#else
	#include <sys/time.h>
#endif


#include <time.h>




L3DTimer::L3DTimer()
{/*
	base = 0;
	initialized = false;
	basetime = 0;*/
}

//---------------------------------------------------------------------------

int L3DTimer::GetCurrentSystemTime(void)
{
	int iReturnValue = -1;

	#ifdef DEVICE_X86_WIN
		iReturnValue = Win_GetCurrentSystemTime();
	#else
		iReturnValue = Linux_GetCurrentSystemTime();
	#endif

	return iReturnValue;
}

//---------------------------------------------------------------------------

#ifdef DEVICE_X86_WIN

int L3DTimer::Win_GetCurrentSystemTime(void)
{
	int curtime;
	static int base;
	static bool initialized = false;
	if(!initialized)
	{
		base = timeGetTime() & 0xffff0000;
		initialized = true;
	}
	curtime = timeGetTime() - base;
	return curtime;      
}

//---------------------------------------------------------------------------

#else

//---------------------------------------------------------------------------

int L3DTimer::Linux_GetCurrentSystemTime(void)
{
	struct timeval tp;
	struct timezone tzp;
	static int basetime;
	gettimeofday(&tp, &tzp);
	if(!basetime)
	{
		basetime = tp.tv_sec;
		return tp.tv_usec / 1000;
	}
	return (tp.tv_sec - basetime) * 1000 + tp.tv_usec / 1000;
}
#endif

//---------------------------------------------------------------------------
