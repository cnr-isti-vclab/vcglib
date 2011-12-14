#ifndef _ATOMIC_INT_GENERIC_H

#define _ATOMIC_INT_GENERIC_H

#include "mt.h"

namespace mt{

class atomicInt
{
public:
	atomicInt()
	{
		value = 0;
	}
	atomicInt( int value )
	{
		value = value;
	}

	// atomic API

	/**
	Reads the current value of this QAtomicInt and then adds valueToAdd
	to the current value, returning the original value.
	*/
	inline int fetchAndAddAcquire( int valueToAdd )
	{
		mutexlocker lock(m);
		int originalValue = value;
		value += valueToAdd;
		return originalValue;
	}

	/**
	Atomically increments the value of this atomicInt. 
	Returns true if the new value is non-zero, false otherwise.*/
	inline bool ref()
	{
		mutexlocker lock(m);
		value++;
		return value == 0;
	}

	/*
	Atomically decrements the value of this QAtomicInt. 
	Returns true if the new value is non-zero, false otherwise.*/
	inline bool deref()
	{
		mutexlocker lock(m);
		value--;
		return value == 0;
	}

	inline bool testAndSetOrdered(int expectedValue, int newValue)
	{
		mutexlocker lock(m);
		if (value == expectedValue) {
			 value = newValue;
			 return true;
		 }
		return false;
	}

    // Non-atomic API
    inline bool operator==(int value) const
    {
        return value == value;
    }

    inline bool operator!=(int value) const
    {
        return value != value;
    }

    inline bool operator!() const
    {
        return value == 0;
    }

    inline operator int() const
    {
        return value;
    }
	   
	inline atomicInt &operator=(int value)
    {
        value = value;
        return *this;
    }

	inline bool operator>(int value) const
	{
		return value > value;
	}

	inline bool operator<(int value) const
	{
		return value < value;
	}

private:
	volatile int value;
	mutex m;

};

}//namespace

#endif

