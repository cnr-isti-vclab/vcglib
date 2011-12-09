#ifndef GC_ATOMIC_INT_APPLE_H

#define GC_ATOMIC_INT_APPLE_H

#include <libkern/OSAtomic.h>


//http://developer.apple.com/library/mac/#documentation/Darwin/Reference/KernelIOKitFramework/OSAtomic_h/index.html

class GCAtomicInt
{
public:
	GCAtomicInt()
	{
		_q_value = 0;
	}
	GCAtomicInt( int value )
	{
		_q_value = value;
	}

	// atomic API

	/**
	Reads the current value of this QAtomicInt and then adds valueToAdd
	to the current value, returning the original value.

	Unfortunately, MacOSX does not provide with fetch-and-add functions,
	only add-and-fetch. Therefore, we have to simulate them.

	Implementation based on SDL:
	//http://lists.libsdl.org/pipermail/commits-libsdl.org/2011-January/003568.html
	*/
	inline int fetchAndAddAcquire( int valueToAdd )
	{
		 //T *originalValue = currentValue;
		 //currentValue += valueToAdd;
		 //return originalValue;

		int originalValue;
		do { 
			originalValue = _q_value;
		} while (!OSAtomicCompareAndSwap32Barrier(originalValue, originalValue+valueToAdd, &_q_value));
		return originalValue;
	}

	/**
	Atomically increments the value of this GCAtomicInt. 
	Returns true if the new value is non-zero, false otherwise.*/
	inline bool ref()
	{
		return OSAtomicIncrement32Barrier(&_q_value) != 0;
	}

	/*
	Atomically decrements the value of this QAtomicInt. 
	Returns true if the new value is non-zero, false otherwise.*/
	inline bool deref()
	{
		return OSAtomicDecrement32Barrier(&_q_value) != 0;
	}

	inline bool testAndSetOrdered(int expectedValue, int newValue)
	{
		//if (currentValue == expectedValue) {
		//	 currentValue = newValue;
		//	 return true;
		// }
		//return false;

                return OSAtomicCompareAndSwap32Barrier(expectedValue, newValue, &_q_value);
	}

    // Non-atomic API
    inline bool operator==(int value) const
    {
        return _q_value == value;
    }

    inline bool operator!=(int value) const
    {
        return _q_value != value;
    }

    inline bool operator!() const
    {
        return _q_value == 0;
    }

    inline operator int() const
    {
        return _q_value;
    }
	   
	inline GCAtomicInt &operator=(int value)
    {
        _q_value = value;
        return *this;
    }

	inline bool operator>(int value) const
	{
		return _q_value > value;
	}

	inline bool operator<(int value) const
	{
		return _q_value < value;
	}

private:
	volatile int _q_value;

};


#endif

