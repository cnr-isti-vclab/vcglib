#ifndef __ML_ATOMIC_GUARD_H
#define __ML_ATOMIC_GUARD_H

class MLAtomicGuard
{
public:
	MLAtomicGuard(bool val)
		:_lock(QReadWriteLock::Recursive),_guard(val) {}
	
	~MLAtomicGuard() {}

	MLAtomicGuard& operator=(bool v)
	{
		QWriteLocker locker(&_lock);
		_guard = v;
		return *this;
	}

	bool operator==(bool v) const
	{
		QReadLocker locker(&_lock);
		return (_guard == v);
	}

	bool operator!=(bool v) const
	{
		QReadLocker locker(&_lock);
		return (_guard != v);
	}

private:
	bool _guard;
	mutable QReadWriteLock _lock;
};

#endif