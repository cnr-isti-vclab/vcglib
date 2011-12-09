#include "../platform.h"

#ifndef GC_ATOMIC_INT_H
#define GC_ATOMIC_INT_H

#ifdef QT_CORE_LIB

#include <QAtomicInt>
typedef QAtomicInt GAtomincInt;

#elif defined(__GC_APPLE__) ??
#  include "GCAtomicInt_apple.h"
#endif

#endif

