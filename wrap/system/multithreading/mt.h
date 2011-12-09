#ifndef MT_MT_H
#define MT_MT_H


#ifdef QT_CORE_LIB
#include <QSemaphore>

tyepdef QSemaphore mt::Semaphore;

#else

#include "base.h"
#include "mutex.h"
#include "rw_lock.h"
#include "semaphore.h"
#include "thread.h"
#include "scoped_mutex_lock.h"
#include "scoped_read_lock.h"
#include "scoped_write_lock.h"

#endif

#endif // MT_MT_H
