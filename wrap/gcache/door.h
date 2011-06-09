/****************************************************************************
* GCache                                                                    *
* Author: Federico Ponchio                                                  *
*                                                                           *
* Copyright(C) 2011                                                         *
* Visual Computing Lab                                                      *
* ISTI - Italian National Research Council                                  *
*                                                                           *
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


#ifndef CACHE_DOOR_H
#define CACHE_DOOR_H

#include <QDebug>


/*
//a door needs to be open for the thread to continue,
//if it is open the thread enter and closes the door
//this mess is to avoid [if(!open.available()) open.release(1)]
#include <QSemaphore>
class QDoor {
 private:
  QSemaphore _open;
  QSemaphore _close;
 public:
  QDoor(): _open(0), _close(1) {} //this means closed
  void open() {
    if(_close.tryAcquire(1)) //check it is not open
      _open.release(1); //open
  }
  void close() {
    if(_open.tryAcquire(1)) //check not already cloed
      _close.release(1);
  }
  void enter(bool close = false) {
    _open.acquire(1);
    if(close)
      _close.release(1); //and close door behind
    else
      _open.release(1); //and leave door opened
  }
  bool isOpen() { return _open.available() == 1; }
};

*/
#include <QMutex>
#include <QWaitCondition>

/**
  A wait condition class that works as a door.
  Should check if the semaphore version is faster.
*/

class QDoor {
 public:

  QDoor(void) : doorOpen(false), waiting(false) {}

  ///opens the door. Threads trying to enter will be awakened
  void open(void) {
    m.lock();
    doorOpen = true;
    m.unlock();
    c.wakeAll();
  }

  ///attempt to enter the door. if the door is closed the thread will wait until the door is opened.
  /** if close is true, the door will be closed after the thread is awakened, this allows to 
     have only one thread entering the door each time open() is called */
  void enter(bool close = false) {
    m.lock();
    waiting = true;
    while (!doorOpen)
      c.wait(&(m));
    
    if(close) 
      doorOpen = false;
    waiting = false;
    m.unlock();
  }
  bool isWaiting() {
    m.lock();
    bool w = waiting;
    m.unlock();
    return w;
  }
  void lock() { //prevend door opening and entering
    m.lock();
  }
  void unlock() { //reverse effect of lock
    m.unlock();
  }
 private:
  QMutex m;
  QWaitCondition c;
  bool doorOpen;
  bool waiting;
};


#endif
