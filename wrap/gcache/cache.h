#ifndef GCACHE_CACHE_H
#define GCACHE_CACHE_H

#include <iostream>
#include <limits.h>
#include <vector>
#include <list>

#include <QThread>
#include "provider.h"

using namespace std;
/* this cache system enforce the rule that the items in a cache are always in all the cache below */
/* two mechanism to remove tokens from the cache:
      1) set token count to something low
      2) set maximum number of tokens in the provider
*/

/** Cache virtual base class. You are required to implement the pure virtual functions get, drop and size.
*/

template <typename Token> class Transfer;

template <typename Token>
class Cache: public Provider<Token> {

public:
  ///true if this is the last cache (the one we use the data from)
  bool final;
  //if true the cache will exit at the first opportunity
  bool quit;
  ///keeps track of changes (if 1 then something was loaded or dropped
  QAtomicInt changed;
  ///callback for changed
  void (*callback)(void *data);

  ///data is fetched from here
  Provider<Token> *input;

  ///threads running over cache...
  std::vector<Transfer<Token> *> transfers;

protected:
  ///max space available
  quint64 s_max;
  ///current space used
  quint64 s_curr;

public:
  Cache(quint64 _capacity = INT_MAX):
    final(false), quit(false), changed(false), input(NULL), s_max(_capacity), s_curr(0) {}
  virtual ~Cache() {}

  void setInputCache(Provider<Token> *p) { input = p; }
  quint64 capacity() { return s_max; }
  quint64 size() { return s_curr; }
  void setCapacity(quint64 c) { s_max = c; }

  ///return true if the cache is waiting for priority to change
  bool isChanged() {
    bool r = changed.testAndSetOrdered(1, 0); //if changed is 1, r is true
    return r;
  }

  ///empty the cache. Make sure no resource is locked before calling this.
  /// Require pause or stop before. Ensure there no locked item
  void flush() {
    //std::vector<Token *> tokens;
    {
      for(int i = 0; i < this->heap.size(); i++) {
        Token *token = &(this->heap[i]);
        //tokens.push_back(token);
        s_curr -= drop(token);
        assert(!(token->count >= Token::LOCKED));
        if(final)
          token->count.testAndSetOrdered(Token::READY, Token::CACHE);
        input->heap.push(token);
      }
      this->heap.clear();
    }
    if(!s_curr == 0) {
      qDebug() << "Cache size after flush is not ZERO!";
      s_curr = 0;
    }
    //assert(s_curr == 0);

/*    {
      for(unsigned int i = 0; i < tokens.size(); i++) {
        input->heap.push(tokens[i]);
      }
    }*/
  }

  ///empty the cache. Make sure no resource is locked before calling this.
  /// Require pause or stop before. Ensure there no locked item
  template <class FUNCTOR> void flush(FUNCTOR functor) {
    std::vector<Token *> tokens;
    {
      int count = 0;
      QMutexLocker locker(&(this->heap_lock));
      for(int k = 0; k < this->heap.size(); k++) {
        Token *token = &this->heap[k];
        if(functor(token)) { //drop it
          tokens.push_back(token);
          s_curr -= drop(token);
          assert(token->count < Token::LOCKED);
          if(final)
            token->count.testAndSetOrdered(Token::READY, Token::CACHE);
        } else
          this->heap.at(count++) = token;
      }
      this->heap.resize(count);
      this->heap_dirty = true;
    }
    {
      QMutexLocker locker(&(input->heap_lock));
      for(unsigned int i = 0; i < tokens.size(); i++) {
        input->heap.push(tokens[i]);
      }
    }
  }

  virtual void abort() {}
protected:
  ///return the space used in the cache by the loaded resource
  virtual int size(Token *token) = 0;
  ///returns amount of space used in cache -1 for failed transfer
  virtual int get(Token *token) = 0;
  ///return amount removed
  virtual int drop(Token *token) = 0;
  ///make sure the get function do not access token after abort is returned.




  ///called in as first thing in run()
  virtual void begin() {}
  ///called in as last thing in run()
  virtual void end() {}

  ///[should be protected]
  void run() {
    assert(input);
    /* basic operation of the cache:
       1) transfer first element of input_cache if
          cache has room OR first element in input as higher priority of last element
       2) make room until eliminating an element would leave space. */
    begin();
    while(!this->quit) {
      input->check_queue.enter();     //wait for cache below to load something or priorities to change
      if(this->quit) break;

      if(unload() || load()) {
        changed.testAndSetOrdered(0, 1);  //if not changed, set as changed
        input->check_queue.open();        //we signal ourselves to check again
      }
      input->check_queue.leave();
    }
    this->quit = false;                   //in case someone wants to restart;
    end();
  }



  /** Checks wether we need to make room in the cache because of:
     size() - sizeof(lowest priority item) > capacity()
  **/
  bool unload() {
    Token *remove = NULL;
    //make room int the cache checking that:
    //1 we need to make room (capacity < current)
    if(size() > capacity()) {

      QMutexLocker locker(&(this->heap_lock));

      //2 we have some element not in the upper caches (heap.size()  > 0
      if(this->heap.size()) {
        Token &last = this->heap.min();
        int itemsize = size(&last);

        //3 after removing the item, we are still full (avoids bouncing items)
        if(size() - itemsize > capacity()) {

          //4 item to remove is not locked. (only in last cache. you can't lock object otherwise)
          if(!final) { //not final we can drop when we want
            remove = this->heap.popMin();
          } else {
            last.count.testAndSetOrdered(Token::READY, Token::CACHE);
            if(last.count <= Token::CACHE) { //was not locked and now can't be locked, remove it.
              remove = this->heap.popMin();
            } else { //last item is locked need to reorder stack
              remove = this->heap.popMin();
              this->heap.push(remove);
              return true;
            }
          }
        }
      }
    }

    if(remove) {
      {
        QMutexLocker input_locker(&(input->heap_lock));
        int size = drop(remove);
        assert(size >= 0);
        s_curr -= size;
        input->heap.push(remove);
      }
      return true;
    }
    return false;
  }

  ///should be protected
  bool load() {
    Token *insert = NULL;
    Token *last = NULL;              //we want to lock only one heap at once to avoid deadlocks.

    /* check wether we have room (curr < capacity) or heap is empty.
       empty heap is bad: we cannot drop anything to make room, and cache above has nothing to get.
       this should not happen if we set correct cache sizes, but if it happens.... */
    {
      QMutexLocker locker(&(this->heap_lock));
      this->rebuild();
      if(size() > capacity() && this->heap.size() > 0) {
        last = &(this->heap.min()); //no room, set last so we might check for a swap.
      }
    }

    {
      QMutexLocker input_locker(&(input->heap_lock));
      input->rebuild();                                  //if dirty rebuild
      if(input->heap.size()) {                           //we need something in input to tranfer.
        Token &first = input->heap.max();
        if(first.count > Token::REMOVE &&
            (!last || last->priority < first.priority)) { //if !last we already decided we want a transfer., otherwise check for a swap
          insert = input->heap.popMax();                 //remove item from heap, while we transfer it.
        }
      }
    }

    if(insert) {                                        //we want to fetch something

      int size = get(insert);

      if(size >= 0) {                                   //success
        s_curr += size;
        {
          QMutexLocker locker(&(this->heap_lock));
          if(final)
            insert->count.ref();                       //now lock is 0 and can be locked

          this->heap.push(insert);
        }
        this->check_queue.open();                      //we should signal the parent cache that we have a new item
        return true;

      } else {                                         //failed transfer put it back, we will keep trying to transfer it...
        QMutexLocker input_locker(&(input->heap_lock));
        input->heap.push(insert);
        return false;
      }
    }
    return false;
  }
};


template<typename Token>
class Transfer: public QThread {
 public:
  Transfer(Cache<Token> *_cache): cache(_cache) {}
 private:
  Cache<Token> *cache;

  void run() {
    cache->loop();
    //end();
  }
};


#endif // GCACHE_H
