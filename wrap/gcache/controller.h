#ifndef GCACHE_CONTROLLER_H
#define GCACHE_CONTROLLER_H

#include <QDebug>
#include "cache.h"

/** Allows to insert tokens, update priorities and generally control the cache.
*/

template <class Token>
class Controller {
 public:
  ///should be private
  std::vector<Token *> tokens;   //tokens waiting to be added
  bool quit;                     //gracefully terminate.
  bool paused;
  bool stopped;

 public:
  ///should be protected
  Provider<Token> provider;
  ///should be protected
  std::vector<Cache<Token> *> caches;

  Controller(): quit(false), paused(false), stopped(true) {}
  ~Controller() { finish(); }

  ///called before the cache is started to add a cache in the chain
  /** The order in which the caches are added is from the lowest to the highest. */      
  void addCache(Cache<Token> *cache) {
    if(caches.size() == 0)
      cache->setInputCache(&provider);
    else
      cache->setInputCache(caches.back());
    assert(cache->input);
    caches.push_back(cache);
  }
  ///insert the token in the cache if not already present (actual insertion is done on updatePriorities)
  bool addToken(Token *token) {
    if(token->count.testAndSetOrdered(Token::OUTSIDE, Token::CACHE)) {
      tokens.push_back(token);
      return true;
    }
    return false;
  }

  ///WARNING: migh stall for the time needed to drop tokens from cache.
  //FUNCTOR has bool operator(Token *) and return true to remove
  template<class FUNCTOR> void removeTokens(FUNCTOR functor) {
    pause(); //this might actually be unnecessary if you mark tokens to be removed
    for(int i = (int)caches.size()-1; i >= 0; i--)
      caches[i]->flush(functor);
    provider.flush(functor);

    resume();
  }

  ///if more tokens than m present in the provider, lowest priority ones will be removed
  void setMaxTokens(int m) {
    QMutexLocker l(&provider.heap_lock);
    provider.max_tokens = m;
  }

  ///ensure that added tokens are processed and existing ones have their priority updated.
  void updatePriorities() {

    if(tokens.size()) {
      QMutexLocker l(&provider.heap_lock);
      for(unsigned int i = 0; i < tokens.size(); i++)
        provider.heap.push(tokens[i]);
      tokens.clear();
    }

    provider.pushPriorities();
    for(unsigned int i = 0; i < caches.size(); i++)
      caches[i]->pushPriorities();
  }

  ///start the various cache threads.
  void start() {
    if(!stopped) return;
    assert(!paused);
    assert(caches.size() > 1);
    caches.back()->final = true;
    for(unsigned int i = 0; i < caches.size(); i++) //cache 0 is a provider, and his thread is not running.
      caches[i]->start();
    stopped = false;
  }
  ///stops the ache threads
  void stop() {
    if(stopped) return;
    if(paused) resume();
    //stop threads
    for(int i = caches.size()-1; i >= 0; i--) {
      caches[i]->quit = true;                      //hmmmmmmmmmmmmmm not very clean.
      if(i == 0)
        provider.check_queue.open();
      else
        caches[i-1]->check_queue.open();           //cache i listens on queue i-1
      caches[i]->wait();
    }
    stopped = true;
  }

  void finish() {
    stop();
  }

  void pause() {
    if(paused) return;
    provider.check_queue.lock();
    for(unsigned int i = 0; i < caches.size()-1; i++)
      caches[i]->check_queue.lock();
/*    provider.heap_lock.lock();
    for(unsigned int i = 0; i < caches.size(); i++)
      caches[i]->heap_lock.lock(); */
    paused = true;
  }

  void resume() {
    if(!paused) return;
    provider.check_queue.unlock();
    for(unsigned int i = 0; i < caches.size()-1; i++)
      caches[i]->check_queue.unlock();


/*    provider.heap_lock.unlock();
    for(unsigned int i = 0; i < caches.size(); i++)
      caches[i]->heap_lock.unlock(); */
    paused = false;
  }
  ///empty all caches AND REMOVES ALL TOKENS!
  void flush() {
    pause();
    for(int i = (int)caches.size()-1; i >= 0; i--)
      caches[i]->flush();
    provider.heap.clear();
    resume();
  }

  bool isWaiting() {
    for(int i = (int)caches.size() -1; i >= 0; i--) {
      if(!caches[i]->input->check_queue.isWaiting()) return false;
    }
    return true;
  }
};


#endif // CONTROLLER_H
