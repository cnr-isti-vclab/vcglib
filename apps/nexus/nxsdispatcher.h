#ifndef NXS_DISPATCHER_H
#define NXS_DISPATCHER_H

#include <ptypes/pinet.h>
#include <ptypes/pasync.h>
#include <vector>
#include <map>
#include <string>

#include "decimate.h"

namespace nxs {

#define MSG_SEND     MSG_USER + 1
#define MSG_RECEIVE  MSG_USER + 2
#define MSG_FAIL  MSG_USER + 3
  
  class Fragment;
  class Nexus;
  class VChain;
  
  class Server;
  class FragIO;
  class Dispatcher;
  

  class Opener: public pt::thread {
  public:
    Opener(Server *s): thread(false), server(s) {}
    ~Opener() { waitfor(); }
    void execute();
    void cleanup() {}
    
    Server *server;
  };
  
  
  class Server: public pt::ipstream {
  public:
    Server(pt::string host, int port): ipstream(host, port), queue(0), 
      connected(false), opener(this) {}
    
    int queue;
    pt::mutex reading;
    pt::mutex writing;
    bool connected;
    Opener opener;
  };
  

  class Dispatcher: public pt::msgqueue {
  public:
    Dispatcher(Nexus *nx, VChain *ch):
      count(0), maxqueue(3), nexus(nx), chain(ch) {}
    ~Dispatcher();
    
    bool Init(const std::string &file);
    void SendFragment(Fragment *frag);
    void ReceiveFragment(Fragment *in, Fragment *out);
    Server *BestServer();

    void msghandler(pt::message &msg);  
    
    int count;
    int maxqueue;
    Nexus *nexus;
    VChain *chain;
    Decimation mode;
    float scaling;
    std::vector<Server *> servers;
    std::map<int, FragIO *> frags;
  };
  
  class FragIO: public pt::thread {
  public:
    FragIO(Server *se, Dispatcher *di, Fragment *frag): 
      thread(false), server(se), dispatcher(di), fragin(frag) {}
    ~FragIO() { waitfor(); }
    void execute();
    void cleanup() {}
    
    Server *server;
    Dispatcher *dispatcher;
    Fragment *fragin;
  };
  


}

#endif
