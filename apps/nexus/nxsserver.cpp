#include <ptypes/ptime.h>
#include <ptypes/pinet.h>
#include <ptypes/pasync.h>


#include "fragment.h"
#include "decimate.h"

using namespace pt;
using namespace nxs;
using namespace vcg;
using namespace std;


class FragOutQueue: public msgqueue {
public:
  FragOutQueue(ipstream &cli): client(cli) {}
  void msghandler(message &msg) {
    if(msg.id != MSG_USER) {
      defhandler(msg);
      return;
    }
    Fragment &fragment = *(Fragment *)(msg.param);

    //    pout.putf("Sending: %d\n", fragment.id);
    outmemory outm;
    outm.open();
    fragment.Write(&outm);
    pt::string a = outm.get_strdata();
    try {
      client.write((const char *)a, length(a));
      client.flush();
      pout.putf("Sent fragment id: %d\n", fragment.id);
    } catch (estream *e) {
      perr.putf("Error: %s\n", pconst(e->get_message()));
      delete e;
      posturgent(MSG_QUIT);
    }
    delete (Fragment *)(msg.param);
  }
  ipstream &client;
};

class FragInQueue: public msgqueue {
public:
  FragInQueue(FragOutQueue &o): out(o) {}
  void msghandler(message &msg) {
    if(msg.id != MSG_USER) {
      if(msg.id == MSG_QUIT)
	out.posturgent(MSG_QUIT);
      defhandler(msg);
      return;
    }
    Fragment &fragin = *(Fragment *)(msg.param);
    //    pout.putf("Processing: %d\n", fragin.id);
    vector<Point3f> newvert;
    vector<unsigned int> newface;
    vector<BigLink> newbord;
    Join(fragin, newvert, newface, newbord);
    
    float error = Decimate(QUADRIC,
    			   (unsigned int)((newface.size()/3) * 0.5),
    			   newvert, newface, newbord);
    
    message *outmsg = new message(MSG_USER);
    outmsg->param = (int)(new Fragment);
    Fragment &fragout = *(Fragment *)(outmsg->param);

    fragout.error = error;
    fragout.id = fragin.id;
    fragout.seeds = fragin.seeds;
    fragout.seeds_id = fragin.seeds_id;
    Split(fragout, newvert, newface, newbord);
    out.post(outmsg);
    delete (Fragment *)(msg.param);
  }

  FragOutQueue &out;
};


class Reader: public thread {  
public:
  Reader(ipstream &cli, FragInQueue &que): 
    thread(false), client(cli), queue(que) {}

  ~Reader() {waitfor(); }

  void execute() {
    while(1) {
      if(get_signaled()) return;
      message *msg = new message(MSG_USER);
      msg->param = (int)(new Fragment);
      Fragment &fragment = *(Fragment *)(msg->param);
      if(!fragment.Read(&client)) {
	pout.putf("Could not read!\n");
	queue.posturgent(MSG_QUIT);
	return;
      }
      queue.post(msg);
      //      pout.putf("Incoming: %d\n", fragment.id);
    }
  }
  void cleanup() {}

  ipstream &client;
  FragInQueue &queue;
};

class Worker: public thread {
public:
  Worker(FragInQueue &que): thread(false), queue(que) {}
  ~Worker() { waitfor(); }

  void execute() {
      queue.run();
  }
  void cleanup() {}
  FragInQueue &queue;
};

class Writer: public thread {  
public:
  Writer(FragOutQueue &que): thread(false), queue(que) {}
  ~Writer() {waitfor(); }
  void execute() {
    queue.run();
  }
  void cleanup() {}
  FragOutQueue &queue;
};


void servermain(ipstmserver& svr) {
  ipstream client;
  
  while(true) {
    // serve() will wait for a connection request and will prepare
    // the supplied ipstream object for talking to the peer.
    svr.serve(client);
    perr.putf("Serving clients!\n");
    if (client.get_active()) {
      try {
	pout.putf("Incoming connection\n");
	FragOutQueue out(client);
	FragInQueue in(out);
	Reader reader(client, in);
	Worker worker(in);
	Writer writer(out);

	reader.start();
	worker.start();
	writer.start();

	reader.waitfor();
	worker.waitfor();
	writer.waitfor();

	client.flush();	
	client.close();
      } catch(estream* e)	{
	perr.putf("Error: %s\n", pconst(e->get_message()));
	delete e;
      }
    }
    perr.putf("Restarting\n");
  }    
}

int main(int argc, char *argv[]) {
  ipstmserver svr;

  int port = 10102;
  if(argc == 2) {
    port = atoi(argv[1]);
    if(port < 1024) {
      perr.putf("Error: invalid port: %s\n", argv[1]);
      return -1;
    }
  }

  try {
    // bind to all local addresses on port 8085
    svr.bindall(port);

    pout.putf("Ready to answer queries on port %d\n", port);
    
    // enter an infinite loop of serving requests
    servermain(svr);
  } catch(estream* e) {
    perr.putf("FATAL: %s\n", pconst(e->get_message()));
    delete e;
  }
  return 0;
}
