#include "nxsdispatcher.h"
#include "fragment.h"
#include "decimate.h"
#include <iostream>
#include <ptypes/ptypes.h>

using namespace std;
using namespace vcg;
using namespace nxs;
using namespace pt;

void SaveFragment(Nexus &nexus, VoronoiChain &chain,
		  Fragment &fragin,
		  Fragment &fragout);


void Opener::execute() {
  server->reading.lock();      
  server->writing.lock();
  while(1) {
    if(get_signaled())
      return;
    try {
      server->open();
      server->connected = true;
      server->queue = 0;
      break;
    } catch(...) {
    }
    sleep(4);
  }
  server->reading.unlock();
  server->writing.unlock();
}

void FragIO::execute() {  
  pincrement(&(server->queue));
  server->writing.lock();
  //  cerr << "Writing frag...: " << fragin->id << "\n";

  outmemory outm;
  outm.open();
  fragin->Write(&outm);
  pt::string a = outm.get_strdata();
  try {
    server->write((const char *)a, length(a));
    server->flush();
  } catch (estream *e) {
    perr.putf("Error reading: %s\n", pconst(e->get_message()));
    delete e;
    server->close();
    server->connected = false;
    server->writing.unlock();
    message *msg = new message(MSG_FAIL, (int)fragin);
    dispatcher->post(msg);
    //TODO restart Server!
    return;
  }

  server->reading.lock();
  server->writing.unlock();

  Fragment *out = new Fragment;
  if(!server->waitfor(10000) || (!out->Read(server))) {
    perr.putf("Error reading!!\n");
    server->close();
    server->connected = false;
    server->reading.unlock();
    message *msg = new message(MSG_FAIL, (int)fragin);
    dispatcher->post(msg);
    return;
  }
  server->reading.unlock();
  pdecrement(&(server->queue));
  //  cerr << "Received frag: " << out->id << endl;

  message *msg = new message(MSG_RECEIVE, (int)fragin);
  msg->result = (int)out;
  dispatcher->post(msg);
}

bool Dispatcher::Init(const std::string &file) {
  FILE *fp = fopen(file.c_str(), "rb");
  if(!fp) return false;
  char host[256];
  int port;
  while(fscanf(fp, "%s %d\n", host, &port) == 2) {
    cerr << "Host: " << host << " port: " << port << endl;
    Server *server = new Server(host, port);
    server->opener.start();
    servers.push_back(server);
  }
  fclose(fp);
  if(servers.size() == 0) {
    cerr << "Empty server file!\n";
    return false;
  }
  return true;
}

Dispatcher::~Dispatcher() {
  for(unsigned int i = 0; i < servers.size(); i++) {
    Server *server = servers[i];
    server->opener.signal();
    server->close();
    delete server;
  }
}

void Dispatcher::SendFragment(Fragment *frag) {
  //WARNING this handles no more than 1<<31 fragments!
  frag->id = count++;
  message *msg = new message(MSG_SEND, (int)frag);
  post(msg);
}

Server *Dispatcher::BestServer() {
  Server *best = NULL;
  for(unsigned int i = 0; i < servers.size(); i++){
    if(servers[i]->connected) {
      if((servers[i]->queue <= maxqueue) && 
	 (!best || servers[i]->queue < best->queue)) {

	best = servers[i];
	//	cerr << "best: " << i << " queue: " << best->queue << endl;
      }
    }
  }
  return best;
}

void Dispatcher::ReceiveFragment(Fragment *in, Fragment *out) {
  //lock nexus if run in thread.
  //  cerr << "Saving: " << in->id << endl;
  SaveFragment(*nexus, *chain, *in, *out);

  if(frags.count(in->id)) {
    FragIO *frag = frags[in->id];
    delete frag;
    frags.erase(frags.find(in->id));
  }
  delete in;
  delete out;
}

void Dispatcher::msghandler(message &msg) {
  switch(msg.id) {
  case MSG_FAIL: 
  case MSG_SEND: {
    //get server!
    Server *best = BestServer();
    Fragment *fragin = (Fragment *)(msg.param);

    if(!best) { //no server process locally....
      //      cerr << "Local: " << fragin->id << endl;
      vector<Point3f> newvert;
      vector<unsigned int> newface;
      vector<BigLink> newbord;
      Join(*fragin, newvert, newface, newbord);
      
      float error = Decimate(QUADRIC,
			     (unsigned int)((newface.size()/3) * 0.5),
			     newvert, newface, newbord);
      
      Fragment *fragout = new Fragment;
      
      fragout->error = error;
      fragout->id = fragin->id;
      fragout->seeds = fragin->seeds;
      fragout->seeds_id = fragin->seeds_id;
      Split(*fragout, newvert, newface, newbord); 
      ReceiveFragment(fragin, fragout);
    } else {
      //      cerr << "Server: " << fragin->id << endl;
      FragIO *frag = new FragIO(best, this, fragin);
      if(msg.id == MSG_SEND)
	assert(!frags.count(fragin->id));
      frags[fragin->id] = frag;
      frag->start();
    }
   } break;
  case MSG_RECEIVE:
    ReceiveFragment((Fragment *)(msg.param), (Fragment *)(msg.result));
    break;
  default:
    defhandler(msg);
  }
}
