#include  <ptypes/pinet.h>

class NxsClientPool {
  vector<NxsClient> clients;
  void Push(NxsRequest *request);
};

class NxsClient: public pt::ipstream {
public:
  queue<NxsRequest *> requests;

  void Push(NxsRequest *request);
};


ipaddress  addr = phostbyname("nb-ponchio.isti.cnr.it");
  ipstream client;
  client.set_ip(addr);
  client.set_port(testport);
  try
    {
      client.open();
      
      pout.put("Sending a request to the server...\n");
      client.write("Hello", 6);
      client.flush();
      
      // receive the response
      string rsp = client.line(maxtoken);
      pout.putf("Received: %s\n", pconst(rsp));
      
      // need to close the socket explicitly to gracefully shutdown 
      // the peer host too. otherwise, ~ipstream() will call cancel()
      // and leave the peer in a waiting state (not forever though).
      client.close();
    }
  catch(estream* e)
    {
      perr.putf("Error: %s\n", pconst(e->get_message()));
      delete e;
    }
  
