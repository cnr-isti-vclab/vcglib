/****************************************************************************
* VCGLib                                                            o o     *
* Visual and Computer Graphics Library                            o     o   *
*                                                                _   O  _   *
* Copyright(C) 2004                                                \/)\/    *
* Visual Computing Lab                                            /\/|      *
* ISTI - Italian National Research Council                           |      *
*                                                                    \      *
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
/****************************************************************************
  History

$Log: not supported by cvs2svn $

****************************************************************************/

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
  
