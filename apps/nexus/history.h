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
Revision 1.6  2005/02/20 00:43:23  ponchio
Less memory x extraction.  (removed frags)

Revision 1.5  2005/02/19 16:22:45  ponchio
Minor changes (visited and Cell)

Revision 1.4  2005/02/17 16:40:35  ponchio
Optimized BuildLevels.

Revision 1.3  2005/02/08 12:43:03  ponchio
Added copyright


****************************************************************************/

#ifndef NXS_HISTORY_H
#define NXS_HISTORY_H

#include <vector>
#include <map>

//TODO fix a bit better the quick <-> updates duality

namespace nxs {
  
  class Nexus;
  class History {
  public:

    enum Mode { QUICK = 1, UPDATES = 2 };
    
    struct Update {
      std::vector<unsigned int> erased;
      std::vector<unsigned int> created;
    };

    struct Node;

    struct Link {
      Node *node;
      unsigned int begin; //begin patch of the fragment
      unsigned int end;   //end patch of the fragment

      typedef unsigned int iterator;
      unsigned int size() { return end - begin; }
    };

    struct Node {
      typedef Link *iterator;
      Link *in_begin, *in_end;
      Link *out_begin, *out_end;
    };

    Node *nodes;
    Link *in_links;
    Link *out_links;

    std::vector<Update> updates;

    History(): nodes(NULL), in_links(NULL), out_links(NULL),// frags(NULL), 
      buffer(NULL) {}
    ~History();
   
    Node *Root() { return &nodes[0]; }

    void Clear();
    void ClearQuick();
    void ClearUpdates();
    //Owns memory afterwards.. do not free mem.
    bool Load(unsigned int size, char *mem);
    bool LoadQuick(unsigned int size, char *mem);
    bool LoadUpdates(unsigned int size, char *mem);

    //after these call history is invalid! and memory returned must be freed...
    char *Save(unsigned int &size); //autodetect
    char *SaveQuick(unsigned int &size);
    char *SaveUpdates(unsigned int &size);

    bool QuickToUpdates();
    bool UpdatesToQuick(Nexus &nexus);
    bool IsQuick() { return buffer != NULL; }

    void BuildLevels(std::vector<int> &levels);

    int &quick() { return ((int *)buffer)[0]; }
    int &n_nodes() { return ((int *)buffer)[1]; }
    int &n_in_links() { return ((int *)buffer)[2]; }
    int &n_out_links() { return ((int *)buffer)[3]; }

    //    typedef Node *iterator;
    //    iterator begin() { return nodes; }
    //    iterator end() { return nodes + n_nodes(); }
  protected:
    unsigned int size;
    char *buffer;

    bool LoadPointers();
  };
  
}
#endif
