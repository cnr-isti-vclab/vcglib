#ifndef NXS_HISTORY_H
#define NXS_HISTORY_H

#include <vector>

//TODO fix a bit better the quick <-> updates duality

namespace nxs {
  
  class History {
  public:

    enum Mode { QUICK = 1, UPDATES = 2 };
    
    struct Update {
      std::vector<unsigned int> erased;
      std::vector<unsigned int> created;
    };

    struct Cell {
      unsigned int patch;
      //      float error;
    };

    struct Node;

    struct Link {
      Node *node;

      typedef Cell *iterator;
      iterator begin() { return frag_begin; }
      iterator end() { return frag_begin + frag_size; }
      unsigned int size() { return frag_size; }

      Cell *frag_begin;
      unsigned int frag_size;
    };

    struct Node {
      typedef Link *iterator;
    
      iterator in_begin() { return in_link_begin; }
      iterator in_end() { return in_link_begin + in_link_size; }
      unsigned int size() { return in_link_size; }

      iterator out_begin() { return out_link_begin; }
      iterator out_end() { return out_link_begin + out_link_size; }
      unsigned int out_size() { return out_link_size; }

      Link *in_link_begin;
      unsigned int in_link_size;
      Link *out_link_begin;
      unsigned int out_link_size;
    };

    Node *nodes;
    Link *in_links;
    Link *out_links;
    Cell *frags;

    std::vector<Update> updates;

    History(): nodes(NULL), in_links(NULL), out_links(NULL), frags(NULL), 
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
    bool UpdatesToQuick();
    bool IsQuick() { return buffer != NULL; }

    int &quick() { return ((int *)buffer)[0]; }
    int &n_nodes() { return ((int *)buffer)[1]; }
    int &n_in_links() { return ((int *)buffer)[2]; }
    int &n_out_links() { return ((int *)buffer)[3]; }
    int &n_frags() { return ((int *)buffer)[4]; }

    typedef Node *iterator;
    iterator begin() { return nodes; }
    iterator end() { return nodes + n_nodes(); }
  protected:
    unsigned int size;
    char *buffer;

    bool LoadPointers();
  };
  
}
#endif
