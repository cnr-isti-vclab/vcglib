#ifndef NXS_NEXUS_MT_H
#define NXS_NEXUS_MT_H


#include "nexus.h"
#include <vector>
#include <queue>

#include <wrap/gui/frustum.h>


namespace nxs {

 class Frag:public std::vector<unsigned int> {};
  
 struct Node {  
   std::vector<Node *> in;
   std::vector<Node *> out;
   std::vector<Frag> frags;    
   float error;
   bool visited;        
 };
 
 class Policy {
 public:
   virtual bool Expand(unsigned int patch, Nexus::Entry &entry) = 0;
   virtual void Visit(Node *node, std::queue<Node *> &qnode);
 };

 class FrustumPolicy: public Policy {
 public:
   vcg::Frustumf frustum;
   float error;
   void GetView() { frustum.GetView(); }
   bool Expand(unsigned int patch, Nexus::Entry &entry);
 };

class NexusMt: public Nexus {
 private:
  std::vector<Node> nodes;

 public:
  void LoadHistory();
  void ClearHistory();

  void ExtractFixed(std::vector<unsigned int> &selected, float error);
  void Extract(std::vector<unsigned int> &selected, Policy *policy);

 protected:
  void Select(std::vector<unsigned int> &selected);
};

}
 
#endif
