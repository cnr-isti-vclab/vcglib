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
   virtual bool Expand(unsigned int patch, Nexus::PatchInfo &entry) = 0;
   virtual void GetView() {}
   virtual void Visit(Node *node, std::queue<Node *> &qnode);
 };

 class FrustumPolicy: public Policy {
 public:
   vcg::Frustumf frustum;
   float error;

   FrustumPolicy(float _err = 4): error(_err) {}
   void GetView() { frustum.GetView(); }
   bool Expand(unsigned int patch, Nexus::PatchInfo &entry);
 };

class NexusMt: public Nexus {
 private:
  std::vector<Node> nodes;

 public:
  //Vertex buffer object mode
  enum Vbo { VBO_AUTO,    //autodetect best size (fallback if not VBO)
	     VBO_AUTO_ON, //autodetect best size must use VBO
	     VBO_OFF,     //no vertex buffer object
	     VBO_FIXED }; //user supplied size
 
  enum PolicyKind { FRUSTUM,    //screen error extraction
		    GEOMETRY }; //geometry error extraction

  enum Mode { POINTS, 
	      SMOOTH,
	      XRAY,
	      HIDDEN_LINE,
	      FLAT_WIRE,
	      FLAT,
              DEBUG };

  enum Component { COLOR   = 0x1, 
		   NORMAL  = 0x2, 
		   TEXTURE = 0x4, 
		   DATA    = 0x8};

  Vbo vbo;
  unsigned int vbo_size;

  Policy *policy;
  float error;
  bool realtime;

  Mode mode;

  unsigned int components;
  bool use_normals;
  bool use_colors;
  bool use_textures;
  bool use_data;
  
  NexusMt();
  ~NexusMt();
  
  bool Load(const std::string &filename);
  bool InitGL();
  
  void Render();
  void SetPolicy(Policy *policy, bool realtime = true);
  void SetPolicy(PolicyKind kind, float error, bool realtime = true);
  void SetVbo(Vbo mode, unsigned int vbo_size = 0, 
	      unsigned int ram_size = 128000000);
  bool SetMode(Mode mode);
  bool SetComponent(Component c, bool on);
  bool SetComponents(unsigned int mask);
  
  void ExtractFixed(std::vector<unsigned int> &selected, float error);
  void Extract(std::vector<unsigned int> &selected, Policy *policy);

 protected:
  void LoadHistory();
  void ClearHistory();
  void Select(std::vector<unsigned int> &selected);
  Patch &LoadPatch(unsigned int p);
};

}
 
#endif
