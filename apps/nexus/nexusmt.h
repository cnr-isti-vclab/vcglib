#ifndef NXS_NEXUS_MT_H
#define NXS_NEXUS_MT_H

#include <vector>
#include <queue>
#include <wrap/gui/frustum.h>

#include "nexusbase.h"
#include "queuepserver.h"
#include "borderserver.h"
#include "prefetch.h"

namespace nxs {

 typedef std::vector<unsigned int> Frag; 
  
 struct Node {  
   std::vector<Node *> in;
   std::vector<Node *> out;
   std::vector<Frag> frags;    
   float error;
   bool visited;
   //   bool pushed;
 };

 struct TNode {
   float error;
   Node *node;
   TNode(Node *n, float e): node(n), error(e) {}
   bool operator<(const TNode &n) { return error < n.error; }
 };

 class Metric {
 public:
   std::vector<PatchInfo> *index;

   float GetError(Node *node);   
   float GetError(Frag &frag);
   virtual float GetError(unsigned int cell) = 0;
   virtual void GetView() {}   
 };

 class FlatMetric: public Metric {
 public:
   void GetView() {}
   float GetError(unsigned int cell);   
 };

 class DeltaMetric: public Metric {
 public:
   vcg::Point3f delta;
   void GetView() {}
   float GetError(unsigned int cell);   
 };

 class FrustumMetric: public Metric {
 public:
   vcg::Frustumf frustum;      

   void GetView() { frustum.GetView(); }
   float GetError(unsigned int cell);   
 };

 /* class Policy {
 public:
   float error;
   int ram_used;
   int ram_size;  

   vector<PatchEntry> *entries;

   void Init();
   bool Expand(TNode &node);
   void NodeVisited(Node *node); 
   };*/

class NexusMt: public NexusBase {
 public:
  //Vertex buffer object mode
  enum Vbo { VBO_AUTO,    //autodetect best size 
	     VBO_OFF,     //no vertex buffer object
	     VBO_FIXED }; //user supplied size
 
  enum MetricKind { FRUSTUM,    //screen error extraction                
		    GEOMETRY,    //geometry error extraction
                    DELTA };    //delta error
  
  enum Mode { POINTS, 
	      SMOOTH,
	      XRAY,
	      HIDDEN_LINE,
	      FLAT_WIRE,
	      FLAT,
              PATCHES };

  enum Component { COLOR   = 0x1, 
		   NORMAL  = 0x2, 
		   TEXTURE = 0x4, 
		   DATA    = 0x8};

  Vbo vbo_mode;

  Metric *metric;
  float target_error;
  
  int extraction_max;
  int extraction_used;


  Mode mode;

  unsigned int components;
  bool use_normals;
  bool use_colors;
  bool use_textures;
  bool use_data;

  //statistics:
  unsigned int tri_rendered;
  unsigned int tri_total;

  std::vector<PServer::Item> visited;

  QueuePServer patches;
  BorderServer borders; 

  Prefetch prefetch;

  NexusMt();
  ~NexusMt();
  
  bool Load(const std::string &filename);
  void Close();

  bool InitGL(Vbo mode = VBO_AUTO, unsigned int vbo_size = 64000000);
  
  void Render();
  
  void SetMetric(MetricKind kind);
  void SetError(float error);
  void SetExtractionSize(unsigned int ram_size);
  void SetPrefetchSize(unsigned int size);
  void SetVboSize(unsigned int vbo_size);


  bool SetMode(Mode mode);
  bool SetComponent(Component c, bool on);
  bool SetComponents(unsigned int mask);
  
  void Draw(std::vector<unsigned int> &selected);
  void Draw(unsigned int cell, QueuePServer::Data &data);
  void Extract(std::vector<unsigned int> &selected);

 protected:
  std::vector<Node> nodes;

  bool Expand(TNode &node);
  void NodeVisited(Node *node); 

  void LoadHistory();
  void ClearHistory();
  void VisitNode(Node *node, std::vector<TNode> &heap);
  void Select(std::vector<unsigned int> &selected);
  Patch &LoadPatch(unsigned int p);
};

}
 
#endif
