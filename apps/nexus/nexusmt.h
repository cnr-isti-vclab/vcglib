#ifndef NXS_MT_H
#define NXS_MT_H

#include <vector>
#include <queue>
#include <wrap/gui/frustum.h>

#include "nexus.h"
#include "history.h"
#include "extraction.h"
#include "metric.h"
#include "preload.h"
#include "watch.h"

namespace nxs {

  struct DrawContest {

    enum Mode { POINTS, SMOOTH,	XRAY, HIDDEN_LINE, FLAT_WIRE, FLAT, PATCHES };
    enum Attr { COLOR = 0x1, NORMAL = 0x2, TEXTURE = 0x4, DATA = 0x8 };

    Mode mode;
    unsigned int attrs;
    
    DrawContest(Mode m = SMOOTH, unsigned int a = 0xf): mode(m), attrs(a) {}
    void SetAttr(Attr attr, bool value);
    bool HasComponent(Attr attr);
  };

  struct Stats {
    float ktri;       //k triangles rendered.	
    float kdisk;      //k readed per frame (mean)
    float kdisk_peak; //k readed peak.
    float fps;
    float fps_peak;   //low fps peaks

    //double last_time;
    unsigned int count;

    Watch watch;

    Stats(): count(0) {}
    void Init();
  };

  class NexusMt: public Nexus {
  public:
    bool use_vbo;
    bool prefetch;

    unsigned int vbo_used; //TODO remember to zero it!

    Preload preload;

    NexusMt();
    ~NexusMt();
    
    bool Load(const std::string &filename);
    //    void Close();
    
    bool InitGL(bool use_vbo = true);
    
    void Render(DrawContest contest = DrawContest());
    void Render(Extraction &extraction, 
		DrawContest &contest,
		Stats *stats = NULL);

    void SetPreload(bool on);

    void Flush(bool all = true);
    Patch &GetPatch(unsigned int patch, float error, bool flush = true);
    bool CanAdd(Item &item);
  protected:
    std::vector<Item> heap;

    void FlushPatch(unsigned int id);
    void LoadVbo(Entry &entry);
    void FlushVbo(Entry &entry);
    void Draw(unsigned int cell, DrawContest &contest);

  };

  /*  class NexusViewer {
  public:
    NexusMt &mt;

    DrawContest contest;
    Extraction extraction;
    Statistics stats;
    
    NexusViewer(NexyusMt &_mt): mt(&mt) {}
    void Render();


    bool SetMode(DrawContest::Mode mode);
    bool SetComponent(DrawContest::Component c, bool on);
    bool SetComponents(unsigned int mask);

    void SetMetric(MetricKind kind);
    void SetError(float error);

    //Uinits expressed in Kb
    void SetExtractionSize(unsigned int ksize);
    void SetDrawSize(unsigned int ksize);
    void SetDiskDelta(unsigned int kdelta);
    //void SetPrefetchSize(unsigned int size);
    };*/
  



}

#endif
