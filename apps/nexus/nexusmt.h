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
Revision 1.24  2005/02/17 15:39:44  ponchio
Reorderes statistics a bit.

Revision 1.23  2005/02/14 17:11:07  ponchio
aggiunta delle sphere

Revision 1.22  2005/02/10 09:18:20  ponchio
Statistics.

Revision 1.21  2005/02/08 12:43:03  ponchio
Added copyright


****************************************************************************/

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
    enum Attr { COLOR = 0x1, NORMAL = 0x2, TEXTURE = 0x4, DATA = 0x8, SPHERES = 0x10 };

    Mode mode;
    unsigned int attrs;
    
    DrawContest(Mode m = SMOOTH, unsigned int a = 0xf): mode(m), attrs(a) {}
    void SetAttr(Attr attr, bool value);
    bool HasComponent(Attr attr);
  };

  struct Stats {
    //per frame data...
    float tri;       //k triangles rendered.	
    float extr;      //k triangles extracted
    float disk_tri;      //k triangles readed from disk

    int log_size;

    deque<float> error;  //max error in extraction (push_front pop_back)
    deque<float> time;
    deque<float> disk;  //kdisk readed per frame

    float fps; //averaged over 8 frames

    Watch watch;

    Stats(): log_size(48), fps(0.0f) {}
    void Start();
    void Disk(float disk);
    void Error(float error);
    void Stop();
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
