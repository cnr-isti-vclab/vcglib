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
Revision 1.1  2004/06/23 17:17:46  ponchio
Created


****************************************************************************/

#ifndef NXS_NET_GRID_H
#define NXS_NET_GRID_H

#include <vector>
#include <set>
#include <string>
#include <stdio.h>

#include <vcg/space/point3.h>
#include <vcg/space/box3.h>
#include <vcg/space/index/grid_static_ptr.h>

namespace nxs {

  class Seed {
  public:
    vcg::Point3f p;  
    float weight;
    typedef float ScalarType;
    bool Dist(const vcg::Point3f & point, float &mindist, vcg::Point3f &res);
    void GetBBox(vcg::Box3f &b) {b.Set(p);}
    bool IsD() { return false; }

    Seed(): weight(1) {}
    Seed(const vcg::Point3f &point): p(point), weight(1) {}
    Seed(const vcg::Point3f &point, const float w): 
      p(point), weight(w) {}
    
    inline float Dist(const vcg::Point3f &q) const { 
      return weight * vcg::Distance(p,q); 
    }
    inline float SquaredDist(const vcg::Point3f &q) const  { 
      return weight * weight *vcg::SquaredDistance(p,q); 
    }  
  };
  

  class VoronoiPartition {
  public:
    enum { MAX_BUF=25 };
    typedef int Key;

    VoronoiPartition() {}  

    void Init(vcg::Box3f &bb) { bbox=bb; ug.SetBBox(bb); }
    int Add(vcg::Point3f p, float weight = 1);
    float Closest(const vcg::Point3f &p, Key &target, float radius = 0);
    
    class iterator {
    public:
      void operator++();
      const Key operator*();
      bool operator==(const iterator &key);
      bool operator!=(const iterator &key);
    private:
      int seed;
      friend class VoronoiPartition;
    };
    iterator begin();
    iterator end();
    int size();
    unsigned int count(Key key);
    Seed &operator[](Key key);
    void clear();
    Key Locate(const vcg::Point3f &p);
    float Priority(const vcg::Point3f &p, Key key);
    
    bool Save(const std::string &file);
    bool Load(const std::string &file);
    unsigned int Save(FILE *fp);
    unsigned int Load(FILE *fp);
  private:    
    vcg::Box3f bbox;
    vcg::GridStaticPtr< std::vector<Seed> > ug;
    std::vector<Seed> all_seeds;
    std::vector<Seed> ug_seeds;
    std::vector<Seed> seedBuf;
  };

}
#endif

