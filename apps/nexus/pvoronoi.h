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
Revision 1.5  2004/08/27 00:39:28  ponchio
Rewrote.

Revision 1.4  2004/07/20 14:17:51  ponchio
*** empty log message ***

Revision 1.3  2004/07/01 21:34:59  ponchio
int -> Key

Revision 1.2  2004/06/25 16:47:13  ponchio
Various debug

Revision 1.1  2004/06/24 14:32:45  ponchio
Moved from wrap/nexus

Revision 1.2  2004/06/24 14:19:20  ponchio
Debugged

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

#include "crude.h"

//TODO provide a Sort function, to sort spatially the seeds.

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

    VoronoiPartition() {}  

    void Init(vcg::Box3f &bb) { bbox=bb; ug.SetBBox(bb); }
    unsigned int Add(const vcg::Point3f &p, float weight = 1);
    float Closest(const vcg::Point3f &p, 
		  unsigned int &target, float radius = 0);
    
    class iterator {
    public:
      void operator++();
      const unsigned int operator*();
      bool operator==(const iterator &key);
      bool operator!=(const iterator &key);
    private:
      unsigned int seed;
      friend class VoronoiPartition;
    };
    iterator begin();
    iterator end();
    int size();
    unsigned int count(unsigned int key);
    Seed &operator[](unsigned int key);
    void clear();
    void reload() { ug_seeds = all_seeds; ug.Set(ug_seeds); }
    unsigned int Locate(const vcg::Point3f &p);
    float Priority(const vcg::Point3f &p, unsigned int key);
    
    bool Save(const std::string &file);
    bool Load(const std::string &file);
    unsigned int Save(FILE *fp);
    unsigned int Load(FILE *fp);

    /** Pass iterators to Point3f container and size 
	to estimate optimal radius.
	At the moment strategy is to campion randomly the file.
    */

    static float OptimalRadius(Crude &crude, unsigned int target);
    template <class T> 
      static std::vector<float> OptimalRadii(unsigned int total, 
				 T begin, T end, 
				 vcg::Box3f &box, 
				 std::vector<unsigned int> target) {
      
      //number of samples
      unsigned int n_points = 20;
      std::vector<vcg::Point3f> samples;
      
      T i;
      unsigned int h;
      for(i = begin, h =0; i != end; ++i, h++) 
	if(!((h+1)%(total/n_points)))                  
	  samples.push_back(*i);
      
      
      float step = box.Diag()/10000;
      
      //for every sample i need to record function distance -> number of points
      std::vector< std::vector<int> > scale;
      scale.resize(samples.size());
      for(unsigned int i = 0; i < samples.size(); i++)
	scale[i].resize(10001, 0);
      
      
      //for every point we check distance from samples
      for(i = begin; i != end; ++i) {
	vcg::Point3f &vp = *i;
	for(unsigned int k = 0; k < samples.size(); k++) {
	  float dist = (vp - samples[k]).Norm();  
	  unsigned int pos = (int)(dist/step);
	  if(pos < 10000)
	    scale[k][pos]++;
	}
      }
      
      
      float count =0;
      unsigned int tcount = 0;
      std::vector<int> counting;
      for(int  j = 0; j < 10000; j++) {
	for(unsigned int k = 0; k < samples.size(); k++) 
	  count += scale[k][j];    
	if(count > samples.size() * target[tcount]) {
	  counting.push_back(j);
	  tcount ++;
	  if(tcount >= target.size())
	    j = 10000;
	}
      } 
      std::vector<float> radius;
      for(unsigned int i = 0; i < counting.size(); i++) 
	radius.push_back(2 * step * (counting[i]));
      return radius;
    }
  
    // private:    
  vcg::Box3f bbox;
    vcg::GridStaticPtr< std::vector<Seed> > ug;
    std::vector<Seed> all_seeds;
    std::vector<Seed> ug_seeds;
    std::vector<Seed> seedBuf;
  };


}
#endif

