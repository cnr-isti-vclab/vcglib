#include <assert.h>

//These header are neede byt tristipper...
#include <list>
#include <vector>
#include <map>

#include "tristripper/tri_stripper.h"
using namespace triangle_stripper;

#include "strip.h"

using namespace std;
using namespace nxs;

void nxs::ComputeTriStrip(unsigned short nfaces, unsigned short *faces, 
			  vector<unsigned short> &strip) {
  
  vector<unsigned int> index;
  index.resize(nfaces*3);
  
  for(int i = 0; i < nfaces*3; i++) 
    index[i] = faces[i];

  int cache_size = 16;
  tri_stripper stripper(index);
  stripper.SetCacheSize(cache_size);		
  // = 0 will disable the cache optimizer
  stripper.SetMinStripSize(0);
  tri_stripper::primitives_vector primitives;
  stripper.Strip(&primitives);

  if(primitives.back().m_Indices.size() < 3) 
    primitives.pop_back();

  //TODO do this when mounting strips together.
  if(primitives.back().m_Type == tri_stripper::PT_Triangles) {
    tri_stripper::primitives p;
    p = primitives.back();
    primitives.pop_back();		
    for(unsigned int i = 0; i < p.m_Indices.size(); i += 3) {
      tri_stripper::primitives s;
      s.m_Type = tri_stripper::PT_Triangle_Strip;
      s.m_Indices.push_back(p.m_Indices[i]);
      s.m_Indices.push_back(p.m_Indices[i+1]);
      s.m_Indices.push_back(p.m_Indices[i+2]);
      primitives.push_back(s);
    }
  }
  
  for(unsigned int i = 0; i < primitives.size(); i++) {
    tri_stripper::primitives &primitive = primitives[i];
    assert(primitive.m_Indices.size() != 0);
    int len = primitive.m_Indices.size();
    for(int l = 0; l < len; l++)  		
      strip.push_back(primitive.m_Indices[l]);
    
    
    if(i < primitives.size()-1) { //not the last primitive.
      strip.push_back(primitive.m_Indices[len-1]);
      //TODO optimize this!
      if((len%2) == 1) 	//do not change orientation....
	strip.push_back(primitive.m_Indices[len-1]);
      strip.push_back(primitives[i+1].m_Indices[0]);
    }			
  }
}
