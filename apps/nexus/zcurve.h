#ifndef VCG_ZCURVE_H

#include <vcg/space/box3.h>

namespace vcg {

class ZCurve: public Box3f {
public:

  unsigned int side;

  ZCurve(): side(1024) {}

  void SetSide(int s) { assert(s <= (1<<10)); side = s; }

  unsigned int Pos(const Point3f &p) const {
    assert(!IsNull());
    unsigned int position = 0;
    
    float x = (p[0] - min[0])/(max[0] - min[0]);
    float y = (p[1] - min[1])/(max[1] - min[1]);
    float z = (p[2] - min[2])/(max[2] - min[2]);

    unsigned int s = side;
    while(s > 1) {
      position *= 8;
      x *= 2;
      y *= 2;
      z *= 2;
      int dx = (int)floor(x);
      int dy = (int)floor(y);
      int dz = (int)floor(z);
      position += dx + 2 * dy + 4 * dz;
      x -= dx;
      y -= dy;
      z -= dz;
      s /= 2;
    }
    return position;
  }
};

}


#endif
