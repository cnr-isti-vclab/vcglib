#ifndef TRACKMODE_H
#define TRACKMODE_H

#include <vcg/space/point3.h>
#include <vcg/math/similarity.h>

namespace vcg {

class TrackMode {
public:
  virtual ~TrackMode() {}
  virtual void Draw() {}
  virtual Similarityf Apply(const Point3f &p, const Similarityf &a) { return Similarityf().SetIdentity(); }
};

class SphereMode: public TrackMode {
public:  
  Similarityf Apply(const Point3f &p, const Similarityf &a);
};

class GravityMode: public TrackMode {
public:
};

class CylinderMode: public TrackMode {
public:
  CylinderMode(const Point3f _axis): axis(_axis) { axis.Normalize(); }
protected:
  Point3f axis;
};

class PlaneMode: public TrackMode {
public:
  PlaneMode(const Point3f _x, const Point3f _y): x(_x), y(_y) { x.Normalize(); y.Normalize(); }
  Similarityf Apply(const Point3f &p, const Similarityf &a);
protected:
  Point3f x;
  Point3f y;
};

class LineMode: public TrackMode {
public:
  LineMode(const Point3f _axis): axis(_axis) { axis.Normalize();}
protected:
  Point3f axis;
};

class ScaleMode: public TrackMode {
public:
};

}//namespace 

#endif