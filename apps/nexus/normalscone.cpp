#include "normalscone.h"
#include <iostream>

using namespace std;
using namespace vcg;
using namespace nxs;


ANCone3f::ANCone3f() {
  scaledNormal = Point3f(0,0,0);
  frontAnchor = Point3f(0,0,0);
  backAnchor = Point3f(0,0,0);
}

bool ANCone3f::Frontface(const Point3f &viewPoint) {
  Point3f d = frontAnchor - viewPoint; //Vector from viewPoint to frontAnchor.
  float f = - d * scaledNormal;
  if (f < 0.001 || f * f < d * d)
    return false;
  return true;
}

bool ANCone3f::Backface(const Point3f &viewPoint) {
  Point3f d = backAnchor - viewPoint; //Vector from viewPoint to frontAnchor.
  float f = d * scaledNormal;
  if (f < 0.001 || f * f < d * d)
    return false;
  return true;
}

void ANCone3f::AddNormals(vector<Point3f> &normal, vector<float> &area, float threshold) {
  assert(normal.size() == area.size());
  scaledNormal = Point3f(0,0,0);
  int count = 0;
  vector<Point3f>::iterator i;
  for(i = normal.begin(); i != normal.end(); i++) {
    scaledNormal += *i;
    count++;
  }
  scaledNormal /= count;
  scaledNormal.Normalize();

  double distr[50];
  for(int k = 0; k < 50; k++)
    distr[k] = 0;
  double tot_area = 0;

  vector<float>::iterator j;
  for(i = normal.begin(), j = area.begin(); i != normal.end(); i++, j++) {
    int pos = (int)(49.0 * Angle(scaledNormal, *i)/M_PI);
    if(pos < 0) continue;
    assert(pos >=0 && pos < 50);
    distr[pos] += *j;
    tot_area += *j;
  }

  float tot = 0;
  int best;
  for(best = 0; best < 50; best++) {
    tot += distr[best];
    if(tot > threshold * tot_area)
      break;
  }
  double alpha = M_PI * (best + 1) / 50;
  if(alpha > M_PI/ 2 - 0.1) 
    scaledNormal = Point3f(0,0,0);
  else 
    scaledNormal /= cos(M_PI/2 - alpha);
}

void ANCone3f::AddNormals(vector<Point3f> &normal, float threshold) {
  //assert(normal.size() > 0);
  scaledNormal = Point3f(0,0,0);
  int count = 0;
  vector<Point3f>::iterator i;
  for(i = normal.begin(); i != normal.end(); i++) {
    Point3f norm = *i;
    if(norm.Norm() < 0.00001) continue;
    norm.Normalize();
    scaledNormal += norm;
    count++;
  }
  scaledNormal /= count;
  scaledNormal.Normalize();

  int distr[50];
  for(int k = 0; k < 50; k++)
    distr[k] =0;

  for(i = normal.begin(); i != normal.end(); i++) {
    int pos = (int)(50.0 * Angle(scaledNormal, *i)/M_PI);
    distr[pos]++;
  }
  int tot = 0;
  int best;
  //  cerr << "Distr: ";
  for(best = 0; best < 50; best++) {
    //    cerr << distr[best] << " ";
    tot += distr[best];
    if(tot >= threshold * normal.size())
      break;
  }
  double alpha = M_PI * (best +1) / 50;
  //  cerr << "best: " << best << " alpha: " << alpha << endl;
  if(alpha > M_PI/ 2) {
    scaledNormal = Point3f(0,0,0);
  } else {
    scaledNormal /= cos(M_PI/2 - alpha);
  }
}

void ANCone3f::AddAnchors(vector<Point3f> &anchors) {
  assert(anchors.size() > 0);
  frontAnchor = anchors[0];
  backAnchor = anchors[0];

  float fa = frontAnchor * scaledNormal;
  float fb = -backAnchor * scaledNormal;

  vector<Point3f>::iterator i;
  for(i = anchors.begin(); i != anchors.end(); i++) {
    Point3f &anchor = *i;
    float na = anchor * scaledNormal;
    if(na < fa) {
      frontAnchor = anchor;
      fa = na;
    }
    if(-na < fb) {
      backAnchor = anchor;
      fb = -na;
    }
  }
}

void NCone3s::Import(const ANCone3f &c) {
  Point3f normal = c.scaledNormal;
  float len = normal.Norm();
  if(len != 0) 
    normal /= len;

  assert(normal[0] <= 1 && normal[1] <= 1 && normal[2] <= 1);
  assert(normal[0] >= -1 && normal[1] >= -1 && normal[2] >= -1);
  n[0] = (short)(normal[0] * 32766);
  n[1] = (short)(normal[1] * 32766);
  n[2] = (short)(normal[2] * 32766);
  //i want to rapresent number from 1 to 10
  if(len > 10.0f) len = 10.0f;
  if(len < -10.0f) len = -10.0f;
  n[3] = (short)(len * 3276);
}




bool NCone3s::Backface(const vcg::Sphere3f &sphere, 
		       const vcg::Point3f &view) const {
  vcg::Point3f norm(n[0]/32766.0f, n[1]/32766.0f, n[2]/32766.0f);
  vcg::Point3f d = (sphere.Center() + norm * sphere.Radius()) - view;
  norm *= n[3]/3276.0f;

  float f = d * norm;
  if (f < 0.001 || f * f < d * d) 
    return false;
  return true;
}

bool NCone3s::Frontface(const vcg::Sphere3f &sphere, 
			const vcg::Point3f &view) const {
  vcg::Point3f norm(n[0]/32766.0f, n[1]/32766.0f, n[2]/32766.0f);
  vcg::Point3f d = (sphere.Center() - norm * sphere.Radius()) - view;
  norm *= n[3]/3276.0f;
  
  float f = -d * norm;
  if (f < 0.001 || f * f < d * d)
    return false;
  return true;
}
