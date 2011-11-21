#ifndef POLYTOQIMAGE_H
#define POLYTOQIMAGE_H
#include <vcg/space/point2.h>
#include <vcg/space/box2.h>
#include <vcg/math/similarity2.h>

void dumpPolySet(const char * imageName, std::vector< std::vector< std::vector<vcg::Point2f> > > &polyVecVec, std::vector<vcg::Similarity2f> &trVec, int width, int height);
void dumpPolySet(const char * imageName, std::vector< std::vector<vcg::Point2f> > &polyVec, std::vector<vcg::Similarity2f> &trVec, int width=1024,int height=1024);
void dumpPolySet(const char * imageName, std::vector< std::vector<vcg::Point2f> > &polyVec, int width=1024,int height=1024);
void rectSetToPolySet(std::vector< vcg::Box2f > &rectVec, std::vector< std::vector<vcg::Point2f> > &polyVec);
#endif // POLYTOQIMAGE_H
