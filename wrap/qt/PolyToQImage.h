#ifndef POLYTOQIMAGE_H
#define POLYTOQIMAGE_H
#include <vcg/space/point2.h>
#include <vcg/math/similarity2.h>

int dumpPolySet(const char * imageName, std::vector< std::vector<vcg::Point2f> > &polyVec, std::vector<vcg::Similarity2f> &trVec,int width=1024,int height=1024);

#endif // POLYTOQIMAGE_H
