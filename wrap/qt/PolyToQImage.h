#ifndef POLYTOQIMAGE_H
#define POLYTOQIMAGE_H
#include <vcg/space/point2.h>
#include <vcg/space/color4.h>
#include <vcg/space/box2.h>
#include <vcg/math/similarity2.h>

class PolyDumperParam
{
public:
  vcg::Color4b backgroundColor;
  int widthPx;
  int heightPx;
  int widthMm;
  int heightMm;
  int dpi;
  bool useDPI;

  PolyDumperParam()
  {
    backgroundColor = vcg::Color4b::Gray;
    widthPx=1024;
    heightPx=1024;
    dpi=72;
    widthMm = 100;
    heightMm = 100;
    useDPI=false;
  }
};

  void dumpPolySet(const char * imageName, std::vector< std::vector< std::vector<vcg::Point2f> > > &polyVecVec, std::vector<vcg::Similarity2f> &trVec, std::vector<std::string> &labelVec, PolyDumperParam &pp);
  void dumpPolySet(const char * imageName, std::vector< std::vector<vcg::Point2f> > &polyVec, std::vector<vcg::Similarity2f> &trVec, PolyDumperParam &pp);
  void dumpPolySet(const char * imageName, std::vector< std::vector<vcg::Point2f> > &polyVec, PolyDumperParam &pp);
  void rectSetToPolySet(std::vector< vcg::Box2f > &rectVec, std::vector< std::vector<vcg::Point2f> > &polyVec);

#endif // POLYTOQIMAGE_H
