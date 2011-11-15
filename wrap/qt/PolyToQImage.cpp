#include <QImage>
#include <QPainter>

#include "PolyToQImage.h"
#include <wrap/qt/col_qt_convert.h>
using namespace vcg;
using namespace std;

int dumpPolySet(const char * imageName, vector< vector<Point2f> > &polyVec, int width, int height)
{
  Box2f bb;
  for(size_t i=0;i<polyVec.size();++i)
    for(int j=0;j<polyVec[i].size();++j)
      bb.Add(polyVec[i][j]);

  Similarity2f sim;
  sim.sca = min(float(width)/bb.DimX(),float(height)/bb.DimY());
  sim.tra = Point2f(width/2.0f,height/2.0f)-bb.Center()*sim.sca;
  vector<Similarity2f> trVec(polyVec.size(),sim);

  return dumpPolySet(imageName,polyVec,trVec,width,height);
}

int dumpPolySet(const char * imageName, vector< vector< vector<Point2f> > > &polyVecVec, vector<Similarity2f> &trVec, int width, int height)
{
  assert(polyVecVec.size() == trVec.size());

  QImage img(width,height,QImage::Format_RGB32);
  img.fill(qRgb(128,128,128));
  QPainter painter(&img);           // paint in picture
  for(size_t i=0;i<polyVecVec.size();++i)
  {
    QVector<QPointF> ppQ;
    for(int j=0;j<polyVecVec[i][0].size();++j)
    {
      Point2f pp=trVec[i]*polyVecVec[i][0][j];
      ppQ.push_back(QPointF(pp[0],pp[1]));
    }

    QBrush bb(vcg::ColorConverter::ToQColor(Color4b::Scatter(polyVecVec.size(),i)));
    painter.setBrush(bb);
    painter.drawPolygon(&*ppQ.begin(),ppQ.size(),Qt::OddEvenFill);
  }
  painter.end();
  img = img.mirrored(false,true);
  img.save(imageName);

  int emptyCnt=0;
  for(int i=0;i<width;++i)
    for(int j=0;j<height;++j)
      if(img.pixel(i,j) == qRgb(128,128,128)) emptyCnt++;

  return emptyCnt;
}

int dumpPolySet(const char * imageName, vector< vector<Point2f> > &polyVec, vector<Similarity2f> &trVec, int width, int height)
{
  vector< vector< vector<Point2f> > > polyVecVec(polyVec.size());
  for(size_t i=0;i<polyVec.size();++i)
  {
    polyVecVec[i].resize(1);
    polyVecVec[i][0]=polyVec[i];
  }
  dumpPolySet(imageName,polyVecVec,trVec,width,height);
}

