#include <QImage>
#include <QSvgGenerator>
#include <QPainter>

#include "PolyToQImage.h"
#include <wrap/qt/col_qt_convert.h>
using namespace vcg;
using namespace std;

void rectSetToPolySet(vector< Box2f > &rectVec, vector< vector<Point2f> > &polyVec)
{
  polyVec.clear();
  for(size_t i=0;i<rectVec.size();++i)
  {
    Box2f &b=rectVec[i];
      polyVec.resize(polyVec.size()+1);
      polyVec.back().push_back(b.min);
      polyVec.back().push_back(Point2f(b.max[0],b.min[1]));
      polyVec.back().push_back(b.max);
      polyVec.back().push_back(Point2f(b.min[0],b.max[1]));

  }
}

void dumpPolySet(const char * imageName, vector< vector<Point2f> > &polyVec, int width, int height)
{
  Box2f bb;
  for(size_t i=0;i<polyVec.size();++i)
    for(int j=0;j<polyVec[i].size();++j)
      bb.Add(polyVec[i][j]);

  Similarity2f sim;
  sim.sca = min(float(width)/bb.DimX(),float(height)/bb.DimY());
  sim.tra = Point2f(width/2.0f,height/2.0f)-bb.Center()*sim.sca;
  vector<Similarity2f> trVec(polyVec.size(),sim);

  dumpPolySet(imageName,polyVec,trVec,width,height);
}

void dumpPolySet(const char * imageName, vector< vector< vector<Point2f> > > &polyVecVec, vector<Similarity2f> &trVec, int width, int height)
{
  assert(polyVecVec.size() == trVec.size());

  QImage img(width,height,QImage::Format_RGB32);
  img.fill(qRgb(128,128,128));

  QSvgGenerator svg;
  svg.setFileName(imageName);
  QPainter painter;

  if(QString(imageName).endsWith("svg",Qt::CaseInsensitive))
        painter.begin(&svg);
  else  painter.begin(&img);

  for(size_t i=0;i<polyVecVec.size();++i)
  {
    painter.resetTransform();
    painter.translate(trVec[i].tra[0],trVec[i].tra[1]);
    painter.rotate(math::ToDeg(trVec[i].rotRad));
    painter.scale(trVec[i].sca,trVec[i].sca);
    QVector<QPointF> ppQ;
    for(int j=0;j<polyVecVec[i][0].size();++j)
    {
      Point2f pp=polyVecVec[i][0][j];
      ppQ.push_back(QPointF(pp[0],pp[1]));
    }

    QBrush bb(vcg::ColorConverter::ToQColor(Color4b::Scatter(polyVecVec.size(),i)));
    painter.setBrush(bb);
    painter.drawPolygon(&*ppQ.begin(),ppQ.size(),Qt::OddEvenFill);
  }
  painter.end();

 if(!QString(imageName).endsWith("svg",Qt::CaseInsensitive))
 {
   img = img.mirrored(false,true);
   img.save(imageName);
 }
}

void dumpPolySet(const char * imageName, vector< vector<Point2f> > &polyVec, vector<Similarity2f> &trVec, int width, int height)
{
  vector< vector< vector<Point2f> > > polyVecVec(polyVec.size());
  for(size_t i=0;i<polyVec.size();++i)
  {
    polyVecVec[i].resize(1);
    polyVecVec[i][0]=polyVec[i];
  }
  dumpPolySet(imageName,polyVecVec,trVec,width,height);
}

