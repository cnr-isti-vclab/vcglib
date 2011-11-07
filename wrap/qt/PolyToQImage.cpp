#include <QImage>
#include <QPainter>

#include "PolyToQImage.h"
#include <wrap/qt/col_qt_convert.h>
using namespace vcg;
using namespace std;

int dumpPolySet(const char * imageName,vector< vector<Point2f> > &polyVec, vector<Similarity2f> &trVec, int width, int height)
{
  QImage img(width,height,QImage::Format_RGB32);
  img.fill(qRgb(128,128,128));
  QPainter painter(&img);           // paint in picture
  for(size_t i=0;i<polyVec.size();++i)
  {
    QVector<QPointF> ppQ;
    for(int j=0;j<polyVec[i].size();++j)
    {
      Point2f pp=trVec[i]*polyVec[i][j];
      ppQ.push_back(QPointF(pp[0],pp[1]));
    }

    QBrush bb(vcg::ColorConverter::ToQColor(Color4b::Scatter(polyVec.size(),i)));
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

