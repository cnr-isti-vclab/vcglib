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

void dumpPolySet(const char * imageName, vector< vector<Point2f> > &polyVec, PolyDumperParam &pp)
{
  Box2f bb;
  for(size_t i=0;i<polyVec.size();++i)
    for(int j=0;j<polyVec[i].size();++j)
      bb.Add(polyVec[i][j]);

  Similarity2f sim;
  sim.sca = min(float(pp.widthPx)/bb.DimX(),float(pp.heightPx)/bb.DimY());
  sim.tra = Point2f(pp.widthPx/2.0f,pp.heightPx/2.0f)-bb.Center()*sim.sca;
  vector<Similarity2f> trVec(polyVec.size(),sim);

  dumpPolySet(imageName,polyVec,trVec,pp);
}


void dumpPolySet(const char * imageName, vector< vector<Point2f> > &polyVec, vector<Similarity2f> &trVec, PolyDumperParam &pp)
{
  vector< vector< vector<Point2f> > > polyVecVec(polyVec.size());
  for(size_t i=0;i<polyVec.size();++i)
  {
    polyVecVec[i].resize(1);
    polyVecVec[i][0]=polyVec[i];
  }
  dumpPolySet(imageName,polyVecVec,trVec,pp);
}

void dumpPolySet(const char * imageName, vector< vector< vector<Point2f> > > &polyVecVec, vector<Similarity2f> &trVec, PolyDumperParam &pp)
{
  vector<string> labelVec(polyVecVec.size());
  dumpPolySet(imageName,polyVecVec,trVec,labelVec,pp);
}

void dumpPolySet(const char * imageName, vector< vector< vector<Point2f> > > &polyVecVec, vector<Similarity2f> &trVec, vector<string> &labelVec, PolyDumperParam &pp)
{
  assert(polyVecVec.size() == trVec.size());
  QFont qf("courier",1);
  QSvgGenerator svg;
  svg.setFileName(imageName);
//  pp.widthPx = pp.widthDot;
//  pp.heightPx = pp.heightDot;
  svg.setSize(QSize(pp.widthMm,pp.heightMm));
  svg.setResolution(pp.dpi); // ?? dpm or dpi

  QImage img(pp.widthPx,pp.heightPx,QImage::Format_RGB32);
  img.fill(vcg::ColorConverter::ToQColor( pp.backgroundColor).rgb());

  QPainter painter;

  if(QString(imageName).endsWith("svg",Qt::CaseInsensitive))
        painter.begin(&svg);
  else  painter.begin(&img);
  QBrush bb;
  QPen qp;
  qp.setWidth(0);
//  printf("polyVecVec.size() = %i \n",polyVecVec.size());
  for(size_t i=0;i<polyVecVec.size();++i)
  {
    painter.resetTransform();
    painter.translate(trVec[i].tra[0],trVec[i].tra[1]);
    painter.rotate(math::ToDeg(trVec[i].rotRad));
    painter.scale(trVec[i].sca,trVec[i].sca);
    QPainterPath QPP;
    Point2f bc(0.f,0.f);
//    printf("polyVecVec[i].size() = %i \n",polyVecVec[i].size());
    for(int jj=0;jj<polyVecVec[i].size();++jj)
    {
      QVector<QPointF> ppQ;
      for(int j=0;j<polyVecVec[i][jj].size();++j)
      {
        Point2f pp=polyVecVec[i][jj][j];
        bc+=pp;
        ppQ.push_back(QPointF(pp[0],pp[1]));
      }
       ppQ.push_back(QPointF(polyVecVec[i][jj][0][0],polyVecVec[i][jj][0][1]));
      QPP.addPolygon(QPolygonF(ppQ));
    }
    bc/=float(polyVecVec[i][0].size());
    bb.setColor(vcg::ColorConverter::ToQColor(Color4b::Scatter(polyVecVec.size(),i)));

    painter.setBrush(bb);
    painter.setPen(qp);
    painter.drawPath(QPP);

//    painter.setFont(qf);
//    painter.drawText(bc[0],bc[1],QString(labelVec[i].c_str()));
  }
  painter.end();

// if(!QString(imageName).endsWith("svg",Qt::CaseInsensitive))
 {
   img = img.mirrored(false,true);
   img.save(imageName);
 }
}

