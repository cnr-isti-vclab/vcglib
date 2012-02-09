/****************************************************************************
* VCGLib                                                            o o     *
* Visual and Computer Graphics Library                            o     o   *
*                                                                _   O  _   *
* Copyright(C) 2004-2009                                           \/)\/    *
* Visual Computing Lab                                            /\/|      *
* ISTI - Italian National Research Council                           |      *
*                                                                    \      *
* All rights reserved.                                                      *
*                                                                           *
* This program is free software; you can redistribute it and/or modify      *
* it under the terms of the GNU General Public License as published by      *
* the Free Software Foundation; either version 2 of the License, or         *
* (at your option) any later version.                                       *
*                                                                           *
* This program is distributed in the hope that it will be useful,           *
* but WITHOUT ANY WARRANTY; without even the implied warranty of            *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
* GNU General Public License (http://www.gnu.org/licenses/gpl.txt)          *
* for more details.                                                         *
*                                                                           *
****************************************************************************/
#include <QtOpenGL/QtOpenGL>
#include<vcg/space/box2.h>
#include<vcg/math/random_generator.h>
#include<wrap/qt/col_qt_convert.h>
#include <vcg/space/rect_packer.h>
#include <vcg/space/poly_packer.h>
#include <wrap/qt/PolyToQImage.h>
#include <time.h>

using namespace vcg;
using namespace std;

static void buildRandRectSet(int rectNum, vector<Box2f> &rectVec)
{
  math::MarsenneTwisterRNG rnd;
  float exp=3.0f;
  float ratioMin=0.2;
  float ratioMax=0.9;
  float sizeMin=0.1;
  float sizeMax=1.0f;
  rnd.initialize(time(0));
  for(int i=0;i<rectNum;++i)
  {
    Box2f bb;
    float ratio=ratioMin+(ratioMax-ratioMin)*rnd.generate01();
    float size= sizeMin+(sizeMax-sizeMin)*pow((float)rnd.generate01(),exp);
    bb.min=Point2f(-size*ratio,-size);
    bb.max=Point2f( size*ratio, size);
    rectVec.push_back(bb);
  }
}

void buildRandRectSetOld(int rectNum, vector<Box2f> &rectVec)
{
  math::MarsenneTwisterRNG rnd;
  float exp=5.0f;
  rnd.initialize(time(0));
  for(int i=0;i<rectNum;++i)
  {
    Box2f bb;
    bb.min=Point2f(-pow((float)rnd.generate01(),exp),-pow((float)rnd.generate01(),exp));
    bb.max=Point2f( pow((float)rnd.generate01(),exp), pow((float)rnd.generate01(),exp));
    rectVec.push_back(bb);
  }
}

void buildRandPolySet(int polyNum, vector< vector<Point2f> > &polyVec)
{
  vcg::math::MarsenneTwisterRNG rnd;
  rnd.initialize(time(0));

  for(int i=0;i<polyNum;++i)
  {
    vector<Point2f> poly;
    for(int j=0;j<10;j++)
      poly.push_back(Point2f(0.5+0.5*rnd.generate01(),2.0f*M_PI*rnd.generate01()));

    std::sort(poly.begin(),poly.end());

    float ratio = rnd.generateRange(0.2,0.9);
    float rot = rnd.generateRange(-M_PI,M_PI);
    float scale = pow(rnd.generateRange(0.3,0.9),1);

    for(size_t j=0;j<poly.size();j++)
    {
      poly[j].Polar2Cartesian();
      poly[j][1]*=ratio;
      poly[j] *= scale;
      poly[j].Cartesian2Polar();
      poly[j][1]+=rot;
      poly[j].Polar2Cartesian();
    }
    polyVec.push_back(poly);
  }
}

int main( int argc, char **argv )
{
  QApplication pippo(argc,argv);

  vector<Box2f> rectVec;
  buildRandRectSet(1000, rectVec);
  vector<Similarity2f> trVec;
  vector< vector<Point2f> > polySet;
  Point2f finalSize;
  buildRandPolySet(100,polySet);
  PolyDumperParam pp;
 /* PolyPacker<float>::PackAsEqualSquares(polySet,Point2f(1024.0f,1024.0f),trVec,finalSize);
  dumpPolySet("testpolyEq.png",polySet,trVec,pp);
  PolyPacker<float>::PackAsAxisAlignedRect(polySet,Point2f(1024.0f,1024.0f),trVec,finalSize);
  dumpPolySet("testpolyAA.png",polySet,trVec,pp);
  PolyPacker<float>::PackAsObjectOrientedRect(polySet,Point2f(1024.0f,1024.0f),trVec,finalSize);
  dumpPolySet("testpolyOO.png",polySet,trVec,pp);*/

  //PolyPacker<float>::PackAsAxisAlignedRect(polySet,Point2f(1024.0f,1024.0f),trVec,finalSize);
  PolyPacker<float>::PackAsObjectOrientedRect(polySet,Point2f(1024.0f,1024.0f),trVec,finalSize);
  //dumpPolySetPNG("testpolyEq.png",polySet,trVec,pp);
  PolyDumper::dumpPolySetSVG("testpolyEq.svg",polySet,trVec,pp);

  /*PolyPacker<float>::PackAsAxisAlignedRect(polySet,Point2f(1024.0f,1024.0f),trVec,finalSize);
  dumpPolySetSVG("testpolyAA.svg",polySet,trVec,pp);
  PolyPacker<float>::PackAsObjectOrientedRect(polySet,Point2f(1024.0f,1024.0f),trVec,finalSize);
  dumpPolySetSVG("testpolyOO.svg",polySet,trVec,pp);*/

  return 0;
}
