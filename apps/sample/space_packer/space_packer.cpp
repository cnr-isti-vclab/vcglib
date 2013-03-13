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
#include<vcg/space/box3.h>
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
  vector<Similarity2f> trVec;
  vector<Similarity2f> trPolyVec;
  vector< vector<Point2f> > polySet;
  vector< vector<Point2f> > multiPolySet;
  Point2f finalSize;
  std::vector<Point2f> finalSizeVec;
  const Point2f containerSize(1000,1000);
    PolyDumperParam pp;
  std::vector<int> contInd;

  vector<Box2f> rectVec;
  buildRandRectSet(10,rectVec);
//  RectPacker<float>::Pack(rectVec,containerSize,trVec,finalSize);
  RectPacker<float>::PackMulti(rectVec,containerSize,3,trVec,contInd,finalSizeVec);
  RectPacker<float>::Stat s = RectPacker<float>::stat();
  printf("RectPacker attempt %i time %5.3f %5.3f\n",s.pack_attempt_num,s.pack_total_time,s.pack_attempt_time);

//  PolyDumper::rectSetToPolySet(rectVec,polySet);

//  PolyDumper::multiRectSetToSinglePolySet(rectVec,trVec,contInd,0,polySet,trPolyVec);
//  PolyDumper::dumpPolySetPNG("testpolyEq0.png",polySet,trPolyVec,pp);
//  PolyDumper::multiRectSetToSinglePolySet(rectVec,trVec,contInd,1,polySet,trPolyVec);
//  PolyDumper::dumpPolySetPNG("testpolyEq1.png",polySet,trPolyVec,pp);
//  PolyDumper::multiRectSetToSinglePolySet(rectVec,trVec,contInd,2,polySet,trPolyVec);
//  PolyDumper::dumpPolySetPNG("testpolyEq2.png",polySet,trPolyVec,pp);


//   buildRandPolySet(100,polySet);
//   PolyPacker<float>::PackMultiAsObjectOrientedRect(polySet,containerSize,3,trVec,contInd,finalSizeVec);

//   PolyDumper::multiPolySetToSinglePolySet(polySet,trVec,contInd,0,multiPolySet,trPolyVec);
//   PolyDumper::dumpPolySetPNG("testpolyEq0.png",multiPolySet,trPolyVec,pp);

//   PolyDumper::multiPolySetToSinglePolySet(polySet,trVec,contInd,1,multiPolySet,trPolyVec);
//   PolyDumper::dumpPolySetPNG("testpolyEq1.png",multiPolySet,trPolyVec,pp);

//   PolyDumper::multiPolySetToSinglePolySet(polySet,trVec,contInd,2,multiPolySet,trPolyVec);
//   PolyDumper::dumpPolySetPNG("testpolyEq2.png",multiPolySet,trPolyVec,pp);

   //  PolyDumper::dumpPolySetPNG("testpolyOO.png",polySet,trVec,pp);


  buildRandPolySet(100,polySet);

  PolyPacker<float>::PackAsEqualSquares(polySet,containerSize,trVec,finalSize);
  PolyDumper::dumpPolySetPNG("testpolyEq.png",polySet,trVec,pp);

  PolyPacker<float>::PackAsAxisAlignedRect(polySet,containerSize,trVec,finalSize);
  PolyDumper::dumpPolySetPNG("testpolyAA.png",polySet,trVec,pp);

  PolyPacker<float>::PackAsObjectOrientedRect(polySet,containerSize,trVec,finalSize);
  PolyDumper::dumpPolySetPNG("testpolyOO.png",polySet,trVec,pp);

  return 0;
}
