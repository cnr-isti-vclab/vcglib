/****************************************************************************
* VCGLib                                                            o o     *
* Visual and Computer Graphics Library                            o     o   *
*                                                                _   O  _   *
* Copyright(C) 2004                                                \/)\/    *
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
/****************************************************************************
  History

$Log: not supported by cvs2svn $
Revision 1.1  2004/09/09 22:38:57  cignoni
Initial Update

****************************************************************************/

#ifndef __VCG_GEN_NORMAL
#define __VCG_GEN_NORMAL
#include <wrap/gl/math.h>
namespace vcg {

  
template <class ScalarType>
class GenNormal 
{
public:
typedef Point3<ScalarType> Point3x;

static void  Random(int vn, std::vector<Point3<ScalarType > > &NN)
{
  NN.clear();
  while(NN.size()<vn)
  {
    Point3x pp(float(rand())/RAND_MAX,float(rand())/RAND_MAX,float(rand())/RAND_MAX);
    pp=pp*2.0-Point3x(1,1,1);
    if(pp.SquaredNorm()<1)
    {
      Normalize(pp);
      NN.push_back(pp);
    }
  }
}
static void UniformCone(int vn, std::vector<Point3<ScalarType > > &NN, ScalarType AngleRad, Point3x dir=Point3x(0,1,0))
{
  std::vector<Point3<ScalarType > > NNT;
  NN.clear();
  // per prima cosa si calcola il volume della spherical cap di angolo AngleRad
  ScalarType Height= 1.0 - cos(AngleRad); // height is measured from top...
  // Surface is the one of the tangent cylinder
  ScalarType CapArea = 2.0*M_PI*Height;
  ScalarType Ratio = CapArea / (4.0*M_PI );

  printf("----------AngleRad %f Angledeg %f ratio %f vn %i vn2 %i \n",AngleRad,math::ToDeg(AngleRad),Ratio,vn,int(vn/Ratio));
  Uniform(vn/Ratio,NNT);
  std::vector<Point3x>::iterator vi;
  
  ScalarType DotProd = cos(AngleRad);
  for(vi=NNT.begin();vi!=NNT.end();++vi)
  {
    if(dir*(*vi) >= DotProd) NN.push_back(*vi);
  }
 }


static void Uniform(int vn, std::vector<Point3<ScalarType > > &NN)
{
  OctaLevel pp;

  int ll=10;
  while(pow(4,ll)+2>vn) ll--;

  pp.Init(ll);
  sort(pp.v.begin(),pp.v.end());
  int newsize = unique(pp.v.begin(),pp.v.end())-pp.v.begin();
  pp.v.resize(newsize);

  NN=pp.v;
  Perturb(NN);
 }

static void Perturb(std::vector<Point3<ScalarType > > &NN)
{
  float width=0.2f/sqrt(float(NN.size()));

  for(vector<Point3x>::iterator vi=NN.begin(); vi!=NN.end();++vi) 
  {
    Point3x pp(float(rand())/RAND_MAX,float(rand())/RAND_MAX,float(rand())/RAND_MAX);
    pp=pp*2.0-Point3x(1,1,1);
    pp*=width;
    (*vi)+=pp;
    (*vi).Normalize();
  }

}

private :
class OctaLevel
  {
  public:
    std::vector<Point3x> v;
    int level;
    int sz;

    Point3f &Val(int i, int j) {
      assert(i>=0 && i<sz);
      assert(j>=0 && j<sz);
      return v[i+j*sz];
    }

    void Init(int lev)
    {
      sz=pow(2,lev+1)+1;
      v.resize(sz*sz);
      if(lev==0)
      {
        Val(0,0)=Point3x( 0, 0,-1); Val(0,1)=Point3x( 0, 1, 0);  Val(0,2)=Point3x( 0, 0,-1);
        
        Val(1,0)=Point3x(-1, 0, 0); Val(1,1)=Point3x( 0, 0, 1);   Val(1,2)=Point3x( 1, 0, 0);

        Val(2,0)=Point3x( 0, 0,-1); Val(2,1)=Point3x( 0,-1, 0);   Val(2,2)=Point3x( 0, 0,-1);
      }
      else
      {
        OctaLevel tmp;
        tmp.Init(lev-1);
        int i,j;
        for(i=0;i<sz;++i)
          for(j=0;j<sz;++j)
          {
            if((i%2)==0 && (j%2)==0) 
              Val(i,j)=tmp.Val(i/2,j/2);
            if((i%2)!=0 && (j%2)==0) 
                Val(i,j)=(tmp.Val(i/2+0,j/2)+tmp.Val(i/2+1,j/2))/2.0;
            if((i%2)==0 && (j%2)!=0) 
                Val(i,j)=(tmp.Val(i/2,j/2+0)+tmp.Val(i/2,j/2+1))/2.0;
            if((i%2)!=0 && (j%2)!=0) 
                Val(i,j)=(tmp.Val(i/2+0,j/2+0)+tmp.Val(i/2+0,j/2+1)+tmp.Val(i/2+1,j/2+0)+tmp.Val(i/2+1,j/2+1))/4.0;
          }
        for(vector<Point3x>::iterator vi=v.begin(); vi!=v.end();++vi) 
            (*vi).Normalize();
  
      }
    }
  };
};
}
#endif
