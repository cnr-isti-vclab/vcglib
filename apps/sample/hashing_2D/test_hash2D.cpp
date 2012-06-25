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
#include <stdio.h>
#include <time.h>
#include <vcg/complex/used_types.h>
#include <vcg/space/distance2.h>
#include<vcg/space/index/index2D/spatial_hashing_2D.h>
#include<vcg/space/intersection2.h>

typedef double MyScalarType;

class MySegmentType:public vcg::Segment2<MyScalarType>
{
public:
	int mark;
	bool deleted;
	bool IsD(){return deleted;}

	MySegmentType(const vcg::Point2<MyScalarType> &_P0,
				  const vcg::Point2<MyScalarType> &_P1)
	{
		P0()=_P0;
		P1()=_P1;
		mark=0;
	}

	void GetBBox(vcg::Box2<ScalarType> &BB2)
	{
		//BB2.SetNull();
		BB2.Set(P0());
		BB2.Add(P1());
	}

	MySegmentType(){}

	MySegmentType(const MySegmentType &s1)
	{
		P0()=s1.P0();
		P1()=s1.P1();
		mark=s1.mark;
		deleted=s1.deleted;
	}
};

//**MARKER CLASSES**//
class MyMarker
{
	
public:
	int mark;

	MyMarker(){mark=0;}
	//MyMarker(	MESH_TYPE *m) {SetMesh(m);}
	void UnMarkAll(){mark++;}

	bool IsMarked(MySegmentType* obj)
	{return(obj->mark==mark);}

	void Mark(MySegmentType* obj)
	{obj->mark=mark;}
	/*void SetMesh(MESH_TYPE *_m)
	{m=_m;}*/
};


vcg::SpatialHashTable2D<MySegmentType,MyScalarType> Hash2D;
std::vector<MySegmentType> Allocated;
MyMarker MyMark;

void RandomSeg(vcg::Point2<MyScalarType> &P0,
				vcg::Point2<MyScalarType> &P1,
				MyScalarType SpaceSize=100,
				MyScalarType maxdim=0.01)
{
	MyScalarType dimAbs=SpaceSize*maxdim;
	int dimension=RAND_MAX;

	int X=rand();
	int Y=rand();
	int dX=rand();
	int dY=rand();
	MyScalarType size=((MyScalarType)(rand()))/(MyScalarType)dimension;

	P0=vcg::Point2<MyScalarType>((MyScalarType)X/dimension,(MyScalarType)Y/dimension);
	P0*=SpaceSize;

	vcg::Point2<MyScalarType> D=vcg::Point2<MyScalarType>((MyScalarType)dX/dimension,(MyScalarType)dY/dimension);
	D.Normalize();
	D*=size*dimAbs;
	P1=P0+D;
}

void InitRandom(int num,
				MyScalarType SpaceSize=100,
				MyScalarType maxdim=0.01)
{
	Allocated.clear();
	Allocated.resize(num);
	srand(clock());
	for (int i=0;i<num;i++)
	{
		vcg::Point2<MyScalarType> P0,P1;
		RandomSeg(P0,P1,SpaceSize,maxdim);
		Allocated[i]=MySegmentType(P0,P1);
		Allocated[i].deleted=false;
	}

}


MyScalarType TestBox(int num_test=100000,
				 MyScalarType SpaceSize=100,
				MyScalarType maxdim=0.02)
{
	//GetInBox(OBJMARKER & _marker,const Box2x _bbox,OBJPTRCONTAINER & _objectPtrs)
	MyMark.UnMarkAll();
	int t0=clock();
	int num=0;
	for (int i=0;i<num_test;i++)
	{
		vcg::Point2<MyScalarType> P0,P1;
		RandomSeg(P0,P1,SpaceSize,maxdim);
		vcg::Box2<MyScalarType> bbox;
		bbox.Add(P0);
		bbox.Add(P1);
		std::vector<MySegmentType*> result;
		num+=Hash2D.GetInBox<MyMarker,std::vector<MySegmentType*> >(MyMark,bbox,result);
	}
	int t1=clock();
	MyScalarType numd=(double)num/(double)num_test;
	return numd;
}

//void TestCorrectNess(int num_sample=1000,
//					int num_test=1000,
//					MyScalarType SpaceSize=100,
//					MyScalarType radius=0.02)
//{
//
//}

void GetIntersectingSegments(MySegmentType *S,
							std::vector<MySegmentType*> &result)
{
	///get the bbox
	result.clear();
	vcg::Box2<MyScalarType> bbox;
	S->GetBBox(bbox);
	///then get into the grid
	std::vector<MySegmentType*> inbox;
	int num=Hash2D.GetInBox<MyMarker,std::vector<MySegmentType*> >(MyMark,bbox,inbox);
	///then test intersection
	for (int j=0;j<num;j++)
	{
		if (inbox[j]==S)continue;
		vcg::Point2<MyScalarType> p_inters;
		if (vcg::SegmentSegmentIntersection<MyScalarType>(*S,*inbox[j],p_inters))
			result.push_back(inbox[j]);
	}
}

void GetCloseSegments(MySegmentType *S,
					const MyScalarType &radius,
					std::vector<MySegmentType*> &result)
{
	///get the bbox
	result.clear();
	vcg::Box2<MyScalarType> bbox;
	S->GetBBox(bbox);
	bbox.Offset(radius);//*1.02);

	///then get into the grid
	std::vector<MySegmentType*> inbox;
	int num=Hash2D.GetInBox<MyMarker,std::vector<MySegmentType*> >(MyMark,bbox,inbox);
	///then test intersection
	for (int j=0;j<num;j++)
	{
		if (inbox[j]==S)continue;
		vcg::Point2<MyScalarType> p_clos;
		MyScalarType dist=vcg::Segment2DSegment2DDistance<MyScalarType>(*S,*inbox[j],p_clos);
		if (dist<radius) result.push_back(inbox[j]);
	}
}

MyScalarType TestIntersection(int num_test=1000000)
{
	int num_inters=0;
	for (int i=0;i<num_test;i++)
	{
		assert(i<Allocated.size());
		std::vector<MySegmentType*> result;
		GetIntersectingSegments(&Allocated[i],result);
		num_inters+=result.size();
	}
	return ((MyScalarType)num_inters/(MyScalarType)num_test);
}

MyScalarType TestClosest(int num_test=1000000,
						MyScalarType radius=0.1)
{
	int num_found=0;
	for (int i=0;i<num_test;i++)
	{
		assert(i<Allocated.size());

		//get the segment
		MySegmentType *S=&Allocated[i];
		
		MyScalarType absRadius=S->Length()*radius;
		
		///get the segments closer than a radius
		std::vector<MySegmentType*> closer;
		GetCloseSegments(S,absRadius,closer);
		
		num_found+=closer.size();
	}
	return ((MyScalarType)num_found/(MyScalarType)num_test);
}

int TestCorrectIntersect(int num_test=1000)
{
	int num=0;
	for (int i=0;i<num_test;i++)
	{
		MySegmentType S0=Allocated[i];
		std::vector<MySegmentType*> result0,result1;
		for (int j=0;j<num_test;j++)
		{
			if (j==i) continue;
			MySegmentType *S1=&Allocated[j];
			vcg::Point2<MyScalarType> p_inters;
			if (vcg::SegmentSegmentIntersection<MyScalarType>(S0,*S1,p_inters))
				result0.push_back(S1);
			/*num+=result0.size();*/
		}
		GetIntersectingSegments(&Allocated[i],result1);
		///then see if equal number
		if (result1.size()==result0.size())num++;
	}
	return (num);
}

int TestCorrectClosest(int num_test=1000,
						MyScalarType radius=0.1)
{
	int num=0;
	for (int i=0;i<num_test;i++)
	{
		MySegmentType *S0=&Allocated[i];
		std::vector<MySegmentType*> result0,result1;
		MyScalarType absRadius=S0->Length()*radius;
		for (int j=0;j<num_test;j++)
		{
			if (j==i) continue;
			MySegmentType *S1=&Allocated[j];
			vcg::Point2<MyScalarType> p_clos;
			MyScalarType dist=vcg::Segment2DSegment2DDistance<MyScalarType>(*S0,*S1,p_clos);
			if (dist<absRadius)
				result0.push_back(S1);
			/*num+=result0.size();*/
		}
		GetCloseSegments(S0,absRadius,result1);
		///then see if equal number
		if (result1.size()==result0.size())num++;
	}
	return (num);
}


int main( int argc, char **argv )
{
  int num_sample=100000;
  int t0=clock();
  InitRandom(100000);
  int t1=clock();

  ///Initialization performance
  printf("** Time elapsed for initialization of %d sample is %d\n \n",num_sample,t1-t0);
  Hash2D.Set(Allocated.begin(),Allocated.end());

  ///Box Query performance
  t0=clock();
  MyScalarType avg_test=TestBox(num_sample);
  t1=clock();
  printf("** Time elapsed for %d BOX queries is %d\n, average found %5.5f \n \n",num_sample,t1-t0,avg_test);
  

  ///Intersecting segment performance
  t0=clock();
  MyScalarType avg_int=TestIntersection(num_sample);
  t1=clock();
  printf("** Time elapsed for %d INTERSECTION queries is %d\n, average found %5.5f \n \n",num_sample,t1-t0,avg_int);
	
  ///closest test
  t0=clock();
  MyScalarType avg_clos=TestClosest(num_sample);
  t1=clock();
  printf("** Time elapsed for %d CLOSEST queries is %d\n, average found %5.5f \n \n",num_sample,t1-t0,avg_clos);
	
 ///reinitialize structure
  MyMark.mark=0;
  Hash2D.Clear();
  int n_test=1000;
  InitRandom(n_test,100,0.1);
  Hash2D.Set(Allocated.begin(),Allocated.end());

  int tested_int=TestCorrectIntersect(n_test);
  printf("** Correct Intersect on %d test are %d \n",n_test,tested_int);
	
  int tested_clos=TestCorrectClosest(n_test);
  printf("** Correct Closest on %d test are %d \n",n_test,tested_clos);

  return 0;
}
