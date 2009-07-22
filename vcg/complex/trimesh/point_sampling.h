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

The sampling Class has a set of static functions, that you can call to sample the surface of a mesh.
Each function is templated on the mesh and on a Sampler object s. 
Each function calls many time the sample object with the sampling point as parameter.
 
Sampler Classes and Sampling algorithms are independent. 
Sampler classes exploits the sample that are generated with various algorithms.
For example, you can compute Hausdorff distance (that is a sampler) using various 
sampling strategies (montecarlo, stratified etc).
 
****************************************************************************/
#ifndef __VCGLIB_POINT_SAMPLING
#define __VCGLIB_POINT_SAMPLING

#include <vcg/math/random_generator.h>
#include <vcg/complex/trimesh/closest.h>
#include <vcg/space/index/spatial_hashing.h>
#include <vcg/complex/trimesh/stat.h>
#include <vcg/complex/trimesh/update/topology.h>
#include <vcg/space/box2.h>
namespace vcg
{
namespace tri
{

/// Trivial Sampler, an example sampler object that show the required interface used by the sampling class. 
/// Most of the sampling classes call the AddFace method with the face containing the sample and its barycentric coord.
/// Beside being an example of how to write a sampler it provides a simple way to use the various sampling classes. 
// For example if you just want to get a vector with positions over the surface You have just to write
//
// vector<Point3f> myVec;
// TrivialSampler<MyMesh> ts(myVec) 
// SurfaceSampling<MyMesh, TrivialSampler<MyMesh> >::Montecarlo(M, ts, SampleNum);
// 
//

template <class MeshType>
class TrivialSampler
{
	public:
		typedef typename MeshType::CoordType			CoordType;
		typedef typename MeshType::VertexType			VertexType;
    typedef typename MeshType::FaceType				FaceType;

	TrivialSampler()
	{
		sampleVec = new std::vector<CoordType>();
		vectorOwner=true;
	};

	TrivialSampler(std::vector<CoordType> &Vec)
	{
		sampleVec = &Vec;
		sampleVec->clear();
		vectorOwner=false;
	};

	~TrivialSampler()
	{
		if(vectorOwner) delete sampleVec;
	}
	
	private:
		std::vector<CoordType> *sampleVec;
		bool vectorOwner;
	public:
	
	void AddVert(const VertexType &p) 
	{
		sampleVec->push_back(p.cP());
	}
	void AddFace(const FaceType &f, const CoordType &p) 
	{
		sampleVec->push_back(f.P(0)*p[0] + f.P(1)*p[1] +f.P(2)*p[2] );
	}
	
	void AddTextureSample(const FaceType &, const CoordType &, const Point2i &)
	{
		// Retrieve the color of the sample from the face f using the barycentric coord p 
		// and write that color in a texture image at position tp[0],tp[1]
	}
}; // end class TrivialSampler

template <class MetroMesh, class VertexSampler>
class SurfaceSampling
{
		typedef typename MetroMesh::CoordType			CoordType;
		typedef typename MetroMesh::ScalarType			ScalarType;
		typedef typename MetroMesh::VertexType			VertexType;
		typedef typename MetroMesh::VertexPointer		VertexPointer;
		typedef typename MetroMesh::VertexIterator		VertexIterator;
		typedef typename MetroMesh::FacePointer			FacePointer;
		typedef typename MetroMesh::FaceIterator		FaceIterator;
		typedef typename MetroMesh::FaceType			FaceType;
		typedef typename MetroMesh::FaceContainer		FaceContainer;

		typedef typename vcg::SpatialHashTable<FaceType, ScalarType> MeshSHT;
		typedef typename vcg::SpatialHashTable<FaceType, ScalarType>::CellIterator MeshSHTIterator;
		typedef typename vcg::SpatialHashTable<VertexType, ScalarType> MontecarloSHT;
		typedef typename vcg::SpatialHashTable<VertexType, ScalarType>::CellIterator MontecarloSHTIterator;
		typedef typename vcg::SpatialHashTable<VertexType, ScalarType> SampleSHT;
		typedef typename vcg::SpatialHashTable<VertexType, ScalarType>::CellIterator SampleSHTIterator;

public:

static math::MarsenneTwisterRNG &SamplingRandomGenerator() 
{
	static math::MarsenneTwisterRNG rnd;
	return rnd;
}

// Returns an integer random number in the [0,i-1] interval using the improve Marsenne-Twister method.
static unsigned int RandomInt(unsigned int i)
{
	return (SamplingRandomGenerator().generate(0) % i);
}

// Returns a random number in the [0,1) real interval using the improved Marsenne-Twister method.
static double RandomDouble01()
{
	return SamplingRandomGenerator().generate01();
}

// Returns a random number in the [0,1] real interval using the improved Marsenne-Twister.
static double RandomDouble01closed()
{
	return SamplingRandomGenerator().generate01closed();
}

static void AllVertex(MetroMesh & m, VertexSampler &ps)
{
	VertexIterator vi;
	for(vi=m.vert.begin();vi!=m.vert.end();++vi)
	{
		if(!(*vi).IsD())
		{
			ps.AddVert(*vi);
		}
	}
}

/// Sample the vertices in a weighted way. Each vertex has a probability of being chosen
/// that is proportional to its quality. 
/// It assumes that you are asking a number of vertices smaller than nv;
/// Algorithm: 
/// 1) normalize quality so that sum q == 1;
/// 2) shuffle vertices.
/// 3) for each vertices choose it if rand > thr;
 
static void VertexWeighted(MetroMesh & m, VertexSampler &ps, int sampleNum)
{
	ScalarType qSum = 0;
	VertexIterator vi;
	for(vi = m.vert.begin(); vi != m.vert.end(); ++vi)
			if(!(*vi).IsD()) 
						qSum += (*vi).Q();
	  
	ScalarType samplePerUnit = sampleNum/qSum;
	ScalarType floatSampleNum =0;
	std::vector<VertexPointer> vertVec;
	FillAndShuffleVertexPointerVector(m,vertVec);

	std::vector<bool> vertUsed(m.vn,false);
	
	int i=0; int cnt=0;
	while(cnt < sampleNum)
		{
			if(vertUsed[i])
				{
						floatSampleNum += vertVec[i]->Q() * samplePerUnit;
						int vertSampleNum   = (int) floatSampleNum;
						floatSampleNum -= (float) vertSampleNum; 
						
						// for every sample p_i in T...
						if(vertSampleNum > 1)
							{
								ps.AddVert(*vertVec[i]);
								cnt++;
								vertUsed[i]=true;
							}
				}
			i = (i+1)%m.vn;				
		}
}

/// Sample the vertices in a uniform way. Each vertex has a probability of being chosen
/// that is proportional to the area it represent. 
static void VertexAreaUniform(MetroMesh & m, VertexSampler &ps, int sampleNum)
{
	VertexIterator vi;
	for(vi = m.vert.begin(); vi != m.vert.end(); ++vi)
		if(!(*vi).IsD()) 
						(*vi).Q() = 0;

	FaceIterator fi;
	for(fi = m.face.begin(); fi != m.face.end(); ++fi)
		if(!(*fi).IsD()) 
		{
			ScalarType areaThird = DoubleArea(*fi)/6.0;
			(*fi).V(0).Q()+=areaThird;
			(*fi).V(1).Q()+=areaThird;
			(*fi).V(2).Q()+=areaThird;
		}
	
	VertexWeighted(m,ps,sampleNum);
}
	
static void	FillAndShuffleFacePointerVector(MetroMesh & m, std::vector<FacePointer> &faceVec)
{
	FaceIterator fi;
	for(fi=m.face.begin();fi!=m.face.end();++fi) 
		if(!(*fi).IsD())	faceVec.push_back(&*fi);
	
	assert((int)faceVec.size()==m.fn);
	
	unsigned int (*p_myrandom)(unsigned int) = RandomInt;
	std::random_shuffle(faceVec.begin(),faceVec.end(), p_myrandom);
}
static void	FillAndShuffleVertexPointerVector(MetroMesh & m, std::vector<VertexPointer> &vertVec)
{
	VertexIterator vi;
	for(vi=m.vert.begin();vi!=m.vert.end();++vi) 
				if(!(*vi).IsD())	vertVec.push_back(&*vi);

	assert((int)vertVec.size()==m.vn);
	
	unsigned int (*p_myrandom)(unsigned int) = RandomInt;
	std::random_shuffle(vertVec.begin(),vertVec.end(), p_myrandom);
}

/// Sample the vertices in a uniform way. Each vertex has the same probabiltiy of being chosen. 
static void VertexUniform(MetroMesh & m, VertexSampler &ps, int sampleNum)
{
	if(sampleNum>=m.vn) {
	  AllVertex(m,ps);
		return;
	} 
	
	std::vector<VertexPointer> vertVec;
	FillAndShuffleVertexPointerVector(m,vertVec);
	
	for(int i =0; i< sampleNum; ++i)
		ps.AddVert(*vertVec[i]);
}


static void FaceUniform(MetroMesh & m, VertexSampler &ps, int sampleNum)
{
	if(sampleNum>=m.fn) {
	  AllFace(m,ps);
		return;
	} 

	std::vector<FacePointer> faceVec;
	FillAndShuffleFacePointerVector(m,faceVec);

	for(int i =0; i< sampleNum; ++i)
		ps.AddFace(*faceVec[i],Barycenter(*faceVec[i]));
}

static void AllFace(MetroMesh & m, VertexSampler &ps)
{
	FaceIterator fi;
	for(fi=m.face.begin();fi!=m.face.end();++fi) 
				if(!(*fi).IsD())
				{
					ps.AddFace(*fi,Barycenter(*fi));
				}
}


static void AllEdge(MetroMesh & m, VertexSampler &ps)
{
    // Edge sampling.
		typedef typename UpdateTopology<MetroMesh>::PEdge SimpleEdge;
		std::vector< SimpleEdge > Edges;
		typename std::vector< SimpleEdge >::iterator ei;
    UpdateTopology<MetroMesh>::FillUniqueEdgeVector(m,Edges);

		for(ei=Edges.begin(); ei!=Edges.end(); ++ei)
				{
					Point3f interp(0,0,0);
					interp[ (*ei).z     ]=.5; 
					interp[((*ei).z+1)%3]=.5;
					ps.AddFace(*(*ei).f,interp);
				}
}

// Regular Uniform Edge sampling
// Each edge is subdivided in a number of pieces proprtional to its lenght
// Sample are choosen without touching the vertices.

static void EdgeUniform(MetroMesh & m, VertexSampler &ps,int sampleNum)
{
		typedef typename UpdateTopology<MetroMesh>::PEdge SimpleEdge;
		std::vector< SimpleEdge > Edges;
    UpdateTopology<MetroMesh>::FillUniqueEdgeVector(m,Edges);
    // First loop compute total edge lenght;
		float edgeSum=0;
		typename std::vector< SimpleEdge >::iterator ei;
		for(ei=Edges.begin(); ei!=Edges.end(); ++ei)
					edgeSum+=Distance((*ei).v[0]->P(),(*ei).v[1]->P());
					
		//qDebug("Edges %i  edge sum %f",Edges.size(),edgeSum);			
		float sampleLen = edgeSum/sampleNum;
		//qDebug("EdgesSamples %i  Sampling Len %f",sampleNum,sampleLen);			
		float rest=0;
		for(ei=Edges.begin(); ei!=Edges.end(); ++ei)
				{
					float len = Distance((*ei).v[0]->P(),(*ei).v[1]->P());
					float samplePerEdge = floor((len+rest)/sampleLen);
					rest = (len+rest) - samplePerEdge * sampleLen;
					float step = 1.0/(samplePerEdge+1);
					for(int i=0;i<samplePerEdge;++i)
					{
						Point3f interp(0,0,0);
						interp[ (*ei).z     ]=step*(i+1); 
						interp[((*ei).z+1)%3]=1.0-step*(i+1);
						ps.AddFace(*(*ei).f,interp);
					}
				}
}

// Generate the barycentric coords of a random point over a single face, 
// with a uniform distribution over the triangle. 
// It uses the parallelogram folding trick. 
static CoordType RandomBaricentric()
{
	CoordType interp;
	interp[1] = RandomDouble01();
	interp[2] = RandomDouble01();
	if(interp[1] + interp[2] > 1.0)
	{
		interp[1] = 1.0 - interp[1];
		interp[2] = 1.0 - interp[2];
		}
	
	assert(interp[1] + interp[2] <= 1.0);
	interp[0]=1.0-(interp[1] + interp[2]);
	return interp;
}

static void StratifiedMontecarlo(MetroMesh & m, VertexSampler &ps,int sampleNum)
{
	ScalarType area = Stat<MetroMesh>::ComputeMeshArea(m);
	ScalarType samplePerAreaUnit = sampleNum/area;
	//qDebug("samplePerAreaUnit %f",samplePerAreaUnit);
	// Montecarlo sampling.
	double  floatSampleNum = 0.0;
	
	FaceIterator fi;	
	for(fi=m.face.begin(); fi != m.face.end(); fi++)
		if(!(*fi).IsD())
		{
			// compute # samples in the current face (taking into account of the remainders)
			floatSampleNum += 0.5*DoubleArea(*fi) * samplePerAreaUnit;
			int faceSampleNum   = (int) floatSampleNum;
			
			// for every sample p_i in T...
			for(int i=0; i < faceSampleNum; i++)
					ps.AddFace(*fi,RandomBaricentric());
			floatSampleNum -= (double) faceSampleNum; 
		}
}

static void Montecarlo(MetroMesh & m, VertexSampler &ps,int sampleNum)
{
	typedef  std::pair<ScalarType, FacePointer> IntervalType;
	std::vector< IntervalType > intervals (m.fn+1);
	FaceIterator fi;	
	int i=0;
	intervals[i]=std::make_pair(0,FacePointer(0));
	// First loop: build a sequence of consecutive segments proportional to the triangle areas.
	for(fi=m.face.begin(); fi != m.face.end(); fi++)
		if(!(*fi).IsD())
		{
			intervals[i+1]=std::make_pair(intervals[i].first+0.5*DoubleArea(*fi), &*fi);
			++i;
		}
	ScalarType meshArea = intervals.back().first;
	for(i=0;i<sampleNum;++i)
		{
			ScalarType val = meshArea * RandomDouble01();
			// lower_bound returns the furthermost iterator i in [first, last) such that, for every iterator j in [first, i), *j < value.
			// E.g. An iterator pointing to the first element "not less than" val, or end() if every element is less than val.
			typename std::vector<IntervalType>::iterator it = lower_bound(intervals.begin(),intervals.end(),std::make_pair(val,FacePointer(0)) );
			assert(it != intervals.end());
			assert(it != intervals.begin());
			assert( (*(it-1)).first <val );
			assert( (*(it)).first >= val);
			ps.AddFace( *(*it).second, RandomBaricentric() );
		}
}
	
static ScalarType WeightedArea(FaceType f)
{
	ScalarType averageQ = ( f.V(0)->Q() + f.V(1)->Q() + f.V(2)->Q() ) /3.0;
	return DoubleArea(f)*averageQ/2.0;
}

/// Compute a sampling of the surface that is weighted by the quality
/// the area of each face is multiplied by the average of the quality of the vertices. 
/// So the a face with a zero quality on all its vertices is never sampled and a face with average quality 2 get twice the samples of a face with the same area but with an average quality of 1;
static void WeightedMontecarlo(MetroMesh & m, VertexSampler &ps, int sampleNum)
{
	assert(tri::HasPerVertexQuality(m));
	
	ScalarType weightedArea = 0;
	FaceIterator fi;
	for(fi = m.face.begin(); fi != m.face.end(); ++fi)
			if(!(*fi).IsD()) 
						weightedArea += WeightedArea(*fi);
	
	ScalarType samplePerAreaUnit = sampleNum/weightedArea;
	//qDebug("samplePerAreaUnit %f",samplePerAreaUnit);
	// Montecarlo sampling.
	double  floatSampleNum = 0.0;
	for(fi=m.face.begin(); fi != m.face.end(); fi++)
		if(!(*fi).IsD())
    {
			// compute # samples in the current face (taking into account of the remainders)
			floatSampleNum += WeightedArea(*fi) * samplePerAreaUnit;
			int faceSampleNum   = (int) floatSampleNum;
			
			// for every sample p_i in T...
			for(int i=0; i < faceSampleNum; i++)
					ps.AddFace(*fi,RandomBaricentric());
							
			floatSampleNum -= (double) faceSampleNum; 
    }
}


// Subdivision sampling of a single face. 
// return number of added samples

static int SingleFaceSubdivision(int sampleNum, const CoordType & v0, const CoordType & v1, const CoordType & v2, VertexSampler &ps, FacePointer fp, bool randSample)
{
    // recursive face subdivision.
    if(sampleNum == 1)
    {
        // ground case.
			CoordType SamplePoint;
			if(randSample) 
			{
				CoordType rb=RandomBaricentric();
				SamplePoint=v0*rb[0]+v1*rb[1]+v2*rb[2];
			}
			else SamplePoint=((v0+v1+v2)/3.0f);
			
			CoordType SampleBary;
			InterpolationParameters(*fp,SamplePoint,SampleBary[0],SampleBary[1],SampleBary[2]);
			ps.AddFace(*fp,SampleBary);
			return 1;
    }

		int s0 = sampleNum /2;
		int s1 = sampleNum-s0;
		assert(s0>0);
		assert(s1>0);
	
	  ScalarType w0 = ScalarType(s1)/ScalarType(sampleNum);
		ScalarType w1 = 1.0-w0;
    // compute the longest edge.
    double  maxd01 = SquaredDistance(v0,v1);
    double  maxd12 = SquaredDistance(v1,v2);
    double  maxd20 = SquaredDistance(v2,v0);
    int     res;
    if(maxd01 > maxd12)
        if(maxd01 > maxd20)     res = 0;
        else                    res = 2;
    else
        if(maxd12 > maxd20)     res = 1;
        else                    res = 2;
				
    int faceSampleNum=0;
    // break the input triangle along the midpoint of the longest edge.
    CoordType  pp;
    switch(res)
    {
     case 0 :    pp = v0*w0 + v1*w1;
                 faceSampleNum+=SingleFaceSubdivision(s0,v0,pp,v2,ps,fp,randSample);
                 faceSampleNum+=SingleFaceSubdivision(s1,pp,v1,v2,ps,fp,randSample);
                 break;
     case 1 :    pp =  v1*w0 + v2*w1;
                 faceSampleNum+=SingleFaceSubdivision(s0,v0,v1,pp,ps,fp,randSample);
                 faceSampleNum+=SingleFaceSubdivision(s1,v0,pp,v2,ps,fp,randSample);
                 break;
     case 2 :    pp = v0*w0 + v2*w1;
                 faceSampleNum+=SingleFaceSubdivision(s0,v0,v1,pp,ps,fp,randSample);
                 faceSampleNum+=SingleFaceSubdivision(s1,pp,v1,v2,ps,fp,randSample);
                 break;
    }
		return faceSampleNum;
}


/// Compute a sampling of the surface where the points are regularly scattered over the face surface using a recursive longest-edge subdivision rule.
static void FaceSubdivision(MetroMesh & m, VertexSampler &ps,int sampleNum, bool randSample)
{
	
	ScalarType area = Stat<MetroMesh>::ComputeMeshArea(m);
	ScalarType samplePerAreaUnit = sampleNum/area;
	//qDebug("samplePerAreaUnit %f",samplePerAreaUnit);
	std::vector<FacePointer> faceVec;
	FillAndShuffleFacePointerVector(m,faceVec);

	double  floatSampleNum = 0.0;
	int faceSampleNum;
    // Subdivision sampling.
	typename std::vector<FacePointer>::iterator fi;
    for(fi=faceVec.begin(); fi!=faceVec.end(); fi++)
    {
        // compute # samples in the current face.
        floatSampleNum += 0.5*DoubleArea(**fi) * samplePerAreaUnit;
        faceSampleNum          = (int) floatSampleNum;
        if(faceSampleNum>0)
            faceSampleNum = SingleFaceSubdivision(faceSampleNum,(**fi).V(0)->cP(), (**fi).V(1)->cP(), (**fi).V(2)->cP(),ps,*fi,randSample);
        floatSampleNum -= (double) faceSampleNum;
    }
}


// Similar Triangles sampling.
// Skip vertex and edges
// Sample per edges includes vertexes, so here we should expect  n_samples_per_edge >=4 

static int SingleFaceSimilar(FacePointer fp, VertexSampler &ps, int n_samples_per_edge)
{
		int n_samples=0;
    int         i, j;
    float segmentNum=n_samples_per_edge -1 ;
		float segmentLen = 1.0/segmentNum;
		// face sampling.
    for(i=1; i < n_samples_per_edge-1; i++)
        for(j=1; j < n_samples_per_edge-1-i; j++)
        {
            //AddSample( v0 + (V1*(double)i + V2*(double)j) );
						CoordType sampleBary(i*segmentLen,j*segmentLen, 1.0 - (i*segmentLen+j*segmentLen) ) ;
            n_samples++;
						ps.AddFace(*fp,sampleBary);
        }
	return n_samples;
}
static int SingleFaceSimilarDual(FacePointer fp, VertexSampler &ps, int n_samples_per_edge, bool randomFlag)
{
		int n_samples=0;
    float         i, j;
    float segmentNum=n_samples_per_edge -1 ;
		float segmentLen = 1.0/segmentNum;
		// face sampling.
    for(i=0; i < n_samples_per_edge-1; i++)
        for(j=0; j < n_samples_per_edge-1-i; j++)
        {
            //AddSample( v0 + (V1*(double)i + V2*(double)j) );
						CoordType V0((i+0)*segmentLen,(j+0)*segmentLen, 1.0 - ((i+0)*segmentLen+(j+0)*segmentLen) ) ;
						CoordType V1((i+1)*segmentLen,(j+0)*segmentLen, 1.0 - ((i+1)*segmentLen+(j+0)*segmentLen) ) ;
						CoordType V2((i+0)*segmentLen,(j+1)*segmentLen, 1.0 - ((i+0)*segmentLen+(j+1)*segmentLen) ) ;
						n_samples++;
						if(randomFlag) 	{
											CoordType rb=RandomBaricentric();
											ps.AddFace(*fp, V0*rb[0]+V1*rb[1]+V2*rb[2]);
							} else  ps.AddFace(*fp,(V0+V1+V2)/3.0);

				if( j < n_samples_per_edge-i-2 )
							{
									CoordType V3((i+1)*segmentLen,(j+1)*segmentLen, 1.0 - ((i+1)*segmentLen+(j+1)*segmentLen) ) ;
									n_samples++;
									if(randomFlag) 	{
														CoordType rb=RandomBaricentric();
														ps.AddFace(*fp, V3*rb[0]+V1*rb[1]+V2*rb[2]);
										} else  ps.AddFace(*fp,(V3+V1+V2)/3.0);								
							}
        }
	return n_samples;
}

// Similar sampling 
// Each triangle is subdivided into similar triangles following a generalization of the classical 1-to-4 splitting rule of triangles. 
// According to the level of subdivision <k> you get 1, 4 , 9, 16 , <k^2> triangles. 
// Depending on the kind of the sampling strategies we can have two different approach to choosing the sample points. 
// 1) you have already sampled both edges and vertices
// 2) you are not going to take samples on edges and vertices. 
// 
// In the first case you have to consider only internal vertices of the subdivided triangles (to avoid multiple sampling of edges and vertices).
// Therefore the number of internal points is ((k-3)*(k-2))/2. where k is the number of points on an edge (vertex included)
// E.g. for k=4 you get 3 segments on each edges and the original triangle is subdivided 
// into 9 smaller triangles and you get (1*2)/2 == 1 only a single internal point.
// So if you want N samples in a triangle you have to solve  k^2 -5k +6 - 2N = 0 
// from which you get:
//
//      5 + sqrt( 1 + 8N ) 
// k = -------------------  
//             2
//
// In the second case if you are not interested to skip the sampling on edges and vertices you have to consider as sample number the number of triangles. 
// So if you want N samples in a triangle, the number <k> of points on  an edge (vertex included) should be simply:
//      k = 1 + sqrt(N)  
// examples: 
// N = 4 -> k = 3
// N = 9 -> k = 4 



//template <class MetroMesh>
//void Sampling<MetroMesh>::SimilarFaceSampling()
static void FaceSimilar(MetroMesh & m, VertexSampler &ps,int sampleNum, bool dualFlag, bool randomFlag)
{	
		ScalarType area = Stat<MetroMesh>::ComputeMeshArea(m);
		ScalarType samplePerAreaUnit = sampleNum/area;

		// Similar Triangles sampling.
    int n_samples_per_edge;
    double  n_samples_decimal = 0.0;
    FaceIterator fi;

    for(fi=m.face.begin(); fi != m.face.end(); fi++)
    {
        // compute # samples in the current face.
        n_samples_decimal += 0.5*DoubleArea(*fi) * samplePerAreaUnit;
        int n_samples          = (int) n_samples_decimal;
        if(n_samples>0)
        {
            // face sampling.
            if(dualFlag) 
							{	
									n_samples_per_edge = (int)((sqrt(1.0+8.0*(double)n_samples) +5.0)/2.0); // original for non dual case
									n_samples = SingleFaceSimilar(&*fi,ps, n_samples_per_edge);
							} else {	
									n_samples_per_edge = (int)(sqrt((double)n_samples) +1.0);
									n_samples = SingleFaceSimilarDual(&*fi,ps, n_samples_per_edge,randomFlag);
						}
        }
        n_samples_decimal -= (double) n_samples;
    }
}


	// Rasterization fuction
	// Take a triangle 
	// T deve essere una classe funzionale che ha l'operatore ()
	// con due parametri x,y di tipo S esempio:
	// class Foo { public void operator()(int x, int y ) { ??? } };



static void SingleFaceRaster(FaceType &f,  VertexSampler &ps, const Point2<ScalarType> & v0, const Point2<ScalarType> & v1, const Point2<ScalarType> & v2)
{
	typedef ScalarType S;
		// Calcolo bounding box
	Box2i bbox;

	if(v0[0]<v1[0]) { bbox.min[0]=int(v0[0]); bbox.max[0]=int(v1[0]); }
	else            { bbox.min[0]=int(v1[0]); bbox.max[0]=int(v0[0]); }
	if(v0[1]<v1[1]) { bbox.min[1]=int(v0[1]); bbox.max[1]=int(v1[1]); }
	else            { bbox.min[1]=int(v1[1]); bbox.max[1]=int(v0[1]); }
	     if(bbox.min[0]>int(v2[0])) bbox.min[0]=int(v2[0]);
	else if(bbox.max[0]<int(v2[0])) bbox.max[0]=int(v2[0]);
	     if(bbox.min[1]>int(v2[1])) bbox.min[1]=int(v2[1]);
	else if(bbox.max[1]<int(v2[1])) bbox.max[1]=int(v2[1]);

		// Calcolo versori degli spigoli
	Point2<S> d10 = v1 - v0;
	Point2<S> d21 = v2 - v1;
	Point2<S> d02 = v0 - v2;

		// Preparazione prodotti scalari
	S b0  = (bbox.min[0]-v0[0])*d10[1] - (bbox.min[1]-v0[1])*d10[0];
	S b1  = (bbox.min[0]-v1[0])*d21[1] - (bbox.min[1]-v1[1])*d21[0];
	S b2  = (bbox.min[0]-v2[0])*d02[1] - (bbox.min[1]-v2[1])*d02[0];
		// Preparazione degli steps
	S db0 =  d10[1];
	S db1 =  d21[1];
	S db2 =  d02[1];
		// Preparazione segni
	S dn0 = -d10[0];
	S dn1 = -d21[0];
	S dn2 = -d02[0];
		// Rasterizzazione

	double de = v0[0]*v1[1]-v0[0]*v2[1]-v1[0]*v0[1]+v1[0]*v2[1]-v2[0]*v1[1]+v2[0]*v0[1];

	for(int x=bbox.min[0];x<=bbox.max[0];++x)
	{
		bool in = false;
		S n0  = b0;
		S n1  = b1;
		S n2  = b2;
		for(int y=bbox.min[1];y<=bbox.max[1];++y)
		{
			if( (n0>=0 && n1>=0 && n2>=0) || (n0<=0 && n1<=0 && n2<=0) )
			{
				CoordType baryCoord;
				baryCoord[0] =  double(-y*v1[0]+v2[0]*y+v1[1]*x-v2[0]*v1[1]+v1[0]*v2[1]-x*v2[1])/de;
				baryCoord[1] = -double( x*v0[1]-x*v2[1]-v0[0]*y+v0[0]*v2[1]-v2[0]*v0[1]+v2[0]*y)/de;
				baryCoord[2] = 1-baryCoord[0]-baryCoord[1];

				ps.AddTextureSample(f, baryCoord, Point2i(x,y));
				in = true;
			} else if(in) break;
			n0 += dn0;
			n1 += dn1;
			n2 += dn2;
		}
		b0 += db0;
		b1 += db1;
		b2 += db2;
	}
}

// Generate a random point in volume defined by a box with uniform distribution
static CoordType RandomBox(vcg::Box3<ScalarType> box)
{
	CoordType p = box.min;
	p[0] += box.Dim()[0] * RandomDouble01();
	p[1] += box.Dim()[1] * RandomDouble01();
	p[2] += box.Dim()[2] * RandomDouble01();
	return p;
}

// generate Poisson-disk sample using a set of pre-generated samples (with the Montecarlo algorithm)
// It always return a point.
static VertexPointer getPrecomputedMontecarloSample(Point3i *cell, MontecarloSHT & samplepool)
{
	MontecarloSHTIterator cellBegin;
	MontecarloSHTIterator cellEnd;
	samplepool.Grid(*cell, cellBegin, cellEnd);
	return *cellBegin;
}

// check the radius constrain
static bool checkPoissonDisk(MetroMesh & vmesh, SampleSHT & sht, const Point3<ScalarType> & p, ScalarType radius) 
{
	// get the samples closest to the given one
	std::vector<VertexType*> closests;
	std::vector<ScalarType> distances;
	std::vector<CoordType> points;

	typedef VertTmark<MetroMesh> MarkerVert;
	MarkerVert mv;
	mv.SetMesh(&vmesh);
	typedef vcg::vertex::PointDistanceFunctor<ScalarType> VDistFunct;
	VDistFunct fn;

	Box3f bb(p-Point3f(radius,radius,radius),p+Point3f(radius,radius,radius));
	int nsamples = GridGetInBox(sht, mv, bb, closests);

	ScalarType r2 = radius*radius;
	for(int i=0; i<closests.size(); ++i)
		if(SquaredDistance(p,closests[i]->cP()) < r2)
			return false;

	return true;
}

struct PoissonDiskParam
{
	PoissonDiskParam()
	{
		adaptiveRadiusFlag = false;
		radiusVariance =1;
		MAXLEVELS = 5;
		invertQuality = false;
	}
	bool adaptiveRadiusFlag;
	float radiusVariance;
	bool invertQuality;
  int MAXLEVELS;
};

static ScalarType ComputePoissonDiskRadius(MetroMesh &origMesh, int sampleNum)
{
	ScalarType meshArea = Stat<MetroMesh>::ComputeMeshArea(origMesh);
	// Manage approximately the PointCloud Case, use the half a area of the bbox. 
	// TODO: If you had the radius a much better approximation could be done.
	if(meshArea ==0) 
		{
					meshArea = (origMesh.bbox.DimX()*origMesh.bbox.DimY() +
											origMesh.bbox.DimX()*origMesh.bbox.DimZ() +
											origMesh.bbox.DimY()*origMesh.bbox.DimZ()); 	
		}
	ScalarType diskRadius = sqrt(meshArea / (0.7 * M_PI * sampleNum)); // 0.7 is a density factor										
	return diskRadius;
}

static int ComputePoissonSampleNum(MetroMesh &origMesh, ScalarType diskRadius)
{
	ScalarType meshArea = Stat<MetroMesh>::ComputeMeshArea(origMesh);
	int sampleNum = meshArea /  (diskRadius*diskRadius *M_PI *0.7)  ; // 0.7 is a density factor
	return sampleNum;
}

static void ComputePoissonSampleRadii(MetroMesh &sampleMesh, ScalarType diskRadius, ScalarType radiusVariance, bool invert)
{
	VertexIterator vi;
	std::pair<float,float> minmax = tri::Stat<MetroMesh>::ComputePerVertexQualityMinMax( sampleMesh);
	float minRad = diskRadius / radiusVariance;
	float maxRad = diskRadius * radiusVariance;
	float deltaQ = minmax.second-minmax.first;
	float deltaRad = maxRad-minRad;
	for (vi = sampleMesh.vert.begin(); vi != sampleMesh.vert.end(); vi++)
	{
	 (*vi).Q() = minRad + deltaRad*((invert ? minmax.second - (*vi).Q() : (*vi).Q() - minmax.first )/deltaQ);
	}
}

/** Compute a Poisson-disk sampling of the surface.
 *  The radius of the disk is computed according to the estimated sampling density.
 *
 * This algorithm is an adaptation of the algorithm of White et al. :
 *
 * "Poisson Disk Point Set by Hierarchical Dart Throwing" 
 * K. B. White, D. Cline, P. K. Egbert,
 * IEEE Symposium on Interactive Ray Tracing, 2007,
 * 10-12 Sept. 2007, pp. 129-132.
 */
static void Poissondisk(MetroMesh &origMesh, VertexSampler &ps, MetroMesh &montecarloMesh, ScalarType diskRadius, const struct PoissonDiskParam pp=PoissonDiskParam())
{
	int cellusedcounter[20];          // cells used for each level
	int cellstosubdividecounter[20];  // cells to subdivide for each level
	int samplesgenerated[20];         // samples generated for each level
	int samplesaccepted[20];          // samples accepted for each level
	int verticescounter[20];          // vertices added to the spatial hash table

	for (int i = 0; i < 20; i++)
	{
		cellusedcounter[i] = 0;
		cellstosubdividecounter[i] = 0;
		samplesgenerated[i] = 0;
		samplesaccepted[i] = 0;
		verticescounter[i] = 0;
	}


	MetroMesh supportMesh;

	// spatial index of montecarlo samples - used to choose a new sample to insert
	MontecarloSHT montecarloSHT;

	// spatial hash table of the generated samples - used to check the radius constrain
	SampleSHT checkSHT;

	// initialize spatial hash table for searching
	ScalarType cellsize = diskRadius / sqrt(3.0);

	// inflating
	origMesh.bbox.Offset(cellsize);

	int sizeX = vcg::math::Max(1.0f,origMesh.bbox.DimX() / cellsize);
	int sizeY = vcg::math::Max(1.0f,origMesh.bbox.DimY() / cellsize);
	int sizeZ = vcg::math::Max(1.0f,origMesh.bbox.DimZ() / cellsize);
	Point3i gridsize(sizeX, sizeY, sizeZ);
#ifdef QT_VERSION
	qDebug("PDS: radius %f Grid:(%i %i %i) ",diskRadius,sizeX,sizeY,sizeZ);
#endif
	// initialize spatial hash to index pre-generated samples
	VertexIterator vi;
	montecarloSHT.InitEmpty(origMesh.bbox, gridsize);
	for (vi = montecarloMesh.vert.begin(); vi != montecarloMesh.vert.end(); vi++)
	{
		montecarloSHT.Add(&(*vi));
		verticescounter[0]++;
	}
#ifdef QT_VERSION
	qDebug("PDS: Completed montercarloSHT, inserted %i vertex in %i cells", montecarloMesh.vn, montecarloSHT.AllocatedCells.size());
#endif
	// initialize spatial hash table for check poisson-disk radius constrain
	checkSHT.InitEmpty(origMesh.bbox, gridsize);


	// sampling algorithm
	// ------------------
	//
	// - generate millions of samples using montecarlo algorithm
	// - extract a cell (C) from the active cell list (with probability proportional to cell's volume)
	// - generate a sample inside C by choosing one of the contained pre-generated samples
	//   - if the sample violates the radius constrain discard it, and add the cell to the cells-to-subdivide list
	// - iterate until the active cell list is empty or a pre-defined number of subdivisions is reached
	//

	std::vector<Point3i *> activeCells;
	std::vector<VertexType *> nextPoints;
	typename std::vector<VertexType *>::iterator nextPointsIt;

	typename std::vector<Point3i>::iterator it;
	Point3i *currentCell;
	vcg::Box3<ScalarType> currentBox;
	int level = 0;
	
  // if we are doing variable density sampling we have to prepare the random samples quality with the correct expected radii.
	if(pp.adaptiveRadiusFlag) 
			ComputePoissonSampleRadii(montecarloMesh, diskRadius, pp.radiusVariance, pp.invertQuality);
	
	do
	{
		// extract a cell (C) from the active cell list (with probability proportional to cell's volume)
		///////////////////////////////////////////////////////////////////////////////////////////////////

		supportMesh.vert.reserve(montecarloMesh.vn);

		// create active cell list
		for (it = montecarloSHT.AllocatedCells.begin(); it != montecarloSHT.AllocatedCells.end(); it++)
		{
			activeCells.push_back(&(*it));
		}

		int ncell = static_cast<int>(activeCells.size());
		cellusedcounter[level] = ncell;
		// shuffle active cells
//		int index,index2;
//		Point3i *temp;
//		for (int i = 0; i < ncell/2; i++)
//		{
//			index = RandomInt(ncell);
//			index2 =  RandomInt(ncell);
//			temp = activeCells[index];
//			activeCells[index] = activeCells[index2];
//			activeCells[index2] = temp;
//		}
		
		unsigned int (*p_myrandom)(unsigned int) = RandomInt;
		std::random_shuffle(activeCells.begin(),activeCells.end(), p_myrandom);


		// generate a sample inside C by choosing one of the contained pre-generated samples
		//////////////////////////////////////////////////////////////////////////////////////////

		for (int i = 0; i < ncell; i++)
		{
			currentCell = activeCells[i];
			//vcg::Point3<ScalarType > s; // current sample

			// generate a sample chosen from the pre-generated one
			VertexPointer sp = getPrecomputedMontecarloSample(currentCell, montecarloSHT);
			samplesgenerated[level]++;
		
			// vr spans between 3.0*r and r / 4.0 according to vertex quality
			ScalarType sampleRadius = diskRadius;
			if(pp.adaptiveRadiusFlag)  sampleRadius = sp->Q();

			if (checkPoissonDisk(*ps.m, checkSHT, sp->cP(), sampleRadius))
			{
				// add sample
				tri::Allocator<MetroMesh>::AddVertices(supportMesh,1);
				supportMesh.vert.back().P() = sp->P();
				supportMesh.vert.back().Q() = sampleRadius;
				//ps.AddVert(supportMesh.vert.back()); Small change, we should call the sampler class with the input mesh. 
				ps.AddVert(*sp);

				// add to control spatial index
				checkSHT.Add(&supportMesh.vert.back());

				samplesaccepted[level]++;
			}
			else
			{
				// subdivide this cell
				///////////////////////////////////////////////////////////////////////
				// pre-generated samples for the next level of subdivision
				MontecarloSHTIterator ptBegin, ptEnd, ptIt;
				montecarloSHT.Grid(*currentCell, ptBegin, ptEnd);

				for (ptIt = ptBegin; ptIt != ptEnd; ++ptIt)
				{
					nextPoints.push_back(*ptIt);
				}

				cellstosubdividecounter[level]++;
			}
		}

		activeCells.clear();

		// proceed to the next level of subdivision
		///////////////////////////////////////////////////////////////////////////

		// cleaning spatial index data structures
		montecarloSHT.Clear();

		// increase grid resolution
		gridsize[0] *= 2;
		gridsize[1] *= 2;
		gridsize[2] *= 2;

		montecarloSHT.InitEmpty(origMesh.bbox, gridsize);

		for (nextPointsIt = nextPoints.begin(); nextPointsIt != nextPoints.end(); nextPointsIt++)
		{
			montecarloSHT.Add(*nextPointsIt);
			verticescounter[level+1]++;
		}

		nextPoints.clear();
#ifdef QT_VERSION
		qDebug("PDS: Completed Level %i, added %i samples",level,samplesaccepted[level]);	
#endif
		level++;

	} while(level < pp.MAXLEVELS);

#ifdef DEBUG_DUMP_STAT
	// write some statistics
	QFile outfile("C:/temp/poissondisk_statistics.txt");
	if (outfile.open(QFile::WriteOnly | QFile::Truncate)) 
	{
		QTextStream out(&outfile);

		for (int k = 0; k < pp.MAXLEVELS; k++)
			out << "Cells used for level " << k << ": " << cellusedcounter[k] << endl;

		for (int k = 0; k < pp.MAXLEVELS; k++)
			out << "Cells to subdivide for level " << k << ": " << cellstosubdividecounter[k] << endl;

		for (int k = 0; k < pp.MAXLEVELS; k++)
			out << "Vertices counter for level " << k << ": " << verticescounter[k] << endl;

		for (int k = 0; k < pp.MAXLEVELS; k++)
			out << "Samples generated for level " << k << ": " << samplesgenerated[k] << endl;

		for (int k = 0; k < pp.MAXLEVELS; k++)
			out << "Samples accepted for level " << k << ": " << samplesaccepted[k] << endl;
	}

	outfile.close();
#endif
}

//template <class MetroMesh>
//void Sampling<MetroMesh>::SimilarFaceSampling()
static void Texture(MetroMesh & m, VertexSampler &ps, int textureWidth, int textureHeight)
{
		FaceIterator fi;

		printf("Similar Triangles face sampling\n");
		for(fi=m.face.begin(); fi != m.face.end(); fi++)
				{
					Point2f ti[3];
					for(int i=0;i<3;++i)
							ti[i]=Point2f((*fi).WT(i).U() * textureWidth, (*fi).WT(i).V() * textureHeight);
					
					SingleFaceRaster(*fi,  ps, ti[0],ti[1],ti[2]);
				}
}

typedef GridStaticPtr<FaceType, ScalarType > TriMeshGrid;

class RRParam
{
public:
float offset;
float minDiag;
tri::FaceTmark<MetroMesh> markerFunctor;
TriMeshGrid gM;
};

static void RegularRecursiveOffset(MetroMesh & m, std::vector<Point3f> &pvec, ScalarType offset, float minDiag)
{
	Box3<ScalarType> bb=m.bbox;
	bb.Offset(offset*2.0);
  
	RRParam rrp;

	rrp.markerFunctor.SetMesh(&m);

	rrp.gM.Set(m.face.begin(),m.face.end(),bb);
	

	rrp.offset=offset;
	rrp.minDiag=minDiag;
	SubdivideAndSample(m, pvec, bb, rrp, bb.Diag());
}

static void SubdivideAndSample(MetroMesh & m, std::vector<Point3f> &pvec, const Box3<ScalarType> bb, RRParam &rrp, float curDiag)
{
	Point3f startPt = bb.Center();
	ScalarType bound; 
	
	ScalarType dist; 
	// Compute mesh point nearest to bb center	
	FaceType   *nearestF=0;
	float dist_upper_bound = curDiag+rrp.offset;
	Point3f closestPt;
	vcg::face::PointDistanceBaseFunctor<ScalarType> PDistFunct;
	dist=dist_upper_bound;
	nearestF =  rrp.gM.GetClosest(PDistFunct,rrp.markerFunctor,startPt,dist_upper_bound,dist,closestPt);
  curDiag /=2;
	if(dist < dist_upper_bound) 
		{
			if(curDiag/3 < rrp.minDiag) //store points only for the last level of recursion (?)
				{
					if(rrp.offset==0) 
							pvec.push_back(closestPt);
					else 
						{
							if(dist>rrp.offset) // points below the offset threshold cannot be displaced at the right offset distance, we can only make points nearer.
							{
								Point3f delta = startPt-closestPt;
								pvec.push_back(closestPt+delta*(rrp.offset/dist));
							}
						}
				}
			if(curDiag < rrp.minDiag) return;
			Point3f hs = (bb.max-bb.min)/2;
			for(int i=0;i<2;i++)
				for(int j=0;j<2;j++)
					for(int k=0;k<2;k++)
						SubdivideAndSample(m,pvec,
																			Box3f(Point3f( bb.min[0]+i*hs[0], bb.min[1]+j*hs[1], bb.min[2]+k*hs[2]),
																						Point3f(startPt[0]+i*hs[0],startPt[1]+j*hs[1],startPt[2]+k*hs[2])),rrp,curDiag);
																			
		}
} 
}; // end class


} // end namespace tri
} // end namespace vcg

#endif
	