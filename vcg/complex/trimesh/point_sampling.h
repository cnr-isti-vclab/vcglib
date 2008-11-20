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

$Log: sampling.h,v $


The sampling Class has a set of static functions, that you can call to sample the surface of a mesh.
Each function is templated on the mesh and on a Sampler object s. 
Each function calls many time the sample object with the sampling point as parameter.
 

****************************************************************************/
#ifndef __VCGLIB_POINT_SAMPLING
#define __VCGLIB_POINT_SAMPLING
#include <vcg/complex/trimesh/stat.h>
#include <vcg/complex/trimesh/update/topology.h>
#include <vcg/space/box2.h>
namespace vcg
{
namespace tri
{

/// Trivial Sampler, an example sampler object that show the required interface used by the sampling class. 
/// Most of the sampling classes call the AddFace method with the face containing the sample and its baricentric coord.


template <class MeshType>
class TrivialSampler
{
	public:
		typedef typename MeshType::CoordType			CoordType;
		typedef typename MeshType::VertexType			VertexType;
    typedef typename MeshType::FaceType				FaceType;

	TrivialSampler(){};
	std::vector<CoordType> sampleVec;
	
	void AddVert(const VertexType &p) 
	{
		sampleVec.push_back(p.cP());
	}
	void AddFace(const FaceType &f, const CoordType &p) 
	{
		sampleVec.push_back(f.P(0)*p[0] + f.P(1)*p[1] +f.P(2)*p[2] );
	}
	
	void AddTextureSample(const FaceType &, const CoordType &, const Point2i &)
	{
		// Retrieve the colorof the sample from the face f using the barycentric coord p 
		// and write that color in a texture image at position tp[0],tp[1]
	}
}; // end class TrivialSampler

template <class MetroMesh, class VertexSampler>
class SurfaceSampling
{
	  typedef typename MetroMesh::CoordType				CoordType;
    typedef typename MetroMesh::ScalarType			ScalarType;
		typedef typename MetroMesh::VertexType			VertexType;
    typedef typename MetroMesh::VertexPointer		VertexPointer;
    typedef typename MetroMesh::VertexIterator	VertexIterator;
    typedef typename MetroMesh::FacePointer		FacePointer;
    typedef typename MetroMesh::FaceIterator		FaceIterator;
    typedef typename MetroMesh::FaceType				FaceType;
    typedef typename MetroMesh::FaceContainer		FaceContainer;
	public:

static 
void AllVertex(MetroMesh & m, VertexSampler &ps)
{
    VertexIterator vi;
    for(vi=m.vert.begin();vi!=m.vert.end();++vi) 
				if(!(*vi).IsD())
						{
								ps.AddVert(*vi);
						}
}

/// Sample the vertices in a weighted way. Each vertex has a probabiltiy of being chosen
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

/// Sample the vertices in a uniform way. Each vertex has a probabiltiy of being chosen
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
	
	std::random_shuffle(faceVec.begin(),faceVec.end());
}
static void	FillAndShuffleVertexPointerVector(MetroMesh & m, std::vector<VertexPointer> &vertVec)
{
	VertexIterator vi;
	for(vi=m.vert.begin();vi!=m.vert.end();++vi) 
				if(!(*vi).IsD())	vertVec.push_back(&*vi);

	assert((int)vertVec.size()==m.vn);
	
	std::random_shuffle(vertVec.begin(),vertVec.end());
}
/// Sample the vertices in a uniform way. Each vertex has the same probabiltiy of being chosen. 
static void VertexUniform(MetroMesh & m, VertexSampler &ps, int sampleNum)
{
	if(sampleNum>=m.vn) 
	{
	  AllVertex(m,ps);
		return;
	} 
	
	std::vector<VertexPointer> vertVec;
	FillAndShuffleVertexPointerVector(m,vertVec);
	
	for(int i =0; i< sampleNum; ++i)
		ps.AddVert(*vertVec[i]);
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
    UpdateTopology<MetroMesh>::FillEdgeVector(m,Edges);

		sort(Edges.begin(), Edges.end());							// Lo ordino per vertici

    typename std::vector< SimpleEdge>::iterator newEnd = unique(Edges.begin(), Edges.end());
		typename std::vector<SimpleEdge>::iterator   ei;

		//qDebug("Edges %i (unique %i) ",(int)Edges.size(), (int)(newEnd-Edges.begin()) );
    Edges.resize(newEnd-Edges.begin());
		for(ei=Edges.begin(); ei!=Edges.end(); ++ei)
				{
					Point3f interp(0,0,0);
					interp[ (*ei).z     ]=.5; 
					interp[((*ei).z+1)%3]=.5;
					ps.AddFace(*(*ei).f,interp);
				}
}
/*
    // sample edges.
		typename std::vector<pvv>::iterator   ei;
    double                  n_samples_per_length_unit;
    double                  n_samples_decimal = 0.0;
    int                     cnt=0;
		
    if(Flags & SamplingFlags::FACE_SAMPLING)	n_samples_per_length_unit = sqrt((double)n_samples_per_area_unit);
															    else        n_samples_per_length_unit = n_samples_per_area_unit;
																	
    for(ei=Edges.begin(); ei!=Edges.end(); ++ei)
    {
        n_samples_decimal += Distance((*ei).first->cP(),(*ei).second->cP()) * n_samples_per_length_unit;
        n_samples          = (int) n_samples_decimal;
        SampleEdge((*ei).first->cP(), (*ei).second->cP(), (int) n_samples);
        n_samples_decimal -= (double) n_samples;
    }
 */

// Generate the baricentric coords of a random point over a single face, with a uniform distribution over the triangle. 
// It uses the parallelgoram folding trick. 

static CoordType RandomBaricentric()
{
		CoordType interp;
    interp[1] = (double)rand() / (double)RAND_MAX;
    interp[2] = (double)rand() / (double)RAND_MAX;
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
	// First loop: build a sequence of consective segments proportional to the triangle areas.
	for(fi=m.face.begin(); fi != m.face.end(); fi++)
		if(!(*fi).IsD())
    {
			intervals[i+1]=std::make_pair(intervals[i].first+0.5*DoubleArea(*fi), &*fi);
			++i;
		}
	ScalarType meshArea = intervals.back().first;
	for(i=0;i<sampleNum;++i)
		{
			ScalarType val = meshArea * (double)rand() / (double)RAND_MAX;
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

/// Similar sampling. Each triangle is subdivided into similar triangles following a generalization of the classical 1-to-4 splitting rule of triangles. 
/// According to the level of subdivision <k> you get 1, 4 , 9, 16 , <k^2> triangles. 
/// Of these triangles we consider only internal vertices. (to avoid multiple sampling of edges and vertices).
/// Therefore the number of internal points is ((k-3)*(k-2))/2. where k is the number of point on an edge (vertex included)
//  e.g. for a k=4 you get (1*2)/2 == 1 e.g. a single point, etc.
/// So if you want N samples in a triangle i have to solve  k^2 -5k +6 - 2N = 0 

//      5 + sqrt( 1 + 8N ) 
// k = -------------------  
//             2



//template <class MetroMesh>
//void Sampling<MetroMesh>::SimilarFaceSampling()
static void FaceSimilar(MetroMesh & m, VertexSampler &ps,int sampleNum)
{
	
		ScalarType area = Stat<MetroMesh>::ComputeMeshArea(m);
		ScalarType samplePerAreaUnit = sampleNum/area;
		//qDebug("samplePerAreaUnit %f",samplePerAreaUnit);

		// Similar Triangles sampling.
    int n_samples_per_edge;
    double  n_samples_decimal = 0.0;
    FaceIterator fi;

    printf("Similar Triangles face sampling\n");
    for(fi=m.face.begin(); fi != m.face.end(); fi++)
    {
        // compute # samples in the current face.
        n_samples_decimal += 0.5*DoubleArea(*fi) * samplePerAreaUnit;
        int n_samples          = (int) n_samples_decimal;
        if(n_samples)
        {
            // face sampling.
            n_samples_per_edge = (int)((sqrt(1.0+8.0*(double)n_samples) +5.0)/2.0);
            //n_samples = 0;
            //SingleFaceSimilar((*fi).V(0)->cP(), (*fi).V(1)->cP(), (*fi).V(2)->cP(), n_samples_per_edge);
						n_samples = SingleFaceSimilar(&*fi,ps, n_samples_per_edge);
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

	if(v0[0]<v1[0]) { bbox.min[0]=v0[0]; bbox.max[0]=v1[0]; }
	else            { bbox.min[0]=v1[0]; bbox.max[0]=v0[0]; }
	if(v0[1]<v1[1]) { bbox.min[1]=v0[1]; bbox.max[1]=v1[1]; }
	else            { bbox.min[1]=v1[1]; bbox.max[1]=v0[1]; }
	     if(bbox.min[0]>v2[0]) bbox.min[0]=v2[0];
	else if(bbox.max[0]<v2[0]) bbox.max[0]=v2[0];
	     if(bbox.min[1]>v2[1]) bbox.min[1]=v2[1];
	else if(bbox.max[1]<v2[1]) bbox.max[1]=v2[1];

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

/// Compute a Poisson-disk sampling of the surface
/// The radius of the disk is computed according to the estimated sampling density.
static void Poissondisk(MetroMesh &m, VertexSampler &ps, int sampleNum)
{
	FaceIterator fi;

	// first of all computed r
	ScalarType meshArea = Stat<MetroMesh>::ComputeMeshArea(m);
	ScalarType r = sqrt(meshArea / (0.7 * 3.1415 * sampleNum)); // 0.7 is a density factor

	for (fi = m.face.begin(); fi != m.face.end(); fi++)
	{
		ps.AddFace(*fi, RandomBaricentric());
	}

	for (int i = 0; i < sampleNum; i++)
	{
		//...TODO...
	}
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



}; // end class


} // end namespace tri
} // end namespace vcg

#endif
