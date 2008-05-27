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
namespace vcg
{
namespace tri
{
template <class MetroMesh, class VertexSampler>
class SurfaceSampling
{
	  typedef typename MetroMesh::CoordType				CoordType;
    typedef typename MetroMesh::ScalarType			ScalarType;
		typedef typename MetroMesh::VertexType			VertexType;
    typedef typename MetroMesh::VertexPointer		VertexPointer;
    typedef typename MetroMesh::VertexIterator	VertexIterator;
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
static void VertexUniform(MetroMesh & m, VertexSampler &ps, int sampleNum)
{
	if(sampleNum>=m.vn) 
	{
	  AllVertex(m,ps);
		return;
	} 
	
	std::vector<VertexPointer> vertVec;
	
	VertexIterator vi;
	for(vi=m.vert.begin();vi!=m.vert.end();++vi) 
				if(!(*vi).IsD())	vertVec.push_back(&*vi);

	assert(vertVec.size()=m.vn);
	
	std::random_shuffle(vertVec.begin(),vertVec.end());
	
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

		qDebug("Edges %i (unique %i) ",Edges.size(),newEnd-Edges.begin());
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

// Get the baricentric coords of a random point over a single face.

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

static void Montecarlo(MetroMesh & m, VertexSampler &ps,int sampleNum)
{
	ScalarType area = Stat<MetroMesh>::ComputeMeshArea(m);
	ScalarType samplePerAreaUnit = sampleNum/area;
	qDebug("samplePerAreaUnit %f",samplePerAreaUnit);
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

};

} // end namespace tri
} // end namespace vcg

#endif



#if 0
struct SamplingFlags{
			 enum{
						VERTEX_SAMPLING									= 0x0002,
						EDGE_SAMPLING										= 0x0004,
						FACE_SAMPLING										= 0x0008,
						MONTECARLO_SAMPLING							= 0x0010,
						SUBDIVISION_SAMPLING						= 0x0020,
						SIMILAR_SAMPLING			          = 0x0040,
						INCLUDE_UNREFERENCED_VERTICES		= 0x0200,
				};
	};
// -----------------------------------------------------------------------------------------------
template <class MetroMesh>
class Sampling
{
public:

		unsigned int n_samples_per_face    ;
		float n_samples_edge_to_face_ratio ;
		float bbox_factor                  ;
		float inflate_percentage			     ;
		unsigned int min_size					     ;
		int n_hist_bins                    ;
		int print_every_n_elements         ;
		int referredBit                    ;
    // parameters
    double          dist_upper_bound;
    double					n_samples_per_area_unit;
    unsigned long   n_samples_target;
    int             Flags;

    // results
    unsigned long   n_total_samples;
    unsigned long   n_total_area_samples;
    unsigned long   n_total_edge_samples;
    unsigned long   n_total_vertex_samples;
    double          max_dist;
    double          mean_dist;
    double          RMS_dist;
    double          volume;
    double          area_S1;

    // globals
    int             n_samples;

    // private methods
    inline double   ComputeMeshArea(MetroMesh & mesh);
    float           AddSample(const Point3x &p);
    inline void     AddRandomSample(FaceIterator &T);
    inline void     SampleEdge(const Point3x & v0, const Point3x & v1, int n_samples_per_edge);
    void            VertexSampling();
    void            EdgeSampling();
    void            FaceSubdiv(const Point3x & v0, const Point3x &v1, const Point3x & v2, int maxdepth);
    void            SimilarTriangles(const Point3x &v0, const Point3x &v1, const Point3x &v2, int n_samples_per_edge);
    void            MontecarloFaceSampling();
    void            SubdivFaceSampling();
    void            SimilarFaceSampling();

public :
    // public methods
    Sampling(MetroMesh &_s1, MetroMesh &_s2);
		~Sampling();

    double          GetDistVolume()             {return volume;}
    unsigned long   GetNSamples()               {return n_total_samples;}
    unsigned long   GetNAreaSamples()           {return n_total_area_samples;}
    unsigned long   GetNEdgeSamples()           {return n_total_edge_samples;}
    unsigned long   GetNVertexSamples()         {return n_total_vertex_samples;}
    double					GetNSamplesPerAreaUnit()    {return n_samples_per_area_unit;}
    unsigned long   GetNSamplesTarget()         {return n_samples_target;}
    
		void            SetFlags(int flags)         {Flags = flags;}
    void            ClearFlag(int flag)         {Flags &= (flag ^ -1);}
    void            SetParam(double _n_samp)    {n_samples_target = _n_samp;}
    void            SetSamplesTarget(unsigned long _n_samp);
    void            SetSamplesPerAreaUnit(double _n_samp);
};

// -----------------------------------------------------------------------------------------------

// constructor
template <class MetroMesh>
Sampling<MetroMesh>::Sampling(MetroMesh &_s1, MetroMesh &_s2):S1(_s1),S2(_s2)
{
    Flags = 0;
    area_S1 = ComputeMeshArea(_s1);
		// set default numbers
		n_samples_per_face             =	10;
		n_samples_edge_to_face_ratio   = 0.1f;
		bbox_factor                    = 0.1f;
		inflate_percentage			       = 0.02f;
		min_size					             = 125;		/* 125 = 5^3 */
		n_hist_bins                    = 256;
		print_every_n_elements         = S1.fn/100;

		if(print_every_n_elements <= 1)
		  print_every_n_elements = 2;

			referredBit = VertexType::NewBitFlag();
			// store the unreferred vertices
			FaceIterator fi; VertexIterator vi; int i;
			for(fi = _s1.face.begin(); fi!= _s1.face.end(); ++fi)
				for(i=0;i<3;++i) (*fi).V(i)->SetUserBit(referredBit);
}

template <class MetroMesh>
Sampling<MetroMesh>::~Sampling()
{
	VertexType::DeleteBitFlag(referredBit);
}


// set sampling parameters
template <class MetroMesh>
void Sampling<MetroMesh>::SetSamplesTarget(unsigned long _n_samp)
{
    n_samples_target        = _n_samp;
    n_samples_per_area_unit =  n_samples_target / (double)area_S1;
}

template <class MetroMesh>
void Sampling<MetroMesh>::SetSamplesPerAreaUnit(double _n_samp)
{
    n_samples_per_area_unit = _n_samp;
    n_samples_target        = (unsigned long)((double) n_samples_per_area_unit * area_S1);
}


// auxiliary functions
template <class MetroMesh>
inline double Sampling<MetroMesh>::ComputeMeshArea(MetroMesh & mesh)
{
    FaceIterator    face;
    double                  area = 0.0;

    for(face=mesh.face.begin(); face != mesh.face.end(); face++)
			if(!(*face).IsD())
        area += DoubleArea(*face);

    return area/2.0;
}

template <class MetroMesh>
float Sampling<MetroMesh>::AddSample(const Point3x &p )
{
	SampleVec.push_back(p);
}


// -----------------------------------------------------------------------------------------------
// --- Vertex Sampling ---------------------------------------------------------------------------

template <class MetroMesh>
void Sampling<MetroMesh>::VertexSampling()
{
    VertexIterator vi;
    for(vi=S1.vert.begin();vi!=S1.vert.end();++vi)
			if(  (*vi).IsUserBit(referredBit) || // it is referred
					((Flags&SamplingFlags::INCLUDE_UNREFERENCED_VERTICES) != 0) ) //include also unreferred
						{
								AddSample((*vi).cP());
						}
}


// -----------------------------------------------------------------------------------------------
// --- Edge Sampling -----------------------------------------------------------------------------

template <class MetroMesh>
inline void Sampling<MetroMesh>::SampleEdge(const Point3x & v0, const Point3x & v1, int n_samples_per_edge)
{
    // uniform sampling of the segment v0v1.
    Point3x     e((v1-v0)/(double)(n_samples_per_edge+1));
    for(int i=1; i <= n_samples_per_edge; i++)
    {
        AddSample(v0 + e*i);
        n_total_edge_samples++;
    }
}


template <class MetroMesh>
void Sampling<MetroMesh>::EdgeSampling()
{
    // Edge sampling.
		typedef std::pair<VertexPointer, VertexPointer> pvv;
		std::vector< pvv > Edges;

    printf("Edge sampling\n");

    // compute edge list.
    FaceIterator fi;
    for(fi=S1.face.begin(); fi != S1.face.end(); fi++)
        for(int i=0; i<3; ++i)
        {
            Edges.push_back(make_pair((*fi).V0(i),(*fi).V1(i)));
            if(Edges.back().first > Edges.back().second)
                swap(Edges.back().first, Edges.back().second);
        }
    sort(Edges.begin(), Edges.end());
		typename std::vector< pvv>::iterator edgeend = unique(Edges.begin(), Edges.end());
    Edges.resize(edgeend-Edges.begin());

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
}


// -----------------------------------------------------------------------------------------------
// --- Face Sampling -----------------------------------------------------------------------------

// Montecarlo sampling.
template <class MetroMesh>
inline void Sampling<MetroMesh>::AddRandomSample(FaceIterator &T)
{
    // random sampling over the input face.
    double      rnd_1, rnd_2;

    // vertices of the face T.
    Point3x p0(T->V(0)->cP());
    Point3x p1(T->V(1)->cP());
    Point3x p2(T->V(2)->cP());
    // calculate two edges of T.
    Point3x v1(p1 - p0);
    Point3x v2(p2 - p0);

    // choose two random numbers.
    rnd_1 = (double)rand() / (double)RAND_MAX;
    rnd_2 = (double)rand() / (double)RAND_MAX;
    if(rnd_1 + rnd_2 > 1.0)
    {
        rnd_1 = 1.0 - rnd_1;
        rnd_2 = 1.0 - rnd_2;
    }

    // add a random point on the face T.
    AddSample (p0 + (v1 * rnd_1 + v2 * rnd_2));
    n_total_area_samples++;
}

template <class MetroMesh>
void Sampling<MetroMesh>::MontecarloFaceSampling()
{
    // Montecarlo sampling.
    int     cnt = 0;
    double  n_samples_decimal = 0.0;
    FaceIterator fi;

    for(fi=S1.face.begin(); fi != S1.face.end(); fi++)
		if(!(*fi).IsD())
    {
        // compute # samples in the current face.
        n_samples_decimal += 0.5*DoubleArea(*fi) * n_samples_per_area_unit;
        n_samples          = (int) n_samples_decimal;

        // for every sample p_i in T...
        for(int i=0; i < n_samples; i++)
            AddRandomSample(fi);

        n_samples_decimal -= (double) n_samples;
    }
}


// Subdivision sampling.
template <class MetroMesh>
void Sampling<MetroMesh>::FaceSubdiv(const Point3x & v0, const Point3x & v1, const Point3x & v2, int maxdepth)
{
    // recursive face subdivision.
    if(maxdepth == 0)
    {
        // ground case.
        AddSample((v0+v1+v2)/3.0f);
        n_total_area_samples++;
        return;
    }

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

    // break the input triangle along the median to the the longest edge.
    Point3x  pp;
    switch(res)
    {
     case 0 :    pp = (v0+v1)/2;
                 FaceSubdiv(v0,pp,v2,maxdepth-1);
                 FaceSubdiv(pp,v1,v2,maxdepth-1);
                 break;
     case 1 :    pp = (v1+v2)/2;
                 FaceSubdiv(v0,v1,pp,maxdepth-1);
                 FaceSubdiv(v0,pp,v2,maxdepth-1);
                 break;
     case 2 :    pp = (v2+v0)/2;
                 FaceSubdiv(v0,v1,pp,maxdepth-1);
                 FaceSubdiv(pp,v1,v2,maxdepth-1);
                 break;
    }
}

template <class MetroMesh>
void Sampling<MetroMesh>::SubdivFaceSampling()
{
    // Subdivision sampling.
    int     cnt = 0, maxdepth;
    double  n_samples_decimal = 0.0;
    typename MetroMesh::FaceIterator fi;

    printf("Subdivision face sampling\n");
    for(fi=S1.face.begin(); fi != S1.face.end(); fi++)
    {
        // compute # samples in the current face.
        n_samples_decimal += 0.5*DoubleArea(*fi) * n_samples_per_area_unit;
        n_samples          = (int) n_samples_decimal;
        if(n_samples)
        {
            // face sampling.
            maxdepth = ((int)(log((double)n_samples)/log(2.0)));
            n_samples = 0;
            FaceSubdiv((*fi).V(0)->cP(), (*fi).V(1)->cP(), (*fi).V(2)->cP(), maxdepth);
        }
        n_samples_decimal -= (double) n_samples;

        // print progress information
        if(!(++cnt % print_every_n_elements))
            printf("Sampling face %d%%\r", (100 * cnt/S1.fn));
    }
    printf("                     \r");
}


// Similar Triangles sampling.
template <class MetroMesh>
void Sampling<MetroMesh>::SimilarTriangles(const Point3x & v0, const Point3x & v1, const Point3x & v2, int n_samples_per_edge)
{
    Point3x     V1((v1-v0)/(double)(n_samples_per_edge-1));
    Point3x     V2((v2-v0)/(double)(n_samples_per_edge-1));
    int         i, j;

    // face sampling.
    for(i=1; i < n_samples_per_edge-1; i++)
        for(j=1; j < n_samples_per_edge-1-i; j++)
        {
            AddSample( v0 + (V1*(double)i + V2*(double)j) );
            n_total_area_samples++;
            n_samples++;
        }
}

template <class MetroMesh>
void Sampling<MetroMesh>::SimilarFaceSampling()
{
    // Similar Triangles sampling.
    int     cnt = 0, n_samples_per_edge;
    double  n_samples_decimal = 0.0;
    FaceIterator fi;

    printf("Similar Triangles face sampling\n");
    for(fi=S1.face.begin(); fi != S1.face.end(); fi++)
    {
        // compute # samples in the current face.
        n_samples_decimal += 0.5*DoubleArea(*fi) * n_samples_per_area_unit;
        n_samples          = (int) n_samples_decimal;
        if(n_samples)
        {
            // face sampling.
            n_samples_per_edge = (int)((sqrt(1.0+8.0*(double)n_samples) +5.0)/2.0);
            n_samples = 0;
            SimilarTriangles((*fi).V(0)->cP(), (*fi).V(1)->cP(), (*fi).V(2)->cP(), n_samples_per_edge);
        }
        n_samples_decimal -= (double) n_samples;

        // print progress information
        if(!(++cnt % print_every_n_elements))
            printf("Sampling face %d%%\r", (100 * cnt/S1.fn));
    }
    printf("                     \r");
}

}
#endif
