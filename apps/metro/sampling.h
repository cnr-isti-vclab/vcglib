/*****************************************************************************
 * VCGLib                                                                    *
 *																																					 *
 * Visual Computing Group                                                o>  *
 * IEI Institute, CNUCE Institute, CNR Pisa                             <|   *
 *                                                                      / \  *
 * Copyright(C) 1999 by Paolo Cignoni, Claudio Rocchini                      *
 * All rights reserved.                                                      *
 *																																					 *
 * Permission  to use, copy, modify, distribute  and sell this  software and *
 * its documentation for any purpose is hereby granted without fee, provided *
 * that  the above copyright notice appear  in all copies and that both that *
 * copyright   notice  and  this  permission  notice  appear  in  supporting *
 * documentation. the author makes  no representations about the suitability *
 * of this software for any purpose. It is provided  "as is" without express *
 * or implied warranty.                                                      *
 *					                                                         *
 *****************************************************************************/
/****************************************************************************
  History
	2003 Dic 17 modifiche per conversione alla versione template	(gano)
	2004 Jan 19 qualche ->P() in ->cP() (Snarf)
****************************************************************************/

// -----------------------------------------------------------------------------------------------
#ifndef _SAMPLING_H
#define _SAMPLING_H
// -----------------------------------------------------------------------------------------------

// standard libraries.
#include <time.h>

// VCG library.
//#include <vcg/tools/xml/xml.h>
#include "min_dist_point.h"
//#include <vcg/tools/Align/Hist.h>
#include <vcg/space/box3.h>
#include <vcg/space/color4.h>
#include <vcg/space/index/grid_static_ptr.h>
using namespace vcg;
// -----------------------------------------------------------------------------------------------

// flags.
#define FLAG_HIST                           0x0001
#define FLAG_VERTEX_SAMPLING                0x0002
#define FLAG_EDGE_SAMPLING                  0x0004
#define FLAG_FACE_SAMPLING                  0x0008
#define FLAG_MONTECARLO_SAMPLING            0x0010
#define FLAG_SUBDIVISION_SAMPLING           0x0020
#define FLAG_SIMILAR_TRIANGLES_SAMPLING     0x0040
#define FLAG_SAVE_ERROR_DISPLACEMENT        0x0080 
#define FLAG_SAVE_ERROR_AS_COLOUR           0x0100

// global constants
#define NO_SAMPLES_PER_FACE             10
#define N_SAMPLES_EDGE_TO_FACE_RATIO    0.1
#define BBOX_FACTOR                     0.1
#define INFLATE_PERCENTAGE			    0.02
#define MIN_SIZE					    125		/* 125 = 5^3 */
#define N_HIST_BINS                     256
#define PRINT_EVERY_N_ELEMENTS          1000
#define FILE_EXT_SMF                    "smf"
#define FILE_EXT_PLY                    "ply"

// -----------------------------------------------------------------------------------------------
template <class MetroMesh>
class Sampling 
{
private:
    typedef GridStaticPtr< typename MetroMesh::FaceContainer > MetroMeshGrid;
		typedef Point3<typename MetroMesh::ScalarType> Point3x;

    typedef typename MetroMesh::CoordType CoordType;
    typedef typename MetroMesh::ScalarType ScalarType;
    typedef typename MetroMesh::VertexPointer  VertexPointer;
    typedef typename MetroMesh::VertexIterator VertexIterator;
    typedef typename MetroMesh::FaceIterator   FaceIterator;
    typedef typename MetroMesh::FaceType   FaceType;



    // data structures
    MetroMesh       &S1; 
    MetroMesh       &S2;
    MetroMeshGrid   gS2;


    // parameters
    double          dist_upper_bound; 
    double					n_samples_per_area_unit;
    unsigned long   n_samples_target;
    int             Flags;
    
    // results
//    Hist            hist;
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
    void            Hausdorff();
    double          GetArea()                   {return area_S1;}
    double          GetDistMax()                {return max_dist;}
    double          GetDistMean()               {return mean_dist;}
    double          GetDistRMS()                {return RMS_dist;}
    double          GetDistVolume()             {return volume;}
    unsigned long   GetNSamples()               {return n_total_samples;}
    unsigned long   GetNAreaSamples()           {return n_total_area_samples;}
    unsigned long   GetNEdgeSamples()           {return n_total_edge_samples;}
    unsigned long   GetNVertexSamples()         {return n_total_vertex_samples;}
    double					GetNSamplesPerAreaUnit()    {return n_samples_per_area_unit;}
    unsigned long   GetNSamplesTarget()         {return n_samples_target;}
//    Hist            &GetHist()                  {return hist;}
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
        area += face->Area();

    return area;
}

template <class MetroMesh>
float Sampling<MetroMesh>::AddSample(const Point3x &p)
{
    FaceType   *f=0;
    Point3x             normf, bestq, ip;
		ScalarType              dist;

    dist = dist_upper_bound;

    // compute distance between p_i and the mesh S2
    MinDistPoint(S2, p, gS2, dist, normf, bestq, f, ip);

    // update distance measures
    if(dist == dist_upper_bound)
        return -1.0;
    
    if(dist > max_dist)
        max_dist = dist;        // L_inf
    mean_dist += dist;	        // L_1
    RMS_dist  += dist*dist;     // L_2
    n_total_samples++;

    //if(Flags & FLAG_HIST)
    //    hist.Add((float)fabs(dist));

    return (float)dist;    
}


// -----------------------------------------------------------------------------------------------
// --- Vertex Sampling ---------------------------------------------------------------------------

template <class MetroMesh>
void Sampling<MetroMesh>::VertexSampling()
{
    // Vertex sampling.
    int   cnt = 0;
    float error;
    
    printf("Vertex sampling\n");
    VertexIterator vi;
    for(vi=S1.vert.begin();vi!=S1.vert.end();++vi)
    {
        error = AddSample((*vi).cP());
        n_total_vertex_samples++;

        // save vertex quality
        if(Flags & (FLAG_SAVE_ERROR_DISPLACEMENT | FLAG_SAVE_ERROR_AS_COLOUR))
            (*vi).Q() = error;

/*
        if(Flags & FLAG_SAVE_ERROR_AS_COLOUR)
        {
            ColorUB  col = ColorUB(ColorUB::White);

            if(error < dist_upper_bound)
                // colour mapped distance
                col.ColorRamp(0, dist_upper_bound, dist_upper_bound-error);
            //else
                // no matching mesh patches -> white

            (*vi).C() = col;
        }
*/

        // print progress information
        if(!(++cnt % PRINT_EVERY_N_ELEMENTS))
            printf("Sampling vertices %d%%\r", (100 * cnt/S1.vn));
    }
    printf("                       \r");
}


// -----------------------------------------------------------------------------------------------
// --- Edge Sampling -----------------------------------------------------------------------------

template <class MetroMesh>
inline void Sampling<MetroMesh>::SampleEdge(const Point3x & v0, const Point3x & v1, int n_samples_per_edge)
{
    // uniform sampling of the segment v0v1.
    Point3x     e((v1-v0)/(double)(n_samples_per_edge+1));
    int         i;
    
    for(i=1; i <= n_samples_per_edge; i++)
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
    if(Flags & FLAG_FACE_SAMPLING)
        n_samples_per_length_unit = sqrt((double)n_samples_per_area_unit);
    else
        n_samples_per_length_unit = n_samples_per_area_unit;
    for(ei=Edges.begin(); ei!=Edges.end(); ++ei)
    {
        n_samples_decimal += Distance((*ei).first->cP(),(*ei).second->cP()) * n_samples_per_length_unit;
        n_samples          = (int) n_samples_decimal;
        SampleEdge((*ei).first->cP(), (*ei).second->cP(), (int) n_samples);
        n_samples_decimal -= (double) n_samples;

        // print progress information
        if(!(++cnt % PRINT_EVERY_N_ELEMENTS))
            printf("Sampling edge %d%%\r", (100 * cnt/Edges.size()));
    }
    printf("                     \r");
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

    srand(clock());
 //   printf("Montecarlo face sampling\n");
    for(fi=S1.face.begin(); fi != S1.face.end(); fi++)
		if(!(*fi).IsD())
    {
        // compute # samples in the current face.
        n_samples_decimal += fi->Area() * n_samples_per_area_unit;
        n_samples          = (int) n_samples_decimal;

        // for every sample p_i in T...
        for(int i=0; i < n_samples; i++)
            AddRandomSample(fi);

        n_samples_decimal -= (double) n_samples;

        // print progress information
//        if(!(++cnt % PRINT_EVERY_N_ELEMENTS))
 //           printf("Sampling face %d%%\r", (100 * cnt/S1.fn));
    }
 //   printf("                     \r");
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
        n_samples++;
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
        n_samples_decimal += fi->Area() * n_samples_per_area_unit;
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
        if(!(++cnt % PRINT_EVERY_N_ELEMENTS))
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
        n_samples_decimal += fi->Area() * n_samples_per_area_unit;
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
        if(!(++cnt % PRINT_EVERY_N_ELEMENTS))
            printf("Sampling face %d%%\r", (100 * cnt/S1.fn));
    }
    printf("                     \r");
}


// -----------------------------------------------------------------------------------------------
// --- Distance ----------------------------------------------------------------------------------

template <class MetroMesh>
void Sampling<MetroMesh>::Hausdorff()
{
		Box3< ScalarType> bbox;

    // set grid meshes.
    gS2.SetBBox(S2.bbox);
	if(S2.face.size() < MIN_SIZE)
		gS2.Set(S2.face, MIN_SIZE);
    else
		gS2.Set(S2.face);

    // set bounding box
    bbox = S2.bbox;
    dist_upper_bound = /*BBOX_FACTOR * */bbox.Diag();
    //if(Flags & FLAG_HIST)
    //	hist.SetRange(0.0, dist_upper_bound, N_HIST_BINS);

    // initialize sampling statistics.
    n_total_area_samples = n_total_edge_samples = n_total_vertex_samples = n_total_samples = n_samples = 0;
		max_dist             = -HUGE_VAL;
		mean_dist = RMS_dist = 0;

    // Vertex sampling.
    if(Flags & FLAG_VERTEX_SAMPLING)
        VertexSampling();
    // Edge sammpling.
    n_samples_target -= (int) n_total_samples;
    if(n_samples_target > 0)
    {
        n_samples_per_area_unit  = n_samples_target / area_S1;
        if(Flags & FLAG_EDGE_SAMPLING)
        {
            EdgeSampling();
            n_samples_target -= (int) n_total_samples;
        }
        // Face sampling.
        if((Flags & FLAG_FACE_SAMPLING) && (n_samples_target > 0))
        {
            n_samples_per_area_unit  = n_samples_target / area_S1;
            if(Flags & FLAG_MONTECARLO_SAMPLING)        MontecarloFaceSampling();
            if(Flags & FLAG_SUBDIVISION_SAMPLING)       SubdivFaceSampling();
            if(Flags & FLAG_SIMILAR_TRIANGLES_SAMPLING) SimilarFaceSampling();
        }
    }

    // compute vertex colour
    if(Flags & FLAG_SAVE_ERROR_AS_COLOUR)
    {
        VertexIterator vi;
        float   error;
        int     cnt = 0;
        for(vi=S1.vert.begin();vi!=S1.vert.end();++vi)
        {
            Color4b  col = Color4b(Color4b::White);
            error = (*vi).Q();
       
            if(error < dist_upper_bound)
                // colour mapped distance
                col.ColorRamp(0, (float)max_dist, (float)max_dist-error);
            //else
            // no matching mesh patches -> white
        
            (*vi).C() = col;

            // print progress information
            if(!(++cnt % PRINT_EVERY_N_ELEMENTS))
                printf("Computing vertex colour %d%%\r", (100 * cnt/S1.vn));
        }
        printf("                       \r");
    }

    // compute statistics
    n_samples_per_area_unit = (double) n_total_samples / area_S1;
    volume     = mean_dist / n_samples_per_area_unit / 2.0;
    mean_dist /= n_total_samples;
    RMS_dist   = sqrt(RMS_dist / n_total_samples);
}
// -----------------------------------------------------------------------------------------------
//#undef FLAG_HIST                           
//#undef FLAG_VERTEX_SAMPLING                
//#undef FLAG_EDGE_SAMPLING                  
//#undef FLAG_FACE_SAMPLING                  
//#undef FLAG_MONTECARLO_SAMPLING            
//#undef FLAG_SUBDIVISION_SAMPLING           
//#undef FLAG_SIMILAR_TRIANGLES_SAMPLING     
//#undef FLAG_SAVE_ERROR_DISPLACEMENT        
//#undef FLAG_SAVE_ERROR_AS_COLOUR           
//
//// global constants
//#undef NO_SAMPLES_PER_FACE           
//#undef N_SAMPLES_EDGE_TO_FACE_RATIO  
//#undef BBOX_FACTOR                   
//#undef INFLATE_PERCENTAGE						
//#undef MIN_SIZE											 
//#undef N_HIST_BINS                   
//#undef PRINT_EVERY_N_ELEMENTS        
//#undef FILE_EXT_SMF                  
//#undef FILE_EXT_PLY                  
//

#endif
// -----------------------------------------------------------------------------------------------
