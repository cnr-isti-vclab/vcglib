#include <vcg/complex/complex.h>

#include <wrap/io_trimesh/import.h>
#include <wrap/io_trimesh/export_ply.h>

#include <vcg/complex/algorithms/geodesic_heat.h>
#include <vcg/complex/algorithms/update/color.h>

#include <iostream>


using namespace vcg;
using namespace std;

class MyEdge;
class MyFace;
class MyVertex;
struct MyUsedTypes : public UsedTypes<	Use<MyVertex>   ::AsVertexType,
                                        Use<MyEdge>     ::AsEdgeType,
                                        Use<MyFace>     ::AsFaceType>{};

class MyVertex  : public Vertex<MyUsedTypes, vertex::Coord3f, vertex::VFAdj, vertex::Color4b, vertex::Qualityf, vertex::BitFlags>{};
class MyFace    : public Face< MyUsedTypes, face::VFAdj, face::FFAdj, face::VertexRef, face::Normal3f, face::BitFlags, vertex::Qualityf> {};
class MyEdge    : public Edge<MyUsedTypes>{};
class MyMesh    : public tri::TriMesh< vector<MyVertex>, vector<MyFace> , vector<MyEdge>  > {};

int main( int argc, char **argv )
{
    if(argc<2)
    {
        printf("Usage trimesh_geodesic <meshfilename> <backEulerStepfloat>\n");
        return -1;
    }
    float param_m = 1.0; 
    if (argc>2){
        param_m = atof(argv[2]);
    }

    MyMesh m;
    if(tri::io::Importer<MyMesh>::Open(m, argv[1])!=0)
    {
        printf("Error reading file  %s\n",argv[1]);
        exit(0);
    }

    Point3f c = m.bbox.Center();
    MyVertex* closest = &*m.vert.begin();
    float minDist = Distance(closest->P(),c);
    for(MyMesh::VertexIterator vi=m.vert.begin();vi!=m.vert.end(); ++vi){
        if(Distance(vi->P(),c)<minDist){
            minDist = Distance(vi->P(),c);
            closest = &*vi;
        }
    }
    vector<MyVertex*> seedVec;
    seedVec.push_back(closest);
    tri::Clean<MyMesh>::RemoveUnreferencedVertex(m);

    bool success;
    success = tri::GeodesicHeat<MyMesh>::Compute(m, seedVec, param_m);
    if (!success){
        printf("computation of HeatGeodesic has failed! %d", success);
        exit(0);
    }
    pair<float,float> minmax = tri::Stat<MyMesh>::ComputePerVertexQualityMinMax(m);
    tri::UpdateColor<MyMesh>::PerVertexQualityRamp(m);
    printf("min %f max %f\n",minmax.first, minmax.second);
    tri::io::ExporterPLY<MyMesh>::Save(m,"base.ply",tri::io::Mask::IOM_VERTCOLOR | tri::io::Mask::IOM_VERTQUALITY);

    int t0=clock();
    tri::GeodesicHeat<MyMesh>::GeodesicHeatCache cache = tri::GeodesicHeat<MyMesh>::BuildCache(m, 1.0);
    tri::GeodesicHeat<MyMesh>::ComputeFromCache(m, seedVec, cache);
    int t1=clock();
    tri::GeodesicHeat<MyMesh>::ComputeFromCache(m, seedVec, cache);
    int t2=clock();
    printf("Non-Cached Time: %6.3f\n",float(t1-t0)/CLOCKS_PER_SEC);
    printf("Cached Time    : %6.3f\n",float(t2-t1)/CLOCKS_PER_SEC);
    tri::io::ExporterPLY<MyMesh>::Save(m,"base_m1.ply",tri::io::Mask::IOM_VERTCOLOR | tri::io::Mask::IOM_VERTQUALITY);

    return 0;
}
