#ifndef __VCGFOREMBREE_H
#define __VCGFOREMBREE_H

#include <iostream>
#include <vcg/complex/complex.h>

#include <vcg/complex/algorithms/update/color.h>
#include <vcg/complex/algorithms/update/quality.h>
#include <wrap/callback.h>
#include <embree3/rtcore.h>
#include <vcg/math/gen_normal.h>
#include <limits>
#include <math.h>
#include <time.h>
#include <omp.h>

class MyVertex; class MyEdge; class MyFace;
struct MyUsedTypes : public vcg::UsedTypes<vcg::Use<MyVertex>   ::AsVertexType,
                                           vcg::Use<MyEdge>     ::AsEdgeType,
                                           vcg::Use<MyFace>     ::AsFaceType>{};

class MyVertex  : public vcg::Vertex< MyUsedTypes, vcg::vertex::Coord3f, vcg::vertex::Normal3f, vcg::vertex::BitFlags, vcg::vertex::VFAdj, vcg::vertex::Qualityf, vcg::vertex::Color4b>{};
class MyFace    : public vcg::Face<   MyUsedTypes, vcg::face::FFAdj ,vcg::face::VFAdj, vcg::face::Normal3f,  vcg::face::VertexRef, vcg::face::BitFlags, vcg::face::Color4b, vcg::face::Qualityf> {};
class MyEdge    : public vcg::Edge<   MyUsedTypes> {};

class MyMesh    : public vcg::tri::TriMesh< std::vector<MyVertex>, std::vector<MyFace> , std::vector<MyEdge>  > {};

using namespace vcg;
using namespace std;



namespace vcg{
    template <class MeshType>
    class EmbreeAdaptor{

        RTCDevice device = rtcNewDevice(NULL);
        RTCScene scene = rtcNewScene(device);
        RTCGeometry geometry = rtcNewGeometry(device, RTC_GEOMETRY_TYPE_TRIANGLE);
        int threads;

        public: 
         EmbreeAdaptor(MeshType &m, int nOfthreads){
                threads = nOfthreads;
                loadVCGMeshInScene(m);  
            }

        public:       
         void computeShadows(MeshType &m, Point3f rayDirection, bool showHit, bool verbose){
            auto tup = loadVCGMeshInScene(m);
            RTCScene scene = get<0>(tup);
            RTCGeometry geometry =std::get<1>(tup);
            RTCDevice device = std::get<2>(tup);
            
            RTCRayHit rayhit; 
            MyMesh hitM; 

            for(int i = 0;i<m.FN(); i++)
            {                    
                Point3f b = vcg::Barycenter(m.face[i]);
                std::vector<Point3f> unifDirVec;
                GenNormal<float>::Fibonacci(1,unifDirVec);
                Point3f dir = unifDirVec[0];

                rayhit.ray.dir_x  = dir[0]; rayhit.ray.dir_y = dir[1]; rayhit.ray.dir_z = dir[2];
                rayhit.ray.org_x  = b[0]; rayhit.ray.org_y = b[1]; rayhit.ray.org_z = b[2];

                rayhit.ray.tnear  = 0.5f;
                rayhit.ray.tfar   = std::numeric_limits<float>::infinity();
                rayhit.hit.geomID = RTC_INVALID_GEOMETRY_ID;
                
                RTCIntersectContext context;
                rtcInitIntersectContext(&context);

                rtcIntersect1(scene, &context, &rayhit);

                if (rayhit.hit.geomID != RTC_INVALID_GEOMETRY_ID) 
                {
                if(verbose)
                    std::cout << "Intersection at t = " << rayhit.ray.tfar << endl;

                if(showHit){
                Point3f p = b+dir*rayhit.ray.tfar;
                tri::Allocator<MyMesh>::AddVertex(hitM,p);
                }

                m.face[i].C() = Color4b::Black;
                } else 
                {
                if(verbose)
                    std::cout << "No Intersection" <<endl;
                m.face[i].C() = Color4b::White;       
                }           
            }
            
            rtcReleaseScene(scene);
            rtcReleaseDevice(device);
            if (showHit)
                tri::io::ExporterOFF<MyMesh>::Save(hitM,"Test2.off",tri::io::Mask::IOM_FACECOLOR);
        }

        public:       
         void loadVCGMeshInScene(MeshType &m){
            //a little mesh preprocessing before adding it to a RTCScene          
            tri::RequirePerVertexNormal(m);
            tri::UpdateNormal<MyMesh>::PerVertexNormalized(m);
            tri::UpdateNormal<MyMesh>::PerFaceNormalized(m);
            tri::UpdateBounding<MyMesh>::Box(m);
            tri::UpdateFlags<MyMesh>::FaceClearV(m);
                       
            float* vb = (float*) rtcSetNewGeometryBuffer(geometry, RTC_BUFFER_TYPE_VERTEX, 0, RTC_FORMAT_FLOAT3, 3*sizeof(float), m.VN());
            for (int i = 0;i<m.VN(); i++){
                vb[i*3]=m.vert[i].P()[0];
                vb[i*3+1]=m.vert[i].P()[1];
                vb[i*3+2]=m.vert[i].P()[2];
            }

            unsigned* ib = (unsigned*) rtcSetNewGeometryBuffer(geometry, RTC_BUFFER_TYPE_INDEX, 0, RTC_FORMAT_UINT3, 3*sizeof(unsigned), m.FN());
            for(int i = 0;i<m.FN(); i++){
                ib[i*3] = tri::Index(m,m.face[i].V(0));
                ib[i*3+1] = tri::Index(m,m.face[i].V(1));
                ib[i*3+2] = tri::Index(m,m.face[i].V(2));
            }

            rtcCommitGeometry(geometry);
            rtcAttachGeometry(scene, geometry);
            rtcReleaseGeometry(geometry);
            rtcCommitScene(scene);                     
        }

        public:
         void computeAmbientOcclusion(MeshType &inputM, int nRay){
            std::vector<Point3f> unifDirVec;
            GenNormal<float>::Fibonacci(nRay,unifDirVec);           
            computeAmbientOcclusion(inputM, nRay, unifDirVec);
        }

        public:
         void computeAmbientOcclusion(MeshType &inputM, int nRay, std::vector<Point3f> unifDirVec){
            tri::UpdateQuality<MeshType>::FaceConstant(inputM,0);
            #pragma omp parallel shared(inputM) 
            {
                #pragma omp for 
                for(int i = 0;i<inputM.FN(); i++)
                {       
                    RTCRayHit rayhit;     
                    Point3f b = vcg::Barycenter(inputM.face[i]);
                    rayhit.ray.org_x  = b[0]; rayhit.ray.org_y = b[1]; rayhit.ray.org_z = b[2];  
                    rayhit.ray.tnear  = 0.00001f;
                    
                    for(int r = 0; r<nRay; r++){
                        Point3f dir = unifDirVec.at(r);                       
                        float scalarP = inputM.face[i].N()*dir;

                        if(scalarP>0){
                            rayhit.ray.dir_x  = dir[0]; rayhit.ray.dir_y = dir[1]; rayhit.ray.dir_z = dir[2];
                            rayhit.ray.tfar   = std::numeric_limits<float>::infinity();
                            rayhit.hit.geomID = RTC_INVALID_GEOMETRY_ID;  

                            RTCIntersectContext context;
                            rtcInitIntersectContext(&context);

                            rtcIntersect1(scene, &context, &rayhit);

                            if (rayhit.hit.geomID == RTC_INVALID_GEOMETRY_ID) 
                                inputM.face[i].Q()+=scalarP;
                                        
                        }
                    }
                }
            }
            tri::UpdateColor<MeshType>::PerFaceQualityGray(inputM);
            rtcReleaseScene(scene);
            rtcReleaseDevice(device);
        }

        public:
         void computeObscurance(MeshType &inputM, int nRay, float tau){          
            std::vector<Point3f> unifDirVec;
                GenNormal<float>::Fibonacci(nRay,unifDirVec);

           computeObscurance(inputM, nRay, unifDirVec, tau);
        }

        public:
         void computeObscurance(MeshType &inputM, int nRay, std::vector<Point3f> unifDirVec, float tau){
            tri::UpdateQuality<MeshType>::FaceConstant(inputM,0);
            
            #pragma omp parallel 
            {
                #pragma omp for
                for(int i = 0;i<inputM.FN(); i++)
                {           
                    RTCRayHit rayhit; 
                    Point3f b = vcg::Barycenter(inputM.face[i]);
                    rayhit.ray.org_x  = b[0]; rayhit.ray.org_y = b[1]; rayhit.ray.org_z = b[2];  
                    rayhit.ray.tnear  = 0.00001f;
                    
                    for(int r = 0; r<nRay; r++){
                        Point3f dir = unifDirVec.at(r);                       
                        float scalarP = inputM.face[i].N()*dir;

                        if(scalarP>0){
                            rayhit.ray.dir_x  = dir[0]; rayhit.ray.dir_y = dir[1]; rayhit.ray.dir_z = dir[2];
                            rayhit.ray.tfar   = std::numeric_limits<float>::infinity();
                            rayhit.hit.geomID = RTC_INVALID_GEOMETRY_ID;

                            RTCIntersectContext context;
                            rtcInitIntersectContext(&context);

                            rtcIntersect1(scene, &context, &rayhit);

                            if (rayhit.hit.geomID == RTC_INVALID_GEOMETRY_ID)                             
                                inputM.face[i].Q()+=scalarP;                             
                            else         
                                inputM.face[i].Q()+=(1-powf(rayhit.ray.tfar,tau));
                                                          
                        }
                    }
                }
            }
            tri::UpdateColor<MeshType>::PerFaceQualityGray(inputM);
            rtcReleaseScene(scene);
            rtcReleaseDevice(device);
        }

        public:
         void computeSDF(MeshType &inputM, int nRay){           
            std::vector<Point3f> unifDirVec;
            GenNormal<float>::Fibonacci(nRay,unifDirVec);
            tri::UpdateQuality<MeshType>::FaceConstant(inputM,0);

            #pragma omp parallel 
            {
                #pragma omp for
                for(int i = 0;i<inputM.FN(); i++)
                {   
                    RTCRayHit rayhit;         
                    Point3f b = vcg::Barycenter(inputM.face[i]);
                    rayhit.ray.org_x  = b[0]; rayhit.ray.org_y = b[1]; rayhit.ray.org_z = b[2];  
                    rayhit.ray.tnear  = 1e-4;    
                    int nHits=0;
                    for(int r = 0; r<nRay; r++){
                        Point3f dir = unifDirVec.at(r);                        
                        float scalarP = inputM.face[i].N()*dir;

                        if(scalarP<0){
                            rayhit.ray.dir_x  = dir[0]; rayhit.ray.dir_y = dir[1]; rayhit.ray.dir_z = dir[2];
                            rayhit.ray.tfar   = std::numeric_limits<float>::infinity();
                            rayhit.hit.geomID = RTC_INVALID_GEOMETRY_ID; 

                            RTCIntersectContext context;
                            rtcInitIntersectContext(&context);

                            rtcIntersect1(scene, &context, &rayhit);

                            if (rayhit.hit.geomID != RTC_INVALID_GEOMETRY_ID) 
                            {
                                nHits++;
                                inputM.face[i].Q()+=rayhit.ray.tfar;
                            } 
                                
                        }
                    }

                    inputM.face[i].Q()/=nHits;
                }
            }
            tri::UpdateColor<MeshType>::PerFaceQualityRamp(inputM);
            rtcReleaseScene(scene);
            rtcReleaseDevice(device);
        }

        public:
         std::vector<Point3f> AOBentNormal(MeshType &inputM, int nRay){
            std::vector<Point3f> unifDirVec;
            std::vector<Point3f> bentNormalDir;            
            GenNormal<float>::Fibonacci(nRay,unifDirVec);
            tri::UpdateQuality<MeshType>::FaceConstant(inputM,0);

            #pragma omp parallel 
            {
                #pragma omp for
                for(int i = 0;i<inputM.FN(); i++)
                {       
                    RTCRayHit rayhit;     
                    Point3f b = vcg::Barycenter(inputM.face[i]);
                    rayhit.ray.org_x  = b[0]; rayhit.ray.org_y = b[1]; rayhit.ray.org_z = b[2];  
                    rayhit.ray.tnear  = 1e-4;
                    
                    Point3f bN;
                    int accRays=0;
                    for(int r = 0; r<nRay; r++){
                        Point3f dir = unifDirVec.at(r);
                        
                        float scalarP = inputM.face[i].N()*dir;

                        if(scalarP>0){
                            rayhit.ray.dir_x  = dir[0]; rayhit.ray.dir_y = dir[1]; rayhit.ray.dir_z = dir[2];
                            rayhit.ray.tfar   = std::numeric_limits<float>::infinity();
                            rayhit.hit.geomID = RTC_INVALID_GEOMETRY_ID;   
                            RTCIntersectContext context;
                            rtcInitIntersectContext(&context);

                            rtcIntersect1(scene, &context, &rayhit);

                            if (rayhit.hit.geomID == RTC_INVALID_GEOMETRY_ID) 
                            {                            
                                bN+=dir;
                                accRays++; 
                            } 
                                
                        }
                    }               
                    bentNormalDir.push_back(bN/accRays);
                
                    //inputM.face[i].Q()/=accRays;//inputM.face[i].N()*bentNormalDir;
                }
            }
            //tri::UpdateColor<MeshType>::PerFaceQualityGray(inputM);
            rtcReleaseScene(scene);
            rtcReleaseDevice(device);

            return bentNormalDir;
        }    

    };
}
#endif