#ifndef __VCGFOREMBREE_H
#define __VCGFOREMBREE_H

#include <iostream>
#include <vcg/complex/complex.h>
#include <vcg/simplex/face/pos.h>
#include <vcg/complex/algorithms/update/topology.h>
#include <vcg/complex/algorithms/closest.h>

#include <vcg/complex/algorithms/update/color.h>
#include <vcg/complex/algorithms/update/quality.h>
#include <wrap/callback.h>
#include <embree4/rtcore.h>
#include <vcg/math/gen_normal.h>
#include <limits>
#include <math.h>
#include <time.h>
#include <omp.h>
#include <tuple>

using namespace vcg;
using namespace std;




/*
    @Author: Paolo Fasano
    @Description: This class aims to integrate intel embree3 with the vcglib giving some basic methods that can be used to build more complex features.
*/
namespace vcg{
    template <class MeshType>
    class EmbreeAdaptor{

        RTCDevice device = rtcNewDevice(NULL);
        RTCScene scene = rtcNewScene(device);
        RTCGeometry geometry = rtcNewGeometry(device, RTC_GEOMETRY_TYPE_TRIANGLE);
        int threads;

        public:
         EmbreeAdaptor(MeshType &m){
                loadVCGMeshInScene(m);
            }

        /*
        @Author: Paolo Fasano
        @Parameter: Point3f rayDirection, direction the rays are shoot towards
        @Description: foreach face the barycenter is found and a single ray is shoot. If the ray intersect with
            something the face color is set to black else is set to white.
        */
        public:
         void selectVisibleFaces(MeshType &m, Point3f rayDirection, bool incrementalSelect){
            
            if (incrementalSelect == false){
                //deselect all previously selected faces
                for(int i = 0;i<m.FN(); i++){
                    if(m.face[i].IsS()){
                        m.face[i].ClearS();
                    }
                }
            }

            Point3f normalizedDir = rayDirection;

            RTCRayHit rayhit = initRayValues();

            for(int i = 0;i<m.FN(); i++)
            {
                Point3f b = vcg::Barycenter(m.face[i]);
                std::vector<Point3f> unifDirVec;
                Point3f dir = normalizedDir;

                rayhit = setRayValues(b, dir, 4);

                RTCRayQueryContext context;
                rtcInitRayQueryContext(&context);

                RTCIntersectArguments intersectArgs;
                rtcInitIntersectArguments(&intersectArgs);
                intersectArgs.context = &context;

                rtcIntersect1(scene, &rayhit, &intersectArgs);

                //if (rayhit.hit.geomID != RTC_INVALID_GEOMETRY_ID)
                if (rayhit.ray.tfar == std::numeric_limits<float>::infinity())
                    m.face[i].SetS();

            }
            rtcReleaseScene(scene);
            rtcReleaseDevice(device);

        }

        

        /*
        @Author: Paolo Fasano
        @Parameter: MeshType &m, reference to a mesh
        @Description: this method apply some preprocessing over it using standard vcglib methods.
            Than the mesh is loaded as a new embree geometry inside a new embree scene. The new embree variables
            are global to the class in order to be used with the other methods.
        */
        public:
         void loadVCGMeshInScene(MeshType &m){
            //a little mesh preprocessing before adding it to a RTCScene
            tri::RequirePerVertexNormal(m);
            tri::UpdateNormal<MeshType>::PerVertexNormalized(m);
            tri::UpdateNormal<MeshType>::PerFaceNormalized(m);
            tri::UpdateBounding<MeshType>::Box(m);
            tri::UpdateFlags<MeshType>::FaceClearV(m);

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

        /*
        @Author: Paolo Fasano
        @Parameter: MeshType &m, reference to a mesh.
        @Parameter: int nRay, number of rays that must be generated and shoot.
        @Description: for each face from the barycenter this method shoots n rays towards a generated direction(to infinity).
            If the ray direction is not pointing inside than the ray is actually shoot.
            If the ray intersect something than the face quality of the mesh is updated with the normal of the fica multiplied by the direction.
        */
        public:
         void computeAmbientOcclusion(MeshType &inputM, int nRay){
            std::vector<Point3f> unifDirVec;
            GenNormal<float>::Fibonacci(nRay,unifDirVec);
            computeAmbientOcclusion(inputM, unifDirVec);
        }

        /*
        @Author: Paolo Fasano
        @Parameter: MeshType &m,reference to a mesh.
        @Parameter: std::vector<Point3f> unifDirVec, vector of direction specified by the user.
        @Description: for each face from the barycenter this method shoots n rays towards some user generated directions(to infinity).
            If the ray direction is not pointing inside than the ray is actually shoot.
            If the ray intersect something than the face quality of the mesh is updated with the normal of the fica multiplied by the direction.

            One more operation done in the AmbientOcclusion is to calculate the bent normal foreach face and save it in an attribute named "BentNormal"
        */
        public:
         void computeAmbientOcclusion(MeshType &inputM, std::vector<Point3f> unifDirVec){
            tri::UpdateQuality<MeshType>::FaceConstant(inputM,0);
            typename MeshType::template PerFaceAttributeHandle<Point3f> bentNormal = vcg::tri::Allocator<MeshType>:: template GetPerFaceAttribute<Point3f>(inputM,string("BentNormal"));

            #pragma omp parallel shared(inputM)
            {
                #pragma omp for
                for(int i = 0;i<inputM.FN(); i++)
                {
                    RTCRayHit rayhit = initRayValues();
                    Point3f b = vcg::Barycenter(inputM.face[i]);
                    updateRayOrigin(rayhit, b);
                    rayhit.ray.tnear  = 0.00001f;

                    Point3f bN;
                    int accRays=0;
                    for(int r = 0; r<unifDirVec.size(); r++){
                        Point3f dir = unifDirVec.at(r);
                        float scalarP = inputM.face[i].N()*dir;

                        if(scalarP>0){

                            updateRayDirection(rayhit, dir);
                            rayhit.ray.tfar   = std::numeric_limits<float>::infinity();
                            rayhit.hit.geomID = RTC_INVALID_GEOMETRY_ID;

                            RTCRayQueryContext context;
                            rtcInitRayQueryContext(&context);

                            RTCIntersectArguments intersectArgs;
                            rtcInitIntersectArguments(&intersectArgs);
                            intersectArgs.context = &context;

                            rtcIntersect1(scene, &rayhit, &intersectArgs);

                            if (rayhit.hit.geomID == RTC_INVALID_GEOMETRY_ID){
                                bN+=dir;
                                accRays++;
                                inputM.face[i].Q()+=scalarP;
                            }

                        }
                    }
                    bentNormal[i] = bN/accRays;
                }
            }
            tri::UpdateColor<MeshType>::PerFaceQualityGray(inputM);
            rtcReleaseScene(scene);
            rtcReleaseDevice(device);
        }

        /*
        @Author: Paolo Fasano
        @Parameter: MeshType &m, reference to a mesh.
        @Parameter: int nRay, number of rays that must be generated and shoot.
        @Parameter: float Tau, the grater this value is the grater the influence of the rays that intersect with some face
        @Description: for each face from the barycenter this method shoots n rays towards a generated direction(to infinity).
            If the ray direction is not pointing inside than the ray is actually shoot.
            If the ray intersect something than the face quality of the mesh is updated with the normal of the fica multiplied by the direction;
            else, if there are no hits, the face get updated of 1-distanceHit^tau
        */
        public:
         void computeObscurance(MeshType &inputM, int nRay, float tau){
            std::vector<Point3f> unifDirVec;
                GenNormal<float>::Fibonacci(nRay,unifDirVec);

           computeObscurance(inputM, unifDirVec, tau);
        }

        /*
        @Author: Paolo Fasano
        @Parameter: MeshType &m, reference to a mesh.
        @Parameter: std::vector<Point3f> unifDirVec, vector of direction specified by the user.
        @Parameter: float Tau, the grater this value is the grater the influence of the rays that intersect with some face
        @Description: for each face from the barycenter this method shoots n rays towards a generated direction(to infinity).
            If the ray direction is not pointing inside than the ray is actually shoot.
            If the ray intersect something than the face quality of the mesh is updated with the normal of the fica multiplied by the direction;
            else, if there are no hits, the face get updated of 1-distanceHit^tau
        */
        public:
         void computeObscurance(MeshType &inputM, std::vector<Point3f> unifDirVec, float tau){
            tri::UpdateQuality<MeshType>::FaceConstant(inputM,0);

            #pragma omp parallel
            {
                #pragma omp for
                for(int i = 0;i<inputM.FN(); i++)
                {
                    RTCRayHit rayhit = initRayValues();
                    Point3f b = vcg::Barycenter(inputM.face[i]);
                    updateRayOrigin(rayhit, b);
                    rayhit.ray.tnear  = 0.00001f;

                    for(int r = 0; r<unifDirVec.size(); r++){
                        Point3f dir = unifDirVec.at(r);
                        float scalarP = inputM.face[i].N()*dir;

                        if(scalarP>0){
                            updateRayDirection(rayhit, dir);
                            rayhit.ray.tfar   = std::numeric_limits<float>::infinity();
                            rayhit.hit.geomID = RTC_INVALID_GEOMETRY_ID;

                            RTCRayQueryContext context;
                            rtcInitRayQueryContext(&context);

                            RTCIntersectArguments intersectArgs;
                            rtcInitIntersectArguments(&intersectArgs);
                            intersectArgs.context = &context;

                            rtcIntersect1(scene, &rayhit, &intersectArgs);

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

        /*
        @Author: Paolo Fasano
        @Parameter: MeshType &m, reference to a mesh.
        @Parameter: int nRay, number of rays that must be generated and shoot.
        @Parameter: float degree, this variable represents the angle of the cone for which we consider a point as a valid direction
        @Description:
            - Use a cone centered around the inward-normal direction of a point on a surface mesh.
            - Send rays inside the cone to the other side of the mesh.
            - Calculate the SDF at a point as the weighted average of all ray lengths within one
                standard deviation from the median of all lengths.
            - Use inverse angle between the ray and the center of the cone as the weight.
            - The SDF is invariant to rigid body transformations of the whole mesh and oblivious to local deformations.
            - Small cone angles are too sensitive to local features, while large opening angles close to 180◦
                expose the SDF measure to noise and errors.

        */
        public:
         void computeSDF(MeshType &inputM, int nRay, float degree){

            if (degree >= 180)
                degree = 120;

            tri::UpdateQuality<MeshType>::FaceConstant(inputM,0);

            //first step is to generate the cone of rays, i will generate a sphear and for each face take only the
            //the ones that fall in the desired degree (in respect to the opposite of face norm)
            std::vector<Point3f> unifDirVec;
            GenNormal<float>::Fibonacci(nRay,unifDirVec);

            for (int i = 0; i < inputM.FN(); i++)
            {
                RTCRayHit rayhit = initRayValues();
                Point3f b = vcg::Barycenter(inputM.face[i]);
                updateRayOrigin(rayhit, b);
                rayhit.ray.tnear  = 1e-4;

                float weight = 0;
                float weight_sum = 0;
                float weighted_sum = 0;

                for(int r = 0; r<unifDirVec.size(); r++){
                    Point3f dir = unifDirVec.at(r);
                    float scalarP = inputM.face[i].N()*dir;

                    float angle_dir_b = Angle(b, dir);

                    if (scalarP < 0 && vcg::math::ToRad(angle_dir_b) <= degree){
                        updateRayDirection(rayhit, dir);
                        rayhit.ray.tfar   = std::numeric_limits<float>::infinity();
                        rayhit.hit.geomID = RTC_INVALID_GEOMETRY_ID;

                        RTCRayQueryContext context;
                        rtcInitRayQueryContext(&context);

                        RTCIntersectArguments intersectArgs;
                        rtcInitIntersectArguments(&intersectArgs);
                        intersectArgs.context = &context;

                        rtcIntersect1(scene, &rayhit, &intersectArgs);

                        if (rayhit.ray.tfar != std::numeric_limits<float>::infinity())
                        {
                            weight = 1/angle_dir_b;
                            //nominator and denominator for weigthed avarage
                            weighted_sum += weight * rayhit.ray.tfar;
                            weight_sum += weight;

                        }
                    }
                }
                //we assign the result of the weighted average to the quality of the face
                inputM.face[i].Q() = (weighted_sum / weight_sum);
            }

            tri::UpdateColor<MeshType>::PerFaceQualityRamp(inputM);
            rtcReleaseScene(scene);
            rtcReleaseDevice(device);
        }

        /*
        @Author: Paolo Fasano
        @Parameter: MeshType &m, reference to a mesh.
        @Parameter: int nRay, number of rays that must be generated and shoot.
        @Description: Given a mesh for each face, for each ray, it detects all the intersections with all the facets
            (i.e., without stopping at the first intersection), and accumulates the number.
            The rays are shoot two times each, one inside and one outside the mesh.
            After shooting all the rays, the facet is flipped if frontHit>backHit.

            For more informations read:
            Kenshi Takayama, Alec Jacobson, Ladislav Kavan, and Olga Sorkine-Hornung, A Simple Method for Correcting Facet Orientations in Polygon Meshes Based on Ray Casting, Journal of Computer Graphics Techniques (JCGT), vol. 3, no. 4, 53-63, 2014
            Available online http://jcgt.org/published/0003/04/02/
        */
        public:
         void computeNormalAnalysis(MeshType &inputM, int nRay,bool parity_computation){

            //bool fast_computation = false;

            std::vector<Point3f> unifDirVec;
            GenNormal<float>::Fibonacci(nRay,unifDirVec);

            std::vector<Point3f> unifDirVec_EXPAND;
            for (Point3f p : unifDirVec)
            {
                unifDirVec_EXPAND.push_back(p);
                unifDirVec_EXPAND.push_back(p*-1);
            }

            tri::UpdateSelection<MeshType>::FaceClear(inputM);


            if (parity_computation){
                paritySampling(inputM,unifDirVec_EXPAND);
            }
            else{
                visibilitySamplig(inputM, unifDirVec_EXPAND);
            }

            // Iterate over the selected faces and flip them
            for (auto& face : inputM.face)
            {
                if (face.IsS())
                {
                    std::swap(face.V(1), face.V(2));  // Swap the second and third vertices
                }
            }
            tri::UpdateSelection<MeshType>::FaceClear(inputM);

        }


        /*
        @Author: Paolo Fasano
        @Parameter: Point3f origin, the origin point to search for intersections to
        @Description: given an origin point this methos counts how many intersection there are starting from there
            (only the number not the positions coordinates)
        */
        public:
         int findInterceptNumber(Point3f origin, Point3f direction){

            int totInterception = 0;
            //float totDistance = 0;
            //float previous_distance = 0;

            RTCRayHit rayhit;
            rayhit = setRayValues(origin, direction, 0.5);
            RTCRayQueryContext context;
            rtcInitRayQueryContext(&context);

            RTCIntersectArguments intersectArgs;
            rtcInitIntersectArguments(&intersectArgs);
            intersectArgs.context = &context;

            while(true){
                rtcIntersect1(scene, &rayhit, &intersectArgs);
                if (rayhit.ray.tfar != std::numeric_limits<float>::infinity()){
                    totInterception++;
                    //totDistance += rayhit.ray.tfar - previous_distance;
                    //previous_distance = rayhit.ray.tfar;

                    // we keep the same origin point and direction but update how far from the face the ray starts shooting
                    rayhit.ray.tnear += rayhit.ray.tfar ;//+ rayhit.ray.tfar * 0.05f;
                    rayhit.ray.tfar = std::numeric_limits<float>::infinity();
                }
                else
                    return totInterception;
            }

             return totInterception;
        }


        public:
            void visibilitySamplig(MeshType &inputM, std::vector<Point3f> unifDirVec_EXPAND){

            #pragma omp parallel
            {
                #pragma omp for
                for(int i = 0;i<inputM.FN(); i++)
                {
                    RTCRayHit rayhit;
                    Point3f b = vcg::Barycenter(inputM.face[i]);
                    Point3f dir(0.0f, 0.0f, 0.0f);
                    rayhit = setRayValues(b, dir, 1e-4f);

                    int frontHit = 0;
                    int backHit = 0;

                    float frontDistance = 0;
                    float backDistance = 0;

                    for(int r = 0; r<unifDirVec_EXPAND.size(); r++){
                        dir = unifDirVec_EXPAND.at(r);
                        float scalarP = inputM.face[i].N()*dir;

                        rayhit = setRayValues(b, dir, 1e-4f);
                        RTCRayQueryContext context;
                        rtcInitRayQueryContext(&context);

                        RTCIntersectArguments intersectArgs;
                        rtcInitIntersectArguments(&intersectArgs);
                        intersectArgs.context = &context;

                        rtcIntersect1(scene, &rayhit, &intersectArgs);

                        if (rayhit.ray.tfar  == std::numeric_limits<float>::infinity()) {

                            if (scalarP > 0){
                                frontHit++;
                                frontDistance += rayhit.ray.tfar;
                            }
                            else{
                                backHit++;
                                backDistance += rayhit.ray.tfar;
                            }
                        }
                    }

                    //std::cout<< "face "<< i <<"front hit: " << frontHit << " backhit "<< backHit << " frontDistance " << frontDistance << " backDistance " << backDistance << endl;
                    if(frontHit  < backHit || (frontHit  == backHit && frontDistance < backDistance))
                        inputM.face[i].SetS();

                }
            }

            //tri::Clean<MeshType>::FlipMesh(inputM,true);
            rtcReleaseScene(scene);
            rtcReleaseDevice(device);
            return;
        }

        public:
        void paritySampling(MeshType &inputM, std::vector<Point3f> unifDirVec_EXPAND){

            #pragma omp parallel
            {
                #pragma omp for
                for(int i = 0;i<inputM.FN(); i++)
                {
                    RTCRayHit rayhit;
                    Point3f b = vcg::Barycenter(inputM.face[i]);
                    Point3f dir(0.0f, 0.0f, 0.0f);
                    rayhit = setRayValues(b, dir, 1e-4f);

                    int frontHit = 0;
                    int backHit = 0;

                    for(int r = 0; r<unifDirVec_EXPAND.size(); r++){
                        dir = unifDirVec_EXPAND.at(r);
                        float scalarP = inputM.face[i].N()*dir;

                        rayhit = setRayValues(b, dir, 1e-4f);
                        RTCRayQueryContext context;
                        rtcInitRayQueryContext(&context);

                        RTCIntersectArguments intersectArgs;
                        rtcInitIntersectArguments(&intersectArgs);
                        intersectArgs.context = &context;

                        rtcIntersect1(scene, &rayhit, &intersectArgs);

                        if (rayhit.ray.tfar  != std::numeric_limits<float>::infinity()) {

                            int n_hits = findInterceptNumber(b,dir);

                            if (scalarP > 0){
                                frontHit += (n_hits % 2);
                            }
                            else{
                                backHit += (n_hits % 2);
                            }
                        }
                    }

                    if(frontHit > backHit )
                        inputM.face[i].SetS();

                }
            }

            //tri::Clean<MeshType>::FlipMesh(inputM,true);
            rtcReleaseScene(scene);
            rtcReleaseDevice(device);
            return;
        }

        public:
            inline RTCRayHit initRayValues(){
                RTCRayHit rayhit;
                rayhit.ray.mask   = -1;
                rayhit.ray.flags = 0;
                rayhit.hit.geomID = RTC_INVALID_GEOMETRY_ID;
                rayhit.hit.instID[0] = RTC_INVALID_GEOMETRY_ID;
                return rayhit;
            }

        //given a ray and a direction expressed as point3f, this method modifies the ray direction of the ray tp the given direction
        public:
            inline void updateRayDirection(RTCRayHit& rayhit, Point3f direction){

                //setting the ray direction
                rayhit.ray.dir_x = direction[0];
                rayhit.ray.dir_y = direction[1];
                rayhit.ray.dir_z = direction[2];
            }


        //given a ray and a point of origin expressed as point3f, this method modifies the origin point of the ray tp the origin point given
        public:
            inline void updateRayOrigin(RTCRayHit& rayhit, Point3f origin){

                //setting the ray point of origin
                rayhit.ray.org_x = origin[0];
                rayhit.ray.org_y = origin[1];
                rayhit.ray.org_z = origin[2];
            }

        public:
            inline RTCRayHit setRayValues(Point3f origin, Point3f direction, float tnear, float tfar = std::numeric_limits<float>::infinity()){

                RTCRayHit rayhit = initRayValues();

                updateRayOrigin(rayhit, origin);
                updateRayDirection(rayhit, direction);
                rayhit.ray.tnear  = tnear;
                rayhit.ray.tfar   = tfar;

                return rayhit;
            }
    };
}
#endif
