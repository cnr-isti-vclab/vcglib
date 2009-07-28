#ifndef OVERLAP_ESTIMATION_H
#define OVERLAP_ESTIMATION_H

#include <vcg/math/gen_normal.h>
#include <vcg/math/random_generator.h>
//#include <vcg/math/point_matching.h>
#include <vcg/space/index/grid_static_ptr.h>
#include <vcg/complex/trimesh/closest.h>
#include <vcg/complex/trimesh/point_sampling.h>
//#include <meshlabplugins/edit_pickpoints/pickedPoints.h>

#include <qdatetime.h>

using namespace std;
using namespace vcg;


template<class MESH_TYPE> class OverlapEstimation
{
    public:

    typedef MESH_TYPE MeshType;
    typedef typename MeshType::ScalarType ScalarType;
    typedef typename MeshType::CoordType CoordType;
    typedef typename MeshType::VertexType VertexType;
    typedef typename MeshType::FaceType FaceType;
    typedef typename MeshType::VertexPointer VertexPointer;
    typedef typename MeshType::VertexIterator VertexIterator;
    typedef typename vector<VertexPointer>::iterator VertexPointerIterator;
    typedef GridStaticPtr<VertexType, ScalarType > MeshGrid;
    typedef tri::VertTmark<MeshType> MarkerVertex;

    private:
    class VertexPointerSampler
    {
        public:

        MeshType* m;  //this is needed for advanced sampling (i.e poisson sampling)

        VertexPointerSampler(){ m = new MeshType(); m->Tr.SetIdentity(); m->sfn=0; }
        ~VertexPointerSampler(){ if(m) delete m; }
        vector<VertexType*> sampleVec;

        void AddVert(VertexType &p){ sampleVec.push_back(&p); }

        void AddFace(const FaceType &f, const CoordType &p){}

        void AddTextureSample(const FaceType &, const CoordType &, const Point2i &){}
    };

    public:
    class Parameters
    {
        public:
        int samples;                                //number of samples
        int bestScore;                              //used to paint mMove: paint only if is a best score
        float consensusDist;                        //consensus distance
        float consensusNormalsAngle;                //holds the the consensus angle for normals, in gradients.
        float threshold;                            //consensus % to win consensus;
        bool normalEqualization;                    //to use normal equalization in consensus
        bool paint;                                 //to paint mMov according to consensus
        void (*log)(int level, const char * f, ... );          //pointer to log function

        Parameters()
        {
            samples = 2500;
            bestScore = 0;
            consensusDist = 2.0f;
            consensusNormalsAngle = 0.965f;   //15 degrees.
            threshold = 0.0f;
            normalEqualization = true;
            paint = false;
            log = NULL;
        }
    };

    MeshType* mFix, *mMov;
    vector<vector<int> >* normBuckets;          //structure to hold normals bucketing. Needed for normal equalized sampling during consensus
    MeshGrid* gridFix;                          //variable to manage uniform grid
    MarkerVertex markerFunctorFix;              //variable to manage uniform grid

    OverlapEstimation() : normBuckets(NULL), gridFix(NULL){}
    ~OverlapEstimation(){
        if(normBuckets) delete normBuckets;
        if(gridFix) delete gridFix;
    }

    void SetFix(MeshType& m){ mFix = &m; }

    void SetMove(MeshType& m){ mMov = &m; }

    void Paint()
    {
        for(VertexIterator vi=mMov->vert.begin(); vi!=mMov->vert.end(); vi++){
            if(!(*vi).IsD()){
                if((*vi).Q()==0.0) (*vi).C() = Color4b::Red;
                if((*vi).Q()==1.0) (*vi).C() = Color4b::Yellow;
                if((*vi).Q()==2.0) (*vi).C() = Color4b::Blue;
            }
        }
    }

    bool Init(Parameters& param){
        //builds the uniform grid with mFix vertices
        gridFix = new MeshGrid();
        SetupGrid();

        //if requested, group normals of mMov into 30 buckets. Buckets are used for Vertex Normal Equalization
        //in consensus. Bucketing is done here once for all to speed up consensus.
        if(normBuckets) {normBuckets->clear(); delete normBuckets; }
        if(param.normalEqualization){
            normBuckets = BucketVertexNormal(mMov->vert, 30);
            assert(normBuckets);
        }
        return true;
    }

    float Compute(Parameters& param)
    {
        return Check(param)/float(param.samples);
    }

    //compute the randomized consensus beetween m1 e m2 (without taking in account any additional transformation)
    //IMPORTANT: per vertex normals of m1 and m2 MUST BE PROVIDED JET NORMALIZED!!!
    int Check(Parameters& param)
    {
        //pointer to a function to compute distance beetween points
        vertex::PointDistanceFunctor<ScalarType> PDistFunct;

        //if no buckets are provided get a vector of vertex pointers sampled uniformly
        //else, get a vector of vertex pointers sampled in a normal equalized manner; used as query points
        vector<VertexPointer> queryVert;
        if(param.normalEqualization){
            assert(normBuckets);
            for(unsigned int i=0; i<mMov->vert.size(); i++) queryVert.push_back(&(mMov->vert[i]));//do a copy of pointers to vertexes
            SampleVertNormalEqualized(queryVert, param.samples);
        }
        else{
            SampleVertUniform(*mMov, queryVert, param.samples);
        }
        assert(queryVert.size()!=0);

        //init variables for consensus
        float consDist = param.consensusDist*(mMov->bbox.Diag()/100.0f);  //consensus distance
        int cons_succ = int(param.threshold*(param.samples/100.0f));      //score needed to pass consensus
        int consensus = 0;                  //counts vertices in consensus
        float dist;                         //holds the distance of the closest vertex found
        VertexType* closestVertex = NULL;   //pointer to the closest vertex
        Point3<ScalarType> queryNrm;        //the query point normal for consensus
        CoordType queryPnt;                 //the query point for consensus
        CoordType closestPnt;               //the closest point found in consensus
        Matrix33<ScalarType> inv33_matMov(mMov->Tr,3);          //3x3 matrix needed to transform normals
        Matrix33<ScalarType> inv33_matFix(Inverse(mFix->Tr),3); //3x3 matrix needed to transform normals

        //consensus loop
        VertexPointerIterator vi; int i;
        for(i=0, vi=queryVert.begin(); vi!=queryVert.end(); vi++, i++)
        {
            dist = -1.0f;
            //set query point; vertex coord is transformed properly in fix mesh coordinates space; the same for normals
            queryPnt = Inverse(mFix->Tr) * (mMov->Tr * (*vi)->P());
            queryNrm = inv33_matFix * (inv33_matMov * (*vi)->N());
            //if query point is bbox, the look for a vertex in cDist from the query point
            if(mFix->bbox.IsIn(queryPnt)) closestVertex = gridFix->GetClosest(PDistFunct,markerFunctorFix,queryPnt,consDist,dist,closestPnt);
            else closestVertex=NULL;  //out of bbox, we consider the point not in consensus...

            if(closestVertex!=NULL && dist < consDist){
                assert(closestVertex->P()==closestPnt); //coord and vertex pointer returned by getClosest must be the same

                //point is in consensus distance, now we check if normals are near
                if(queryNrm.dot(closestVertex->N())>param.consensusNormalsAngle)  //15 degrees
                {
                    consensus++;  //got consensus
                    if(param.paint) (*vi)->Q() = 0.0f;  //store 0 as quality
                }
                else{
                    if(param.paint) (*vi)->Q() = 1.0f;  //store 1 as quality
                }
            }
            else{
                if(param.paint) (*vi)->Q() = 2.0f;  //store 2 as quality
            }
        }

        //apply colors if consensus is the best ever found.
        //NOTE: we got to do this here becouse we need a handle to sampler. This is becouse vertex have been shuffled
        //and so colors have been stored not in order in the buffer!
        if(param.paint){
            if(consensus>=param.bestScore && consensus>=cons_succ) Paint();
        }

        return consensus;
    }

    private:

    void SampleVertUniform(MESH_TYPE& m, vector<typename MESH_TYPE::VertexPointer>& vert, int sampleNum)
    {
        VertexPointerSampler sampler;
        tri::SurfaceSampling<MeshType, VertexPointerSampler>::VertexUniform(m, sampler, sampleNum);
        for(unsigned int i=0; i<sampler.sampleVec.size(); i++) vert.push_back(sampler.sampleVec[i]);
    }

    vector<vector<int> >* BucketVertexNormal(typename MESH_TYPE::VertContainer& vert, int bucketDim = 30)
    {
        static vector<Point3f> NV;
        if(NV.size()==0) GenNormal<float>::Uniform(bucketDim,NV);

        // Bucket vector dove, per ogni normale metto gli indici
        // dei vertici ad essa corrispondenti
        vector<vector<int> >* BKT = new vector<vector<int> >(NV.size()); //NV size is greater then bucketDim, so don't change this!

        int ind;
        for(int i=0;i<vert.size();++i){
            ind=GenNormal<float>::BestMatchingNormal(vert[i].N(),NV);
            (*BKT)[ind].push_back(i);
        }

        return BKT;
    }

    bool SampleVertNormalEqualized(vector<typename MESH_TYPE::VertexPointer>& vert, int SampleNum)
    {
        assert(normBuckets);
        // vettore di contatori per sapere quanti punti ho gia' preso per ogni bucket
        vector<int> BKTpos(normBuckets->size(),0);

        if(SampleNum >= int(vert.size())) SampleNum= int(vert.size()-1);

        int ind;
        for(int i=0;i<SampleNum;){
            ind=LocRnd(normBuckets->size()); // Scelgo un Bucket
            int &CURpos = BKTpos[ind];
            vector<int> &CUR = (*normBuckets)[ind];

            if(CURpos<int(CUR.size())){
                swap(CUR[CURpos], CUR[ CURpos + LocRnd((*normBuckets)[ind].size()-CURpos)]);
                swap(vert[i],vert[CUR[CURpos]]);
                ++BKTpos[ind];
                ++i;
            }
        }

        vert.resize(SampleNum);
        return true;
    }

    static math::SubtractiveRingRNG &LocRnd(){
        static math::SubtractiveRingRNG myrnd(time(NULL));
        return myrnd;
    }

    static int LocRnd(int n){
        return LocRnd().generate(n);
    }

    inline void SetupGrid()
    {
        gridFix->Set(mFix->vert.begin(),mFix->vert.end());
        markerFunctorFix.SetMesh(mFix);
    }
};

#endif // OVERLAP_ESTIMATION_H
