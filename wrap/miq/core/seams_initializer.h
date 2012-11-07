#ifndef MIQ_SEAMS_INTIALIZER
#define MIQ_SEAMS_INTIALIZER
#include <vcg/complex/complex.h>
#include <vcg/simplex/face/pos.h>
#include <vcg/simplex/face/jumping_pos.h>
#include <vcg/simplex/face/topology.h>
#include <wrap/io_trimesh/import_field.h>

template <class MeshType>
class SeamsInitializer
{

private:
    typedef typename MeshType::ScalarType ScalarType;
    typedef typename MeshType::FaceType FaceType;
    typedef typename MeshType::VertexType VertexType;
    typedef typename MeshType::CoordType CoordType;
    typedef typename MeshType::FaceIterator FaceIterator;

    MeshType *mesh;

    ///per face per edge of mmatch in the solver
    typename  MeshType::template PerFaceAttributeHandle<vcg::Point3i> Handle_MMatch;
    ///per vertex singular or not
    typename  MeshType::template PerVertexAttributeHandle<bool> Handle_Singular;
    ///per vertex degree of a singularity
    typename  MeshType::template PerVertexAttributeHandle<int> Handle_SingularDegree;
    ///seam per face
    typename  MeshType::template PerFaceAttributeHandle<vcg::Point3<bool> > Handle_Seams;

    bool IsRotSeam(const FaceType *f0,const int edge)
    {
        unsigned char MM=Handle_MMatch[f0][edge];//MissMatch(f0,edge);
        return (MM!=0);
    }

    ///return true if a vertex is singluar by looking at initialized missmatches
    bool IsSingularByMMatch(const VertexType &v,int &missmatch)
    {
        ///check that is on border..
        if (v.IsB())return false;

        std::vector<FaceType*> faces;
        std::vector<int> edges;
        vcg::face::VFOrderedStarVF_FF<FaceType>(v,faces,edges);

        missmatch=0;
        for (unsigned int i=0;i<faces.size();i++)
        {
            FaceType *curr_f=faces[i];
            int currMM=Handle_MMatch[curr_f][edges[i]];
            missmatch+=currMM;
        }
        missmatch=missmatch%4;
        return(missmatch!=0);
    }

    ///initialized mapping structures if are not already initialized
    void AddAttributesIfNeeded()
    {

        bool HasHandleMMatch=vcg::tri::HasPerFaceAttribute(*mesh,std::string("MissMatch"));
        if (!HasHandleMMatch)
            Handle_MMatch = vcg::tri::Allocator<MeshType>::template AddPerFaceAttribute<vcg::Point3i>(*mesh,std::string("MissMatch"));
        else
            Handle_MMatch = vcg::tri::Allocator<MeshType>::template GetPerFaceAttribute<vcg::Point3i>(*mesh,std::string("MissMatch"));

        bool HasHandleSingularity=vcg::tri::HasPerVertexAttribute(*mesh,std::string("Singular"));
        if (!HasHandleSingularity)
            Handle_Singular=vcg::tri::Allocator<MeshType>::template AddPerVertexAttribute<bool>(*mesh,std::string("Singular"));
        else
            Handle_Singular=vcg::tri::Allocator<MeshType>::template GetPerVertexAttribute<bool>(*mesh,std::string("Singular"));

        bool HasHandleSingularityDegree=vcg::tri::HasPerVertexAttribute(*mesh,std::string("SingularityDegree"));
        if (!HasHandleSingularityDegree)
            Handle_SingularDegree=vcg::tri::Allocator<MeshType>::template AddPerVertexAttribute<int>(*mesh,std::string("SingularityDegree"));
        else
            Handle_SingularDegree=vcg::tri::Allocator<MeshType>::template GetPerVertexAttribute<int>(*mesh,std::string("SingularityDegree"));

        bool HasHandleSeams=vcg::tri::HasPerFaceAttribute(*mesh,std::string("Seams"));
        if (!HasHandleSeams)
            Handle_Seams=vcg::tri::Allocator<MeshType>::template AddPerFaceAttribute<vcg::Point3<bool> >(*mesh,std::string("Seams"));
        else
            Handle_Seams=vcg::tri::Allocator<MeshType>::template GetPerFaceAttribute<vcg::Point3<bool> >(*mesh,std::string("Seams"));
    }


    void FloodFill(FaceType* start)
    {
        std::deque<FaceType*> d;
        ///clean the visited flag
        start->SetV();
        d.push_back(start);

        while (!d.empty()){
            FaceType *f = d.at(0); d.pop_front();
            for (int s = 0; s<3; s++)
            {
                FaceType *g = f->FFp(s);
                int j = f->FFi(s);
                if ((!(IsRotSeam(f,s))) && (!(IsRotSeam(g,j)))  && (!g->IsV()) )
                {
                    Handle_Seams[f][s]=false;
                    Handle_Seams[g][j]=false;
                    g->SetV();
                    d.push_back(g);
                }
            }
        }
    }

    void Retract(){
        std::vector<int> e(mesh->vert.size(),0); // number of edges per vert
        VertexType *vb = &(mesh->vert[0]);
        for (FaceIterator f = mesh->face.begin(); f!=mesh->face.end(); f++) if (!f->IsD()){
            for (int s = 0; s<3; s++){
                //if (f->seam[s])
                if (Handle_Seams[f][s])
                    if (f->FFp(s)<=&*f)  {
                        e[ f->V(s) - vb ] ++;
                        e[ f->V1(s) - vb ] ++;
                    }
            }
        }
        bool over=true;
        int guard = 0;
        do {
            over = true;
            for (FaceIterator f = mesh->face.begin(); f!=mesh->face.end(); f++) if (!f->IsD()){
                for (int s = 0; s<3; s++){
                    //if (f->seam[s])
                    if (Handle_Seams[f][s])
                        if (!(IsRotSeam(&(*f),s))) // never retract rot seams
                            //if (f->FFp(s)<=&*f)
                        {
                            if (e[ f->V(s) - vb ] == 1) {
                                // dissolve seam
                                //f->seam[s] = false;
                                Handle_Seams[f][s]=false;
                                //f->FFp(s)->seam[(int)f->FFi(s)] = false;
                                Handle_Seams[f->FFp(s)][(int)f->FFi(s)]=false;
                                e[ f->V(s) - vb ] --;
                                e[ f->V1(s) - vb ] --;
                                over = false;
                            }
                        }
                }
            }

            if (guard++>10000) over = true;

        } while (!over);
    }

    void AddSeamsByMM()
    {
        for (unsigned int i=0;i<mesh->face.size();i++)
        {
            FaceType *f=&mesh->face[i];
            if (f->IsD())continue;
            for (int j=0;j<3;j++)
            {
                if (IsRotSeam(f,j))
                    Handle_Seams[f][j]=true;
                    //f->SetSeam(j);
            }
        }
    }

    void SelectSingularityByMM()
    {
        vcg::tri::UpdateFlags<MeshType>::VertexBorderFromNone(*mesh);
        for (unsigned int i=0;i<mesh->vert.size();i++)
        {
            if (mesh->vert[i].IsD())continue;
            int missmatch;
            bool isSing=IsSingularByMMatch(mesh->vert[i],missmatch);
            if (isSing)
            {
                //mesh->vert[i].SetS();
                Handle_Singular[i]=true;
                Handle_SingularDegree[i]=missmatch;
            }
            else
            {
                //mesh->vert[i].ClearS();
                Handle_Singular[i]=false;
                Handle_SingularDegree[i]=0;
            }
        }
    }


    int InitTopologycalCuts(){
        vcg::tri::UpdateFlags<MeshType>::FaceClearV(*mesh);

        for (FaceIterator f = mesh->face.begin(); f!=mesh->face.end(); f++)
            if (!f->IsD())
            {
                Handle_Seams[f][0]=true;
                Handle_Seams[f][1]=true;
                Handle_Seams[f][2]=true;
            }

        int index=0;
        for (FaceIterator f = mesh->face.begin(); f!=mesh->face.end(); f++)
            if (!f->IsD())
            {
                if (!f->IsV())
                {
                    index++;
                    FloodFill(&*f);
                }
            }

        Retract();
        return index;
    }

    void InitMMatch()
    {
        for (unsigned int i=0;i<mesh->face.size();i++)
        {
            FaceType *curr_f=&mesh->face[i];
            for (int j=0;j<3;j++)
            {
                FaceType *opp_f=curr_f->FFp(j);
                if (curr_f==opp_f)
                    Handle_MMatch[curr_f][j]=0;
                else
                    Handle_MMatch[curr_f][j]=vcg::tri::CrossField<MeshType>::MissMatchByCross(*curr_f,*opp_f);
            }
        }
    }

public:


    void Init(MeshType *_mesh,
              bool orient_globally,
              bool initMM,
              bool initCuts)
    {
        mesh=_mesh;
        AddAttributesIfNeeded();
        if (orient_globally)
            vcg::tri::CrossField<MeshType>::MakeDirectionFaceCoherent(*mesh);
        if (initMM)
            InitMMatch();
        SelectSingularityByMM();
        if (initCuts)
        {
            InitTopologycalCuts();
            AddSeamsByMM();
        }
    }

    SeamsInitializer(){mesh=NULL;}
};

#endif
