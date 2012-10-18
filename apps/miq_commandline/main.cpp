#include "mesh_type.h"
#include <wrap/miq/MIQ.h>
#include <wrap/miq/quadrangulator.h>
#include <wrap/io_trimesh/import_ply.h>
#include <wrap/io_trimesh/import_off.h>
#include <wrap/io_trimesh/import_obj.h>

using namespace std;
CMesh trimesh;
MyPolyMesh polymesh;

bool OpenTriMesh(std::string PathMesh)
{
    int position;
    position=PathMesh.find(".ply");
    if (position!=-1)
    {
        int err=vcg::tri::io::ImporterPLY<CMesh>::Open(trimesh,PathMesh.c_str());
        return (err==vcg::ply::E_NOERROR);
    }
    position=PathMesh.find(".obj");
    if (position!=-1)
    {
        int mask;
        bool readed=vcg::tri::io::ImporterOBJ<CMesh>::LoadMask(PathMesh.c_str(),mask);
        if (!readed)return false;
        int err=vcg::tri::io::ImporterOBJ<CMesh>::Open(trimesh,PathMesh.c_str(),mask);
        return (err==vcg::tri::io::ImporterOBJ<CMesh>::E_NOERROR);
    }
    position=PathMesh.find(".off");
    assert(position!=-1);
    int mask;
    bool readed=vcg::tri::io::ImporterOFF<CMesh>::LoadMask(PathMesh.c_str(),mask);
    if (!readed)return false;
    int err=vcg::tri::io::ImporterOFF<CMesh>::Open(trimesh,PathMesh.c_str(),mask);
    return (err==0);
}

int main(int argc, const char * argv[])
{
    const char* configfile;

    if (argc != 2)
    {
        cout << "Not enough parameters: ./MIQ configfile" << endl;
        exit(EXIT_FAILURE);
    }

    configfile = argv[1];
    printf("configuration file %s",configfile);
    fflush(stdout);

    // Read Config File
    FILE *f= fopen(configfile,"r");
    if (f==NULL)
    {
        printf("Cannot open config file\n");
        return -1;
    }
    char buff[200];

    // Mesh name
    fscanf(f,"mesh=%s\n",&buff[0]);
    std::string filename = std::string(buff);
    printf("FILENAME %s",filename.c_str());
    fflush(stdout);
    bool opened=OpenTriMesh(filename);
    if (!opened)
    {
        printf("error loading mesh file \n");
        exit(0);
    }

    vcg::tri::UpdateBounding<CMesh>::Box(trimesh);
    vcg::tri::UpdateNormal<CMesh>::PerVertexNormalizedPerFace(trimesh);
    vcg::tri::UpdateNormal<CMesh>::PerFaceNormalized(trimesh);
    vcg::tri::UpdateTopology<CMesh>::FaceFace(trimesh);
    vcg::tri::UpdateTopology<CMesh>::VertexFace(trimesh);
    vcg::tri::UpdateFlags<CMesh>::FaceBorderFromFF(trimesh);
    vcg::tri::UpdateFlags<CMesh>::VertexBorderFromFace(trimesh);

    // Field name
    fscanf(f,"field=%s\n",&buff[0]);
    std::string fieldname = std::string(buff);
    int position=fieldname.find(".ffield");
    if (position==-1)
    {
        printf("error loading mesh file \n");
        exit(0);
    }
    bool field_loaded=vcg::tri::io::ImporterFIELD<CMesh>::LoadFFIELD(trimesh,fieldname.c_str());
    if (!field_loaded)return false;

    int scalegradient;
    fscanf(f,"scalegradient=%d\n",&scalegradient);

    // Gradient Size
    float GradientSize;
    fscanf(f,"gradient=%f\n",&GradientSize);
    if (scalegradient!=0)
    GradientSize*=1.0/trimesh.bbox.Diag();

    // Stiffness
    float Stiffness=4;

    // DirectRounding
    int DirectRound;
    fscanf(f,"directround=%d\n",&DirectRound);

    // Number of iterations
    int iter=10;

    // Number of local iterations
    int localIter=5;

    // Output name
    fscanf(f,"out=%s\n",&buff[0]);
    std::string out = std::string(buff);
    position=out.find(".obj");
    if (position==-1)
    {
        printf("error output mesh file \n");
        exit(0);
    }

    bool isvalid=MIQ_parametrization<CMesh>::IsValid(trimesh);
    if (!isvalid)
    {
        printf("mesh not valid for parametrization \n");
        exit(0);
    }
    MIQ_parametrization<CMesh>::InitSeamsSing(trimesh,true,true,true);
    MIQ_parametrization<CMesh>::Parametrize(trimesh,MIQ_parametrization<CMesh>::ITERATIVE,Stiffness,GradientSize,(bool)DirectRound,iter,localIter,true);

    Quadrangulator<CMesh,MyPolyMesh> Quad;
    Quad.TestIsProper(trimesh);
    Quad.Quadrangulate(trimesh,polymesh);
    vcg::tri::io::ExporterOBJ<MyPolyMesh>::Save(polymesh,out.c_str(),0);

    fclose(f);
}

